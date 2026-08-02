# LAMBDA & Integer Ambiguity Resolution

## Overview

Two apparently different problems in this codebase reduce to the same mathematics. When a
PPP filter needs to fix carrier-phase ambiguities to integers, it has a float estimate
`a_hat` and its covariance `Q_ahat`, and it wants the most probable integer vector. When a
cycle-slip detector wants to *repair* a slip, it has a float estimate of the cycle jump and
its covariance, and it wants the most probable integer vector. Both are the mixed-integer
least-squares problem `min ||a_hat - a||²_Q` over `a ∈ Z^n`, and both are solved here by the
same subsystem: `com.gnssAug.helper.lambdaNew`, a Java port of the TU Delft **LAMBDA 4.0**
toolbox (Massarweh, Verhagen & Teunissen). The port keeps the toolbox's structure — a single
dispatching entry point, one class per integer estimator, and shared decorrelation /
success-rate machinery — so the MATLAB reference implementation remains a usable cross-check.
The theoretical derivation and the validated comparison of these estimators live in
[thesis/06-lambda-estimators.md](thesis/06-lambda-estimators.md); this document is the
engineering reference for the code as it stands.

Callers today are the Kalman filters in `Rinex/estimation/` (`EKF_PPP`, `EKF_PPP_DF`,
`EKF_PPP_LowCostRx`) and `Android/estimation/KalmanFilter/` (`EKF_PPP`, `EKF_PPP3`,
`EKF_TDCP_ambFix2`, `EKF_TDCP_ambFix_allEst`). See [ppp-engine.md](ppp-engine.md) for how the
PPP filters build `a_hat`/`Q_ahat` in the first place, [kalman-filters.md](kalman-filters.md)
for the filter structure around it, and
[thesis/05-cycle-slip-detection.md](thesis/05-cycle-slip-detection.md) for the cycle-slip
side of the same problem.

## Class inventory (`helper/lambdaNew/`)

| Class | Responsibility |
| --- | --- |
| `LAMBDA.java` | Main entry point: decorrelate, dispatch to one `EstimatorType`, back-transform, return a `LambdaResult`. |
| `LAMBDA_all.java` | Same pipeline, but runs *every* estimator on the same input and returns a `HashMap<EstimatorType, LambdaAllResult>` — used for side-by-side comparison runs. |
| `EstimatorType.java` | Enum of the dispatchable estimators: `ILS`, `IA_FFRT`, `BIE`, `PAR`, `PAR_FFRT`, `ALL`. |
| `Estimators/EstimatorILS.java` | Integer Least-Squares by search-and-shrink; also owns the `ILSResult` type and the candidate-sorting helper shared by its subclass. |
| `Estimators/EstimatorILS_GHAH.java` | Faster ILS search (Ghasemmehdi–Agrell recursions); subclass of `EstimatorILS`, and the implementation actually wired in at every call site. |
| `Estimators/EstimatorIA_FFRT.java` | Integer Aperture with a Fixed Failure-rate Ratio Test — accepts or rejects the ILS fix. |
| `Estimators/EstimatorBIE.java` | Best Integer Equivariant: weighted average over integer candidates inside a search ellipsoid. |
| `Estimators/EstimatorPAR.java` | Partial Ambiguity Resolution: fix the largest high-success-rate subset, condition the rest. |
| `Estimators/EstimatorPAR_FFRT.java` | PAR whose subset-acceptance gate is an FFRT ratio test rather than success rate alone. |
| `Estimators/EstimatorIR.java` | Integer Rounding — componentwise rounding. Ported for completeness; no dispatch path. |
| `Estimators/EstimatorIB.java` | Integer Bootstrapping — sequential conditioned rounding. Ported for completeness; no dispatch path. |
| `DecomposeLtDL.java` | In-place LtDL decomposition `Q = Lᵀ·diag(d)·L` of the ambiguity vc-matrix. |
| `TransformZ.java` | Builds the admissible (unimodular) Z-transformation by integer Gauss reduction and reordering of conditional variances. |
| `DecorrelateVC.java` | Composes the two above into the decorrelation step every estimator runs on; also defines the package-private `DecorrelateVCResult` carrier. |
| `ComputeSR_IBexact.java` | Closed-form Integer-Bootstrapping success rate — full, per-component, and cumulative over subsets. |
| `ComputeFFRTCoefficient.java` | Looks up / interpolates the ratio-test critical value `mu` for a target maximum failure rate (embedded coefficient tables, dimensions up to 66). |
| `ComputeVariance.java` | Monte Carlo estimation of the *fixed* solution's covariance and of empirical success/failure/abstain rates; multi-threaded, one method per estimator and one that does all of them at once. |
| `OptimizedVarCalc.java` | Reduces a matrix of Monte Carlo fixed-solution samples into a symmetric covariance plus success / failure / abstain fractions. |
| `GammaIncompleteInverse.java` | Newton–Raphson inverse of the regularized incomplete gamma function — supplies the chi-squared radius for BIE and PAR. |
| `Utilities.java` | Input validation (symmetry, positive-definiteness, dimensions) and the Chebyshev-based Monte Carlo sample-count formula. |
| `BIEvariance.java` | Diagnostic: cross-checks three independent ways of computing the BIE covariance. Not on any active path. |

## How `LAMBDA.java` dispatches

`LAMBDA.computeLambda(aHat, QaHat, method, estimateVar, varArgs...)` is the single entry
point, and every estimator sees the same pre- and post-processing:

1. **Origin translation.** Each float ambiguity is shifted by its floor so the working vector
   lies in `(-1, 1)`. This is purely numerical; the offset is added back at the end.
2. **Decorrelation.** `DecorrelateVC` runs `DecomposeLtDL` then `TransformZ`, producing the
   decorrelated float vector `zHat`, the transformed `L`/`d` factors, the Z-matrix and its
   inverse-transpose. Every estimator works in Z-space, where the search ellipsoid is close
   to spherical and the search is dramatically cheaper.
3. **Baseline success rate.** `ComputeSR_IBexact` gives the bootstrapped success rate from the
   conditional variances `d`, reported in `LambdaResult` regardless of method.
4. **Dispatch** on `EstimatorType` via a `switch`. Optional `varArgs` are interpreted
   per-method: candidate count and minimum success rate for `ILS`/`PAR`/`PAR_FFRT`, maximum
   failure rate for `IA_FFRT`, and the tail probability `alphaBIE` for `BIE`.
5. **Covariance of the fix.** When `estimateVar` is true the branch calls `ComputeVariance`
   (or reuses the stats the PAR estimators already computed) to get a Monte Carlo covariance
   for the fixed solution; when false, a near-zero diagonal is substituted, which asserts the
   fix is deterministic. An `IA_FFRT` rejection (`nFixed == 0`) instead returns the float
   covariance unchanged.
6. **Back-transformation.** The fixed vector and its covariance are mapped out of Z-space with
   `iZt` and the origin offset is restored. For `PAR`, `BIE` and `PAR_FFRT` the squared norm is
   recomputed as a Mahalanobis distance by back-substitution, because those estimators do not
   produce a directly comparable norm during the search.

The returned `LambdaResult` carries the fixed vector, its covariance, the squared norm, the
number of fixed components, the bootstrapped success rate, and the Z-matrix and decorrelated
covariance for callers that want to inspect the transformation.

`LAMBDA_all.java` mirrors this exactly but skips the `switch`: it runs ILS, PAR, IA-FFRT, BIE
and PAR-FFRT on one input, back-transforms each, and returns them keyed by `EstimatorType`.
Its Monte Carlo pass is shared — `ComputeVariance.computeVarianceAll` draws one sample set and
pushes every sample through all three sampling-based estimators — so comparison runs cost far
less than five separate `LAMBDA` calls. `EKF_TDCP_ambFix_allEst` is the caller that uses it.

## What distinguishes each estimator

- **ILS** (`EstimatorILS`, `EstimatorILS_GHAH`) — the optimal integer estimator. A depth-first
  search descends the LtDL levels choosing conditioned nearest integers, and each complete
  candidate shrinks the ellipsoid radius so the remaining search space contracts. Returns the
  `nCands` best candidates ordered by squared norm. Always fixes all `n` components.
- **IA-FFRT** (`EstimatorIA_FFRT`) — ILS plus an accept/reject decision. It asks for the two
  best ILS candidates and compares their squared-norm ratio against a critical value `mu`
  obtained from `ComputeFFRTCoefficient` for the requested maximum failure rate. If the best
  candidate is not decisively better than the runner-up, the fix is rejected and the float
  solution is returned with `nFixed == 0`. If the bootstrapped failure rate is already below
  the threshold the test is skipped and the ILS fix is accepted outright.
- **BIE** (`EstimatorBIE`) — abandons the requirement that the answer be an integer. It
  enumerates every integer inside a chi-squared ellipsoid (radius from
  `GammaIncompleteInverse` at probability `1 - alphaBIE`) by level-wise recursion and returns
  the Gaussian-weighted average of them. The result is real-valued and minimum-mean-squared-
  error by construction; it degrades smoothly toward the float solution when no integer
  dominates. If the ellipsoid comes up empty it falls back to weighting the nearest
  `1 + 2(2ⁿ - 1)` ILS candidates.
- **PAR** (`EstimatorPAR`) — uses the *cumulative* success rates from `ComputeSR_IBexact`,
  which are ordered most-precise-last after decorrelation, to find the largest trailing subset
  whose joint success rate clears `minSR`. That subset is fixed with ILS (optionally BIE); the
  remaining components are conditioned on the fix rather than dropped. If no subset qualifies,
  the float solution is returned with `nFixed == 0`.
- **PAR-FFRT** (`EstimatorPAR_FFRT`) — same subset-shrinking loop, but a subset is accepted as
  soon as *either* the cumulative success rate clears `minSR` *or* an FFRT ratio test on that
  subset passes. Accepting on the ratio test lets it fix subsets that the success-rate gate
  alone would reject, and its conditioning step propagates the fixed subset's covariance into
  the conditioned block rather than assuming it is exact.
- **IR** (`EstimatorIR`) — rounds each float component independently. Cheapest and weakest;
  ignores correlation entirely.
- **IB** (`EstimatorIB`) — rounds sequentially last-to-first, each component conditioned on the
  ones already fixed via the `L` factor. Strictly better than IR, and its analytical success
  rate (`ComputeSR_IBexact`) is what the whole subsystem uses as a cheap proxy for ILS
  reliability.

`EstimatorIR` and `EstimatorIB` are ported but not reachable from `LAMBDA.java`: they have no
`EstimatorType` constant and no call sites. They are kept because the analytical IB success
rate is central to PAR and IA-FFRT, and because a future estimator (see *Extending this*) would
use them as block estimators.

## `EstimatorILS_GHAH` vs `EstimatorILS`

`EstimatorILS_GHAH` subclasses `EstimatorILS` and overrides `estimatorILS()` only — same
inputs, same `ILSResult`, same candidate ordering. It implements the faster sphere-decoding
recursions of Ghasemmehdi & Agrell (2011), with two changes to the search:

1. **The `k0` shortcut.** When only one candidate is requested, the search returns to level 2
   after storing a candidate instead of level 1. At level 1 the best integer is always
   `round(acond[1])`, so enumerating alternatives there cannot improve the solution.
2. **Path-tracked, column-wise `S` updates.** The base class recomputes a full row of the
   accumulator matrix `S` on every descent. GHAH keeps a `path` array recording the highest
   level each level last descended from, and updates only the single diagonal element it
   actually reads, accumulating column-by-column along that path. Levels skipped during an
   ascent are therefore not re-traversed on the next descent.

Both are exact — GHAH returns the same candidates and squared norms — so the change is purely
one of cost, and it compounds because ILS is called in the innermost loop of IA-FFRT, PAR,
PAR-FFRT, the BIE fallback, and every Monte Carlo sample in `ComputeVariance`. That is why it
is the default: `LAMBDA.java`, `LAMBDA_all.java`, `EstimatorIA_FFRT`, `EstimatorPAR_FFRT`,
`EstimatorBIE` and `ComputeVariance` all instantiate `EstimatorILS_GHAH`. Reverting is a
mechanical substitution of the class name at those call sites, as the class header notes; the
parent is retained precisely so that remains possible.

## Monte Carlo variance and success-rate machinery

`ComputeSR_IBexact` gives a closed-form success rate, but only for bootstrapping, and it says
nothing about the *covariance* of a fixed solution. Once a fix is fed back into a Kalman
filter, that covariance matters: treating a fix as exact when it may be wrong makes the filter
overconfident. `ComputeVariance` estimates it by simulation:

1. Take the (decorrelated) ambiguity covariance, force exact symmetry, Cholesky-factor it, and
   draw `nSamples` zero-mean Gaussian float vectors from it. Zero mean means the *correct*
   integer solution is the zero vector, which makes scoring trivial.
2. Push every sample through the chosen estimator — `1` = ILS, `2` = IA-FFRT, `3` = BIE — with
   the work split across `Runtime.availableProcessors()` threads.
3. Hand the resulting matrix of fixed solutions to `OptimizedVarCalc`, which forms the sample
   second moment (symmetrised), and classifies each sample: an all-integer zero vector is a
   success, an all-integer non-zero vector is a failure, a `NaN` column is an abstention (only
   IA-FFRT abstains). The return is `{covariance, successRate, failureRate, abstainRate}`.

Sample count comes from `GnssDataConfig.nSamplesMC`; passing zero or null instead derives it
from the bootstrapped success rate using the Chebyshev bound in `Utilities.computeNumSamples`.
`computeVarianceAll` is the batched variant that shares one sample set across ILS, IA-FFRT and
BIE for `LAMBDA_all`. `computeVariance2` and `computeVarianceAll2` are the earlier
single-threaded versions of the same two methods, kept as references and not called.

This is the expensive part of the subsystem — it runs a full integer search per sample — so
callers pass `estimateVar = false` when they only need the fix itself.

## Legacy `lambda/` package

`com.gnssAug.helper.lambda` is the earlier, pre-4.0 LAMBDA implementation (Jama-based rather
than EJML). It is not dead, but its live surface is narrow, and the following was verified by
searching call sites rather than assumed:

- **One caller remains.** `Android/estimation/LLS_TDCP_ambFix.java` is the only file outside
  the package that imports it. It constructs `Lambda` with method `6` (`ILS_SSEARCH_RT`,
  ILS with a ratio test) and falls back to method `5` (`PAR`) when nothing gets fixed. That
  class is the least-squares (non-Kalman) TDCP velocity and cycle-slip path, reached from
  `Android/Android.java` and `Android/Android_Static.java`. The Kalman-filter TDCP path
  (`EKF_TDCP_ambFix2`, `EKF_TDCP_ambFix_allEst`) uses `lambdaNew` instead — see
  [android-pipeline.md](android-pipeline.md).
- **`lambdaNew` does not depend on it.** `LAMBDA.java` still carries
  `import com.gnssAug.helper.lambda.Decorrel;`, but the only reference to `Decorrel` in that
  file is inside a commented-out block preserved as a cross-check against `DecorrelateVC`.
  The active decorrelation is entirely `DecorrelateVC` / `TransformZ` / `DecomposeLtDL`. The
  import is vestigial and removable.
- **Not all of the package is reachable even internally.** `Lambda` uses `Decorrel` (which
  uses `Ldldecom`), `Ssearch`, `Isearch`, `Bootstrap` and `Parsearch3`. `Parsearch.java`,
  `Parsearch2.java` and `MyNormalDistribution.java` have no callers anywhere in the tree
  (`Parsearch2` survives only in a commented-out line of `Lambda.java`).

Treat `lambda/` as frozen. New code should call `lambdaNew`; the remaining use in
`LLS_TDCP_ambFix` is the one migration still outstanding, and it is a contained one — the
call sites need `Jama.Matrix` swapped for `SimpleMatrix` and the method constants mapped onto
`EstimatorType`.

### `helper/IntegerLeastSquares.java`

Unrelated to both packages and unused. It is a standalone brute-force ILS: it enumerates every
integer within ±3 cycles of each float ambiguity, scores the full Cartesian product by
Mahalanobis distance, and applies a hardcoded ratio test of 3. It predates the LAMBDA ports,
has no decorrelation step, and scales as `7ⁿ`. Nothing references it. Use `lambdaNew`.

## Extending this

The estimators are not behind a common Java interface — `LAMBDA.java` dispatches through a
`switch` over `EstimatorType` and each estimator is a plain class with its own result type.
The pattern is conventional rather than enforced, so a new estimator should follow it
deliberately. Taking Vectorial Integer Bootstrapping (VIB, method #6 in the toolbox
documentation reproduced in the `LAMBDA.java` header, not implemented here) as the worked
example:

1. **Add the class** as `Estimators/EstimatorVIB.java`. Accept the decorrelated inputs the
   other estimators take — `zHat` (float vector), `LMat` and `dVec` (the LtDL factors) — plus
   any method-specific options. Do *not* re-decorrelate: `LAMBDA.java` has already done it, and
   estimators that need the covariance rebuild it locally as `Lᵀ·diag(d)·L` (see
   `EstimatorPAR_FFRT`). VIB specifically would partition the levels into blocks and apply IR
   or ILS within each, so it would call `EstimatorIR` / `EstimatorILS_GHAH` per block.
2. **Define a result carrier** as a nested class, following `ILSResult` / `PARResult` /
   `IAFFRTResult`. At minimum expose the fixed vector; expose `nFixed` if the estimator can
   partially fix or abstain, and a squared norm if the search produces a comparable one. If it
   computes its own Monte Carlo statistics, expose them as the `Object[] {covariance,
   successRate, failureRate, abstainRate}` that `OptimizedVarCalc` already returns, as the PAR
   estimators do — `LAMBDA.java` reads that array positionally.
3. **Add the enum constant** to `EstimatorType`.
4. **Add the `switch` case** in `LAMBDA.computeLambda`. The case must set `zFix`, and set
   `nFixed`, `sr` and `sqNorm` if they differ from the defaults (`nFixed = n`, `sr` from
   `ComputeSR_IBexact`, `sqNorm = 0`). Honour `estimateVar`: either call
   `ComputeVariance.computeVariance` with a new estimator code, or set the near-zero diagonal
   placeholder. If the estimator can abstain, return `nFixed == 0` and let the existing
   early-return path hand back the float solution.
5. **Parse options** in the `varArgs` block, keyed on the new `EstimatorType`, and give the
   parameter a default at the top of the method alongside `nCands`, `minSR`, `maxFR` and
   `alphaBIE`.
6. **Extend the Monte Carlo path** if the estimator is to be evaluated: add a `case` to the
   `switch` inside `ComputeVariance.computeVariance` (and to `computeVarianceAll` if it should
   appear in comparison runs). Whatever it returns must be classifiable by `OptimizedVarCalc`
   — all-integer for a fix, `NaN` column for an abstention.
7. **Wire it into `LAMBDA_all`** if it belongs in side-by-side comparisons: add it to the
   estimator arrays there and to the `aFixMap` / `nFixedMap` / `qFixMap` population.

If the estimator produces a real-valued (non-integer) answer, as BIE does, note that
`OptimizedVarCalc` counts only all-integer samples toward success or failure, so its empirical
rates will read as zero — that is by design, not a bug, and `LAMBDA_all` handles it as such.

## See also

- [thesis/06-lambda-estimators.md](thesis/06-lambda-estimators.md) — full theoretical
  derivation of these estimators and the validated comparison results between them.
- [thesis/05-cycle-slip-detection.md](thesis/05-cycle-slip-detection.md) — the cycle-slip
  detection that produces the integer problem this subsystem solves.
- [ppp-engine.md](ppp-engine.md) — the PPP filters that call `LAMBDA` for BSD ambiguity fixing.
- [kalman-filters.md](kalman-filters.md) — the filter implementations around those calls.
- [android-pipeline.md](android-pipeline.md) — the Android TDCP paths, including the one that
  still uses legacy `lambda/`.
- [architecture.md](architecture.md) — how this subsystem sits in the wider codebase.
- Massarweh, L., Verhagen, S., Teunissen, P.J.G. (2024), *New LAMBDA toolbox for mixed-integer
  models: Estimation and Evaluation*, and the references listed in the `LAMBDA.java` header.
