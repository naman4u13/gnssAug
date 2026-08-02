# PPP Engine

## Overview

The PPP engine is the codebase's undifferenced, uncombined Precise Point Positioning (UU-PPP) implementation: a family of Extended Kalman Filters that estimate receiver position from pseudorange, carrier phase, and Doppler, using precise orbit and clock products instead of a differential base station. Every member of the family shares one architecture — a **velocity filter** driven by time-differenced carrier phase (TDCP) and Doppler that also performs cycle-slip detection and repair, feeding a **position filter** that carries the full PPP state vector (position, receiver clocks, velocity, clock drifts, residual wet troposphere, per-satellite float ambiguities, per-satellite ionosphere). The variants differ only in which receiver class they target and which ionosphere/clock parameterisation that implies. `com.gnssAug.Rinex.estimation.EKF_PPP` is the reference implementation and carries by far the best Javadoc in the repository; read its class-level and method-level comments alongside this document, and lean on them for method signatures and parameter semantics rather than duplicating them here.

## Classes

| Class | Responsibility |
| --- | --- |
| `com.gnssAug.Rinex.estimation.EKF_PPP` | Reference UU-PPP filter for geodetic receivers on RINEX data; single- or multi-signal, GIM-constrained ionosphere, BSD ambiguity resolution. |
| `com.gnssAug.Rinex.estimation.RinexPPPStateLayout` | Package-private index arithmetic for `EKF_PPP`'s position-filter state vector. |
| `com.gnssAug.Rinex.estimation.EKF_PPP_DF` | Dual-frequency variant: ionosphere estimated purely from the two-frequency observation geometry, with no GIM pseudo-observations. |
| `com.gnssAug.Rinex.estimation.EKF_PPP_LowCostRx` | Low-cost-receiver variant on RINEX input: adds a separate receiver *phase* clock block alongside the code clock block. |
| `com.gnssAug.Android.estimation.KalmanFilter.EKF_PPP` | Android raw-GNSS-log equivalent of `EKF_PPP_LowCostRx`, consuming `Android.models.Satellite`. |
| `com.gnssAug.Android.estimation.KalmanFilter.EKF_PPP3` | Most configurable Android variant: optional single phase clock, single clock drift, time-sliced filter resets, and API-flag-driven cycle-slip handling. |
| `com.gnssAug.Android.estimation.KalmanFilter.PPPStateLayout` | Index arithmetic for `EKF_PPP3`'s state vector, including the phase-clock block. |
| `com.gnssAug.Android.models.CycleSlipDetect` | Per-satellite, per-epoch carrier of TDCP/Doppler/pseudorange delta-ranges plus slip flags and repair results; the hand-off object between the two filters. |
| `com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig` | Builds the state-transition and process-noise matrices (`configTDCP`, `configPPP_IGS`, `configPPP_Android`). |
| `com.gnssAug.Android.estimation.KalmanFilter.EKFParent` | Base class holding the shared diagnostic maps and their accessors. |

## The dual-filter architecture

Each PPP class instantiates two independent `KFconfig` objects and runs them in sequence every epoch.

**`vel_kfObj` — the velocity / TDCP filter.** State vector `[ vel(3) | clkDrift(nSys) | driftRate(nSys) ]`. It is updated in three passes per epoch (documented in detail on `runFilter_vel`):

1. A **Doppler pass** updates velocity and clock drift from Doppler-derived delta-range, optionally after chi-squared + w-test outlier pruning.
2. A **TDCP detection pass** runs the same statistical test on TDCP residuals; satellites the test rejects are flagged as cycle slips on their `CycleSlipDetect` object.
3. A **TDCP update + repair pass** re-forms the design matrix with one extra float state per slipped satellite (initial variance 1e12, column coefficient = wavelength), updates, and — when `repairCS` is on — hands the float slip vector and its covariance to LAMBDA for integer fixing. The temporary slip states are stripped from the state vector after the update regardless of whether the fix succeeded.

**`pos_kfObj` — the position / PPP filter.** State vector as laid out below. It runs after the velocity filter, consumes the slip flags and repair values the velocity filter wrote onto the `CycleSlipDetect` list, and produces the epoch's position estimate.

The two filters are deliberately decoupled: the velocity filter is small, always full-rank, and re-estimated per epoch, which makes it a stable place to run slip detection; the position filter is large and resized every epoch as satellites come and go. Only the `CycleSlipDetect` list crosses between them. Both filters share the predict/update implementation in `KF` and the process-model builders in `KFconfig` — see [kalman-filters.md](kalman-filters.md).

## The position-filter state vector

`RinexPPPStateLayout` (for `EKF_PPP` and, in inline arithmetic form, `EKF_PPP_DF`) encodes:

```
[ pos(3) | clkOff(clkOffNum) | vel(3) | clkDrift(clkDriftNum) | tropo(1) | amb(n) | iono(ionoParamNum) ]
```

`PPPStateLayout` (for `EKF_PPP3`, and matching the inline layout in `EKF_PPP_LowCostRx` and Android `EKF_PPP`) inserts a phase-clock block:

```
[ pos(3) | codeClk(codeClkOffNum) | phaseClk(phaseClkOffNum) | vel(3) | clkDrift(clkDriftNum) | tropo(1) | amb(n) | iono(ionoParamNum) ]
```

Both layout classes are package-private, final, and carry only `final int` boundary fields computed in the constructor (`clkOffStart`, `velStart`, `clkDriftStart`, `tropoIdx`, `coreStateNum`, `ambStart`, `ionoStart`, `totalStateNum`). They exist purely to replace scattered expressions like `6 + clkOffNum + clkDriftNum + 1` with named indices; introducing one where a filter still does the arithmetic inline is a safe, mechanical refactor.

The block sizes are data-driven and change per epoch:

- `clkOffNum` / `codeClkOffNum` equals the number of **observation codes** being processed (e.g. `{"G1C","G5Q","E1C"}` → 3). Index 0 is the base receiver clock; the rest are inter-signal / inter-system biases relative to it.
- `clkDriftNum` equals the number of distinct **constellations** among those codes.
- `n` is the number of satellite-signal observations in the current epoch — one ambiguity per satellite *per observation code*.
- `ionoParamNum` is the number of **unique satellites** (constellation + PRN) in the epoch — one VTEC state per physical satellite, shared across all its signals.

That last distinction is the mechanism by which a multi-frequency run becomes ionosphere-observable: two signals from the same satellite contribute two ambiguities but map into the *same* ionosphere state with frequency-dependent coefficients `±40.3·10¹⁶/f²`, so the frequency-dependence of the ionospheric delay is what separates it from the geometry.

The troposphere is a single state: the *residual* zenith wet delay. The modelled dry delay was already removed by `IGS.filterSat`, which also stored each satellite's wet mapping function on the `Satellite` object for use as the partial derivative.

### Rank deficiency and the implicit datum

Undifferenced uncombined PPP is rank-deficient in its raw form: receiver clock, satellite clock, receiver and satellite hardware biases, ionospheric delays, and ambiguities are not all separately estimable from code and phase alone. The code resolves this with a specific, mostly implicit set of choices, and understanding them is essential before modifying the state vector or the measurement model:

1. **Satellite-side terms are removed upstream, not estimated.** `IGS.SingleFreq` applies the precise clock product (with DCB folded in by the `Clock` parser) and the OSB code and phase biases before the filter ever sees a measurement. Nothing satellite-side remains in the state vector.
2. **The receiver code clock is the datum.** Observation code index 0 defines the base clock. In `get_z_ze_H` the base clock enters *both* pseudorange and carrier phase for every satellite, while entries for code index `j > 0` add an inter-signal bias that enters that signal's pseudorange only. All other receiver clock parameters are differential with respect to that base.
3. **Receiver phase biases are absorbed into the float ambiguities** in the geodetic parameterisation. `KFconfig.configPPP_IGS` states this explicitly: no separate phase-clock states are included. The direct consequence is that the estimated float ambiguities are *not* integers — they are integer ambiguity plus a receiver-dependent bias plus whatever residual satellite phase bias the OSB product left behind.
4. **Integer recovery therefore happens in the between-satellite-differenced (BSD) domain.** Differencing two satellites observed on the same observation code cancels the receiver-dependent part, leaving a quantity that is integer-valued if the satellite phase biases have been correctly removed. This is why the ambiguity-resolution step builds a `Z` matrix of single differences rather than handing the undifferenced ambiguities to LAMBDA directly.
5. **The ionosphere block needs external support in the single-frequency case.** With one signal per satellite, ionosphere is not separable from the ambiguity and clock terms by geometry alone, so `EKF_PPP` and the low-cost variants append `ionoParamNum` **GIM VTEC pseudo-observations** — one direct observation of each satellite's ionosphere state, weighted by `GnssDataConfig.GIM_TECU_variance`. `EKF_PPP_DF` omits this block entirely because dual-frequency data makes the ionosphere estimable from the observations themselves. The pseudo-observations act as a soft datum constraint: strengthen or loosen them by changing that one variance.
6. **Low-cost variants add a phase-clock block instead of absorbing.** `configPPP_Android` appends `phaseClkOffNum` phase-clock states with their own large process noise, on the grounds that consumer hardware exhibits code and phase clock behaviour that a single clock state cannot track. This re-introduces a near-deficiency between the phase clock and the ambiguity block, which is held in check only by the very different process noise assigned to each (phase clock: `1e4·dt`; ambiguities: `1e-16·dt`, i.e. effectively constant between slips). If you change either constant, expect the split between "phase clock" and "ambiguity" to change with it.

The process-noise model also encodes several structural assumptions worth knowing: position/velocity noise is built in ENU and rotated to ECEF; clock-drift index 0 (base oscillator) gets much larger noise than the inter-system drift states (slow thermal variation); the troposphere is a slow random walk; ambiguities get a tiny non-zero noise purely to keep `Q` positive definite; and ionosphere noise is scaled by `1/sin²(elevation)` so that a constant slant-domain noise level is achieved. `KFconfig` asserts positive-definiteness of `Q` and throws if it fails — a useful early warning when a tuning change goes wrong.

### Per-epoch state rebuilding

The single most intricate mechanic in `runFilter_pos` is not the filter maths but the bookkeeping: because `n` and `ionoParamNum` change every epoch, the whole state vector and covariance matrix are **rebuilt from scratch each epoch** before the predict step, with `ambIndexMap` and `ionoIndexMap` (satellite ID → previous state index) carrying continuity.

For each satellite in the current epoch, the ambiguity state is carried forward if the satellite was tracked cleanly (`!isCS()`) or was slipped-but-repaired-and-not-reset; in the carried-forward case the repaired integer slip and its variance are *added* to the state and its variance. Otherwise a fresh ambiguity is initialised from `(phase − pseudorange)/λ` with variance 1e16. Critically, the rebuild also copies the **cross-covariances** — ambiguity↔core, ambiguity↔ambiguity, and ambiguity↔ionosphere — rather than just diagonal terms, which is what preserves the filter's convergence across satellite set changes. Ionosphere states are handled the same way, falling back to zero with variance 1e6 for a newly seen satellite. This nested index-mapping code is verbose and easy to break; if you change the state layout, this is the section to update most carefully.

## Cycle-slip detection and repair

Detection is layered:

1. **Pre-filter geometry-free screen.** While building the `CycleSlipDetect` list, `iterate` compares TDCP against Doppler delta-range; a discrepancy exceeding `csThresh` (100) wavelengths means the observation is beyond repair, so the satellite is excluded from the epoch outright and counted in `resetCount_gfTest`.
2. **Statistical detection in the velocity filter.** `performTesting` computes the global test statistic `T = vᵀ(HPHᵀ+R)⁻¹v` against a chi-squared distribution at α = 0.01, and while the test fails removes the satellite with the largest normalised w-statistic, stopping when fewer than half the satellites remain. On the TDCP pass, removed satellites are flagged `setCS(true)`.
3. **Consecutive-slip hard reset.** `consecutiveCSmap` tracks back-to-back detections; a satellite exceeding `consecutiveSlips` (1) is marked `setReset(true)`, excluded from the TDCP update, and forced to re-initialise its ambiguity in the position filter rather than being repaired.
4. **Receiver-flag-driven detection (Android only).** `EKF_PPP3` additionally consumes the platform's own measurement flags — ADR validity, ADR reset, ADR cycle-slip, half-cycle resolution — and maintains separate counters (`resetCount_invalidADR`, `resetCount_adrReset`, `resetCount_halfCycle`, `csDetected_api`) so the contribution of each detection source can be attributed in the run log.

Repair happens inside the velocity filter's third pass, via LAMBDA in PAR-FFRT mode. On success, the fixed integers and their variances are written back onto each `CycleSlipDetect` object (`setRepaired`, `setIntAmb`, `setIntAmbCov`), which is exactly what the position filter's state rebuild reads. `EKF_PPP` also contains an `injectArtificialCS` switch and `applyArtificialCS` method that deterministically injects growing synthetic slips into the first six satellites — a validation harness for the detection/repair chain, off by default.

The theory and the validated detection/repair performance are covered in [thesis/05-cycle-slip-detection.md](thesis/05-cycle-slip-detection.md).

## Ambiguity resolution

When `fixAmb` is set, the position filter attempts integer resolution after every float update:

1. Satellites are grouped by observation code (`"G1C"`, `"E1C"`, …).
2. Within each group the highest-weight satellite — weight from `Weight.computeCoVariance(CN0, elevation)`, slipped satellites excluded — becomes the reference. If any group has no eligible reference, AR is skipped for the epoch.
3. A `Z` matrix of between-satellite single differences is built, giving `sd_a = Z·â` and `sd_Q = Z·Q_â·Zᵀ` (symmetrised).
4. `LAMBDA.computeLambda(sd_a, sd_Q, EstimatorType.PAR, …)` with a 99.999 % success-rate threshold resolves as many single-differenced integers as it can. The integration surface is deliberately narrow: the filter passes a float vector and its covariance, and receives back `nFixed`, a success rate, the fixed vector, and its covariance via `LambdaResult`. Swapping estimators is a one-line change of the `EstimatorType` argument. See [lambda-ambiguity-resolution.md](lambda-ambiguity-resolution.md) for what happens inside.
5. If anything was fixed, the SD fix is mapped back to the undifferenced ambiguities (`sd_Q` inverted with a small Tikhonov `ε·I` regularisation, since short arcs and poor geometry make it near-singular), and a **conditional mean and covariance update** propagates the fix through the rest of the state. `EKF_PPP` does this properly across three blocks — core states, ambiguities, and ionosphere — recomputing `Pbb`, `Pcc`, `Pbc`, `Pba`, `Pca` so that no cross-covariance is silently dropped, then symmetrises the result before writing it back. The velocity filter's slip repair uses the same conditional-update formula over a two-block (core, slip) partition.

Note that this fix is applied to the filter state itself, not held as a separate "fixed solution" alongside a float solution — subsequent epochs propagate from the fixed state.

## What distinguishes the variants

| Variant | Input | Clock parameterisation | Ionosphere | Distinguishing features |
| --- | --- | --- | --- | --- |
| `Rinex…EKF_PPP` | RINEX, geodetic | Code clocks only (`configPPP_IGS`); phase biases absorbed into ambiguities | Per-satellite VTEC + GIM pseudo-observations | Reference implementation. Uses `RinexPPPStateLayout`; per-observable-group outlier testing via variance inflation; full analysis map set including per-state 1-σ and GIM residual/innovation maps; artificial-slip injection harness. |
| `Rinex…EKF_PPP_DF` | RINEX, dual-frequency | Same as `EKF_PPP` (`configPPP_IGS`) | Per-satellite VTEC estimated from the two-frequency geometry — **no** GIM block | Measurement vector is `3n` rows (code / phase / Doppler) rather than `3n + ionoParamNum`. Otherwise the closest sibling to `EKF_PPP`; still uses inline index arithmetic instead of a layout class. |
| `Rinex…EKF_PPP_LowCostRx` | RINEX, low-cost receiver | Code **and** phase clock blocks (`configPPP_Android`, one phase clock per code) | Per-satellite VTEC + GIM pseudo-observations | Same measurement stacking as `EKF_PPP`; splits its clock diagnostic maps by `Measurement.Pseudorange` / `CarrierPhase` and its drift maps by `Doppler` / `TDCP`, reflecting the separately-tracked code and phase clocks. |
| `Android…EKF_PPP` | Android raw GNSS log | Code and phase clocks (`configPPP_Android`), optional `predictPhaseClock` coupling | Per-satellite VTEC + GIM pseudo-observations | Functionally the Android-input twin of `EKF_PPP_LowCostRx`, consuming `Android.models.Satellite` and its bootstrap from `Android.estimation.LinearLeastSquare`. |
| `Android…EKF_PPP3` | Android raw GNSS log | Configurable: `singlePhaseClock` collapses the phase-clock block to one state; `singleClockDrift` collapses drifts to one state | Per-satellite VTEC + GIM pseudo-observations | The most parameterised variant. Uses `PPPStateLayout`; forces GPS to be constellation index 0; consumes receiver ADR/half-cycle flags for slip detection with per-cause reset counters; and supports `doTimeSlice`, which reinitialises the whole filter at a fixed interval to produce repeated independent convergence runs from one dataset. Filter initialisation is factored into a reusable `initializePPPFilter` so the reset path and the startup path share code. |

The `predictPhaseClock` flag, threaded through several variants, couples the base phase-clock state to the base drift state (`phi[phase][drift] = dt`) so the phase clock follows a Wiener-process model instead of a pure random walk.

## Quality metrics

When `doAnalyze` is on, `performAnalysis` computes per-observable-group post-fit residuals, **redundancy numbers** (trace of the group's diagonal block of `I − HK`, i.e. the effective number of redundant observations in that group), and the **post-variance of unit weight** (weighted squared residuals over redundancy, ideally ≈ 1.0 when the stochastic model is correctly calibrated). It also stores DOP in ENU, and 1-σ values for position, troposphere, per-satellite ambiguity, and per-satellite ionosphere extracted from the post-update covariance. These maps are the primary tool for diagnosing a mis-tuned `GnssDataConfig`: a post-variance far from 1 for one observable group points straight at that group's prior variance.

Position-filter outlier handling differs from the velocity filter's: `performObservableTesting` does **not** remove observations, it inflates their measurement variance to 1e16, keeping the vector dimensions stable, and it uses the effective redundancy number rather than the raw observation count as the chi-squared degrees of freedom so the test stays calibrated as variances are progressively inflated.

## Extending this

**Adding a PPP variant.** The variants are copy-and-modify siblings rather than subclasses, so start by copying the closest existing one. Then: (a) decide the clock parameterisation and add or reuse a `configPPP_*` method on `KFconfig` — the shared `configPPP_core` already takes all receiver-class noise constants as parameters, so a new variant usually needs only a new thin wrapper with different constants; (b) create or reuse a state-layout class rather than writing inline index arithmetic; (c) update `get_z_ze_H` for any new state block, remembering it is called twice per epoch (once pre-update to build `z`/`ze`/`H`, once post-update with `doAnalyze` to compute post-fit quantities); (d) update the per-epoch state-rebuild section in `runFilter_pos` including all cross-covariance copies; (e) add an `EstimatorMode` constant, include it in `isPPPMode()` so `IGS.filterSat` stops pre-correcting the ionosphere, and add a dispatch block in the relevant pipeline entry point.

**Changing the ionosphere handling.** Removing the GIM constraint means dropping the last `ionoParamNum` rows from `z`, `ze`, `H`, and `R` — compare `EKF_PPP` against `EKF_PPP_DF` for the exact diff. Adding a different external constraint (a regional model, a different map product) is the same block with a different value and variance source.

**Changing the ambiguity-resolution strategy.** The reference-satellite selection, the `Z` construction, and the LAMBDA call are contiguous in `runFilter_pos`. Different grouping (per constellation instead of per observation code), a different reference criterion, or a different `EstimatorType` are all local edits. The conditional mean/covariance update below them is generic and should not need changing.

**Retuning.** Measurement priors live in `Rinex.constants.GnssDataConfig` / `Android.constants.GnssDataConfig`; process noise lives in `KFconfig.configPPP_core` and its wrappers; initial variances live in each filter's `process` method. `IGS` echoes the `GnssDataConfig` values into every run header, so keep new constants there rather than inline if you want them recorded.

**Reducing duplication.** The five variants share most of their code, and the layout classes were an early step toward consolidation. `EKF_PPP_DF` and `EKF_PPP_LowCostRx` still use inline index arithmetic where a layout class would do; unifying them is the highest-value cleanup available in this package.

## See also

- [kalman-filters.md](kalman-filters.md) — `KF`, `KFconfig`, `EKFParent`, and the non-PPP filters
- [lambda-ambiguity-resolution.md](lambda-ambiguity-resolution.md) — LAMBDA decorrelation, search, and the ILS/PAR/BIE/FFRT estimators
- [rinex-igs-pipeline.md](rinex-igs-pipeline.md) — how `IGS.IGS` builds and dispatches the satellite map these filters consume
- [android-pipeline.md](android-pipeline.md) — the Android entry points for the Android variants
- [corrections-and-models.md](corrections-and-models.md) — solid Earth tide, troposphere, ionosphere, and wind-up models applied before the filter
- [architecture.md](architecture.md) — package-level overview
- [thesis/07-ppp-engine.md](thesis/07-ppp-engine.md) — UU-PPP theory, state-vector derivation, and validated results
- [thesis/05-cycle-slip-detection.md](thesis/05-cycle-slip-detection.md) — cycle-slip detection and repair theory and validation
- [thesis/06-lambda-estimators.md](thesis/06-lambda-estimators.md) — ambiguity estimator theory and validation
