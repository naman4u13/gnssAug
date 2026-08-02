# RINEX / IGS Pipeline

## Overview

The `com.gnssAug.IGS` package is the entry point for positioning a geodetic (IGS-station-class) receiver from RINEX 3 observation files plus MGEX precise products. `IGS.IGS` owns the whole run: it parses every input file once, walks the observation file epoch by epoch, hands each epoch to `IGS.SingleFreq` for satellite-position and hardware-bias correction, applies quality masks and atmospheric models via `IGS.filterSat`, accumulates the result into an epoch-keyed `TreeMap<Long, ArrayList<Satellite>>`, and then feeds that map to whichever estimator the caller selected through `EstimatorMode`. Everything downstream of parsing — least squares, the pseudorange EKF, and the PPP filters — consumes that same satellite map, which is what makes the estimators interchangeable. The classes described here cover parse/correct/filter/estimate/analyze orchestration and the non-PPP estimators; the PPP filter family itself is documented separately in [ppp-engine.md](ppp-engine.md).

## Classes

| Class | Responsibility |
| --- | --- |
| `com.gnssAug.IGS.IGS` | Top-level run orchestrator: parse all inputs, iterate epochs, dispatch to an estimator, compute error statistics, and emit plots and a text log. |
| `com.gnssAug.IGS.SingleFreq` | Per-epoch pre-processor: turns raw RINEX observables into corrected `Satellite` objects (satellite position/velocity, clock, PCO, wind-up, OSB biases, relativity). |
| `com.gnssAug.IGS.models.IGSOrbit` | Immutable record of one SP3 epoch: GPS timestamp plus satellite ECEF coordinates keyed by constellation character and PRN. |
| `com.gnssAug.IGS.models.IGSClock` | Immutable record of one CLK epoch: GPS timestamp plus satellite clock biases keyed by constellation character and PRN. |
| `com.gnssAug.IGS.models.IGSAntenna` | One ANTEX antenna block: eccentricity (PCO), zenith grid, azimuth-dependent and azimuth-independent PCV tables, plus a validity-window check. |
| `com.gnssAug.Rinex.constants.GnssDataConfig` | Geodetic-receiver stochastic constants: prior variances of unit weight per observable, random-walk spectral densities, GIM variance. |
| `com.gnssAug.Rinex.estimation.LinearLeastSquare` | Single-epoch (W)LS position and velocity solver with iterative outlier detection; also used for filter bootstrapping. |
| `com.gnssAug.Rinex.estimation.EKF` | Pseudorange-only Extended Kalman Filter producing a smoothed code solution. |
| `com.gnssAug.Rinex.models.*` | Plain data carriers for parsed RINEX and product content — see the table further down. |
| `com.gnssAug.Rinex.fileParser.*` | RINEX/product readers that populate the model classes — see [file-formats-and-parsers.md](file-formats-and-parsers.md). |
| `com.gnssAug.Android.Android_Static_Rinex` | Variant of the same pipeline for consumer receivers that log RINEX; reuses the entire `Rinex.*` stack. |

## `IGS.IGS` — the run orchestrator

`IGS.posEstimate` is a single static method that takes file paths and a large set of boolean/enum switches, and produces a run's worth of output. Its class Javadoc lists every parameter; the parts worth understanding structurally are the ordering and the side effects.

**1. Output redirection and the run header.** Before any work happens, `System.out` is redirected to a `.txt` file whose name is built from the RINEX site ID, epoch tag, signal list, and estimator mode. `printRunHeader` then writes a structured block containing the UTC timestamp, the current git commit short hash (with a `(dirty)` marker if the tree has uncommitted changes), the RINEX header's marker/receiver/antenna strings, the ground-truth ARP if SINEX was loaded, all estimator switches, the product file names, and every prior variance from `GnssDataConfig`. Every saved log is therefore self-documenting and traceable to an exact code revision — worth preserving if you refactor this method.

**2. Parsing.** RINEX NAV is parsed first (yielding broadcast ephemerides, Klobuchar coefficients, and time-system corrections), then RINEX OBS (yielding the observation messages, the antenna reference point, receiver PCOs, and the sampling interval). Precise products are loaded conditionally: `Orbit` (SP3) and `Clock` (CLK, with DCBs folded in) when `useIGS` is set, `OSB_Bias` and `DCB_Bias` when `useBias` is set, `IONEX` when `useGIM` is set. If SINEX supplied an ARP but not a receiver PCO for some signal, `IGS` falls back to the ANTEX table via `Antenna.getRxPCO_ENU`, converting the ENU offset to an ECEF delta about the ARP.

**3. Static Orekit setup.** `buildGeoid()` builds a 50×50 EGM96 geoid over WGS-84 in ITRF 2014 / IERS 2010 and caches it in a static field. `ComputeTropoCorr` uses it to get orthometric station height for the pressure/temperature model. Note the two hard-coded absolute paths in this class — the Orekit data directory and the input-file base path — which are the first thing to parameterise if the code is ever run on another machine.

**4. Per-epoch loop.** For every `ObservationMsg`: call `SingleFreq.process` to build the epoch's `Satellite` list, drop the epoch if fewer than `minSat` satellites survive, compute a rough position with `LinearLeastSquare` for use as the next epoch's wind-up geometry reference, attach elevation/azimuth via `ComputeEleAzm`, and call `filterSat`. In least-squares modes the epoch is solved immediately inside the loop; in EKF and PPP modes the epoch is only accumulated into `satMap` and `timeList`, and the estimator runs once over the whole session after the loop.

**5. `filterSat` — masking and atmospheric handling.** This static method mutates the satellite list in place. It removes satellites below the elevation cut-off and below the C/N0 mask, then computes ionospheric delay (from the IONEX GIM when loaded, otherwise from the broadcast Klobuchar model) and the full slant tropospheric delay with its wet mapping function. The important branch is the PPP special case: when `EstimatorMode.isPPPMode()` is true, the ionospheric delay is *stored* on the satellite (`setIonoErr`) but **not** removed from the measurements, because the PPP filter estimates ionosphere as a state and uses the stored value only as a pseudo-observation and as an initial value. In non-PPP modes the iono is subtracted normally. Tropospheric delay is always applied, with the wet mapping function stored on the satellite so the PPP filter can use it as a partial derivative for the residual zenith wet delay. Sign conventions are `PR -= iono + tropo` and `CP += iono - tropo`.

**6. Estimator dispatch.** Modes are grouped by `EstimatorMode`:

- `LLS_CODE` / `WLS_CODE` / `LLS_WLS_BOTH` — per-epoch `LinearLeastSquare` for both position (pseudorange) and velocity (Doppler).
- `RINEX_EKF_CODE` (and, historically, `LLS_WLS_BOTH`) — `Rinex.estimation.EKF` over the full session.
- `IGS_PPP_FLOAT` / `IGS_PPP_AR` — `Rinex.estimation.EKF_PPP`; the only difference between the two is that `IGS_PPP_AR` sets `fixAmb = true`.

**7. Analysis and output.** When `doAnalyze` is set, `IGS` pulls the estimator's per-epoch diagnostic maps (residuals, innovations, post-variance of unit weight, redundancy numbers, DOP, satellite counts) and reshapes them into the per-satellite `SatResidual` structures that `GraphPlotter` expects. It prints ENU error statistics against the ARP, cycle-slip statistics, and a data summary; renders the plot set (`graphENU`, `graphSatRes`, `graphPostUnitW`, `graphDOP`, `graphRedundancyPPP`, `createPPPplots`); and optionally writes a per-epoch JSONL file and CSV exports. In PPP mode the output `.txt` is renamed at the end to embed the session's UTC start–end tag.

## `IGS.SingleFreq` — per-epoch pre-processing

`SingleFreq.process` converts one `ObservationMsg` into the `Satellite` list the rest of the pipeline consumes. It has two modes. In **broadcast mode** (`useIGS = false`) it delegates to `ComputeSatPos` and applies no PCO, wind-up, or bias correction; multi-constellation processing is rejected in this mode. In **IGS mode** it runs the full precise chain per satellite: estimate transmit time from the pseudorange, look up the precise clock offset and drift (the DCB is already folded into the clock product by the `Clock` parser), look up OSB code and phase biases in metres, interpolate SP3 position and velocity with a 10th-order Lagrange polynomial, apply the relativistic clock term and recompute the transmit time, then call `Antenna.getSatPC_windup` for the phase-centre-corrected position and the cumulative carrier-phase wind-up.

Two pieces of state deserve attention. First, `phase_windup_map` is a **static** accumulator keyed by signal + PRN that carries wind-up continuity across epochs — meaning `SingleFreq` is not reentrant and cannot process two sessions concurrently. Second, eclipse handling: when `getSatPC_windup` returns `null` (satellite in umbra or inside the post-eclipse recovery window), the satellite is dropped from the epoch *and* its wind-up entry is removed, so continuity restarts cleanly when it re-emerges rather than resuming from a stale value.

## `IGS.models` — precise-product records

These three classes are thin, immutable-by-convention carriers with no logic beyond field access; they exist so the parsers in `Rinex.fileParser` have a typed record to emit and the interpolators have a typed record to consume.

- `IGSOrbit` holds one SP3 epoch: a GPS timestamp and `HashMap<Character, HashMap<Integer, double[]>>` of ECEF coordinates (constellation character → PRN → XYZ). `Rinex.fileParser.Orbit` builds a time-ordered collection of these and Lagrange-interpolates across them.
- `IGSClock` mirrors that shape for clock biases (`Character → Integer → Double`), consumed by `Rinex.fileParser.Clock`.
- `IGSAntenna` holds one ANTEX antenna block — eccentricity vector (PCO, converted from millimetres), zenith-angle grid, validity window, and both the azimuth-independent (`PCV_NOAZI`) and azimuth-dependent (`PCV_AZI`) phase-centre-variation tables. `checkValidity(double[] time)` selects the right block when an antenna has several dated entries. `Rinex.fileParser.Antenna` owns the lookup and applies the result.

Note the package boundary: these live under `IGS.models` but are used exclusively by `Rinex.fileParser`. That is historical rather than architectural; treat them as product-record types for the shared parser layer.

## `Rinex.constants.GnssDataConfig`

A small class of `public static final` constants defining the geodetic-receiver stochastic model. It is deliberately parallel to `Android.constants.GnssDataConfig`, which carries looser values for smartphone hardware — the geodetic pseudorange and phase priors are considerably tighter. The fields cover prior variances of unit weight for pseudorange, Doppler, TDCP, and carrier phase; the GIM VTEC variance used to weight the ionospheric pseudo-observations; and the ENU random-walk spectral densities for position and velocity. Changing the filter's stochastic tuning almost always means editing this class rather than the filters, and every value is echoed into the run header so a log always records the tuning that produced it.

## Non-PPP estimators

### `Rinex.estimation.LinearLeastSquare`

A static-method solver for single-epoch positioning. `getEstPos` (pseudorange → position + clock) and `getEstVel` (Doppler → velocity + drift) are thin wrappers over `process`, which builds the weight matrix (identity for LS, elevation/C/N0-based via `Weight` for WLS), calls `estimate` to iterate the linearised normal equations, and optionally runs detection-identification-adaptation. Two DIA strategies exist behind an internal `DIA_type` switch; the active one removes the single worst satellite at a time via `qualityControl2` and re-estimates until no satellite fails, marking rejected satellites with `setOutlier(true)`.

Its results (residuals, post-variance of unit weight, state covariance, DOP, the surviving satellite list) are stashed in **static** maps keyed by `Measurement` and read back through the `getResidual` / `getPostVarOfUnitW` / `getCxx_hat` / `getDop` / `getTestedSatList` accessors. This makes the class effectively single-threaded and order-dependent: you must read the accessors immediately after the corresponding solve. Beyond its role as a standalone estimator, `LinearLeastSquare` is the bootstrap for every filter in the codebase — the PPP filters call `getEstPos`/`getEstVel` on the first epoch to initialise their state vectors, and `IGS` calls it every epoch to get a reference position for wind-up geometry.

### `Rinex.estimation.EKF`

A pseudorange-only Extended Kalman Filter with state vector `[pos(3) | clkOff(m) | clkDrift(m)]` for `m` observation codes. Position is initialised from a first-epoch WLS solution; clock states start with large variance. The filter predicts and updates once per epoch starting at the second epoch, checks the covariance for positive semi-definiteness, and (when analysing) rotates the position covariance into ENU before storing it. It shares the `KF` / `KFconfig` machinery described in [kalman-filters.md](kalman-filters.md) and exposes the same diagnostic-map accessors as the other estimators, so `IGS` can treat it uniformly. Use it when you want a smoothed code-only solution with no carrier-phase handling; it is much simpler than the PPP filters and a good place to start reading the filter code.

## `Rinex.models` — what flows through the pipeline

| Class | Content |
| --- | --- |
| `Observable` | One signal's raw measurements for one satellite: pseudorange, carrier phase and cycles, Doppler / pseudorange rate, C/N0, lock and LLI flags, plus the derived carrier frequency and wavelength for its observation code. |
| `Satellite` | `Observable` plus everything computed for it: satellite ECEF/ECI position and velocity, clock offset and drift, reception time, elevation/azimuth, iono and tropo errors, wet mapping function, and outlier flags. This is the object every estimator consumes. |
| `ObservationMsg` | One RINEX epoch: timestamp, GPS week, satellite count, and a nested map of observables keyed by constellation → PRN → band → signal. |
| `NavigationMsg` | One broadcast ephemeris record (Keplerian elements, clock polynomial, correction terms, TOE/IODE). |
| `IonoCoeff` | The eight Klobuchar alpha/beta coefficients from the RINEX NAV header. |
| `TimeCorrection` | GPS-UTC time-system correction (`a0`, `a1`, reference epoch, week) and leap seconds from the NAV header. |
| `SatResidual` | One residual or innovation sample tagged with time, elevation, C/N0, and outlier flag — the unit of data `GraphPlotter` renders. |

The parser classes in `Rinex.fileParser` that produce these are covered in depth in [file-formats-and-parsers.md](file-formats-and-parsers.md); briefly: `ObservationRNX` reads RINEX 3 OBS files, `NavigationRNX` reads mixed NAV files, `Orbit` reads and interpolates SP3, `Clock` reads CLK (folding in DCBs), `Antenna` reads the ANTEX/CSV antenna tables and computes satellite PCO and phase wind-up, `OSB_Bias` and `DCB_Bias` read SINEX-BIAS observable-specific and differential code biases, `IONEX` reads global ionosphere maps and interpolates slant delay, and `SINEX` reads station solutions for the precise ARP.

## Reuse by `Android.Android_Static_Rinex`

Some consumer devices and low-cost receivers can log RINEX rather than platform-specific raw measurements. `com.gnssAug.Android.Android_Static_Rinex` exists for exactly that case: it is structurally a copy of `IGS.posEstimate` — same parse order, same `SingleFreq`-equivalent pre-processing, same `filterSat` masking and atmospheric handling, same `Rinex.models.Satellite` objects — but with the ground truth supplied as an explicit ECEF array instead of read from SINEX, and with a different set of estimator modes: `RINEX_EKF_CODE` (`Rinex.estimation.EKF`), `RINEX_PPP_LOWCOST` (`EKF_PPP_LowCostRx`), and `RINEX_PPP_DF` (`EKF_PPP_DF`). The practical consequence is that the `Rinex.*` stack is receiver-agnostic: the only thing that changes between an IGS station and a low-cost receiver is which estimator variant you select and which stochastic constants it uses. The Android-native entry points that consume raw GNSS logs rather than RINEX — `Android.Android`, `Android.Android_Static` — are described in [android-pipeline.md](android-pipeline.md).

## Extending this

**Adding an estimator mode.** Add a constant to `Android.constants.EstimatorMode` with an unused legacy code, update the relevant predicate (`isPPPMode()`, `isLLSMode()`, …) so `filterSat` and the plotting code classify it correctly, and add a dispatch block in `IGS.posEstimate` after the existing ones. Estimators are expected to accept the epoch-keyed `TreeMap<Long, ArrayList<Satellite>>` and expose the `EKFParent` diagnostic-map accessors; following that contract means the analysis and plotting sections need no changes.

**Adding a precise-product source.** Add a parser under `Rinex.fileParser` that emits per-epoch record objects in the shape of `IGSOrbit` / `IGSClock` (timestamp plus constellation → PRN → value), give it a `findPts`-style windowing method plus an interpolating accessor if the product is time-sampled, load it conditionally in `IGS.posEstimate` behind a new boolean flag, and consume it in `SingleFreq.process`. Add its filename to `printRunHeader` so runs stay traceable.

**Adding a correction model.** Corrections that apply per satellite per epoch belong in `filterSat` (see [corrections-and-models.md](corrections-and-models.md)); corrections that need transmit-time geometry or product lookups belong in `SingleFreq.process`. If a correction should be estimated rather than applied in PPP mode, follow the ionosphere pattern: store the modelled value on the `Satellite` and skip the subtraction when `estimatorMode.isPPPMode()`.

**Retargeting paths.** The absolute Orekit data path in `buildGeoid()` and the `base_path` / output-path strings in `posEstimate` are the only machine-specific items in the pipeline; hoisting them into a config object is the minimal change needed to run elsewhere.

## See also

- [architecture.md](architecture.md) — how the packages fit together overall
- [ppp-engine.md](ppp-engine.md) — the PPP filter family this pipeline dispatches to
- [kalman-filters.md](kalman-filters.md) — the shared `KF` / `KFconfig` / `EKFParent` machinery
- [file-formats-and-parsers.md](file-formats-and-parsers.md) — RINEX, SP3, CLK, ANTEX, IONEX, SINEX, and bias parsers
- [corrections-and-models.md](corrections-and-models.md) — iono, tropo, solid Earth tide, wind-up, satellite position
- [android-pipeline.md](android-pipeline.md) — the consumer-device entry points
- [utilities.md](utilities.md) — `GraphPlotter`, `Analyzer`, `MakeCSV`, matrix and time helpers
- [thesis/02-architecture.md](thesis/02-architecture.md) and [thesis/09-file-formats-and-parsers.md](thesis/09-file-formats-and-parsers.md) — background and validation
