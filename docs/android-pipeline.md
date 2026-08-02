# Android Pipeline

The `com.gnssAug.Android` package processes raw GNSS measurements logged by an Android
smartphone (via the `GnssLogger` app) into position and velocity solutions. It owns the
three top-level driver classes that `MainApp` dispatches to, the per-epoch satellite
assembly step that turns raw logs plus orbit/clock products into `Satellite` objects, the
tuning constants that every estimator reads, the snapshot (non-Kalman) least-squares
estimators, the Android-specific file parsers and data models, and the IMU/magnetometer
support classes used by the GNSS/INS fusion mode. The recursive filters live in a
sub-package and are documented separately in [kalman-filters.md](kalman-filters.md);
the PPP filters are documented in [ppp-engine.md](ppp-engine.md).

## Classes

| Class | Responsibility |
| --- | --- |
| `Android` | Kinematic driver: runs an Android log against a time-series ground truth track and every non-RINEX estimator mode. |
| `Android_Static` | Static driver: same pipeline against a single fixed reference ECEF, plus the ambiguity-estimator comparison and PPP-AR modes. |
| `Android_Static_Rinex` | Driver for Android data already converted to RINEX; routes through the RINEX/IGS parsers and estimators instead of the Android ones. |
| `SingleFreq` | Per-epoch satellite assembly: pairs raw observations with broadcast (derived CSV) or precise (IGS) orbits/clocks and emits corrected `Satellite` objects. |
| `constants.EstimatorMode` | Enum naming every supported estimator configuration; carries the legacy integer code and mode-family predicates. |
| `constants.GnssDataConfig` | Prior variances of unit weight per measurement type and ENU process-noise (random-walk) vectors. |
| `constants.Constellation` | Carrier-frequency lookup keyed by constellation letter and band number. |
| `constants.Measurement` | Enum of measurement types (`Pseudorange`, `Doppler`, `CarrierPhase`, `TDCP`, `GIM_Iono`) used to key analysis maps. |
| `constants.State` | Enum distinguishing `Position` from `Velocity` state blocks in analysis output. |
| `constants.AndroidSensor` | Enum of logged sensor streams (`GNSS`, `Accelerometer`, `Gyroscope`, `Magnetometer`). |
| `constants.ClockAllanVar` | Allan-variance coefficients converted into clock phase/frequency PSDs for a low-quality TCXO. |
| `constants.ImuDataSheets` | Per-device IMU noise parameters (bias, random walk, correlation time) taken from vendor datasheets. |
| `estimation.LinearLeastSquare` | Snapshot LS/WLS position from pseudorange and velocity from Doppler, with statistical testing and quality control. |
| `estimation.LLS_TDCP` | Snapshot velocity from time-differenced carrier phase between two consecutive epochs. |
| `estimation.LLS_TDCP_ambFix` | Same as `LLS_TDCP` but detects cycle slips and repairs them by integer ambiguity resolution. |
| `fileParser.GNSS_Log` | Parses a `gnss_log.txt` into per-epoch raw GNSS records and a time-ordered IMU list. |
| `fileParser.DerivedCSV` | Parses Google's `device_gnss.csv` / derived CSV (broadcast satellite states and correction terms). |
| `fileParser.GroundTruth` | Parses the Google Decimeter Challenge ground-truth CSV (lat/lon/alt). |
| `fileParser.GroundTruth_GSA` | Parses whitespace-delimited GSA-style ground truth (already in ECEF). |
| `models.GNSSLog` | One raw Android GNSS measurement row, with derived receive time, observation code and week number. |
| `models.Satellite` | A `GNSSLog` enriched with satellite position/velocity, corrected observables, geometry and per-satellite flags. |
| `models.Derived` | One row of the derived CSV: satellite ECEF/velocity, clock bias/drift, raw pseudorange and correction terms. |
| `models.TDCP` | A single-satellite time-differenced carrier-phase observation with its unit line-of-sight and satellite-motion correction. |
| `models.CycleSlipDetect` | Per-satellite record used by the cycle-slip detection and repair path; carries both delta-range flavours and the resolved ambiguity. |
| `models.IMUsensor` | One accelerometer / gyroscope / magnetometer sample, with bias and a GPS-aligned timestamp. |
| `helper.INS.IMUconfigure` | Downsamples and interpolates the raw IMU list onto a uniform grid aligned with the GNSS epochs. |
| `helper.INS.StateInitialization` | Derives the initial body-to-ENU attitude from averaged accelerometer and magnetometer samples. |
| `helper.GeomagneticField` | WMM-2020 implementation used to convert magnetic north to true north. |

## Entry points

All three drivers expose a single static `posEstimate(...)` method that takes file paths,
masking thresholds, an `EstimatorMode`, and a set of analysis flags. They share the same
overall shape — parse inputs, loop over epochs building a `TreeMap<Long, ArrayList<Satellite>>`
keyed by receive time in milliseconds, then hand that map to whichever estimator the mode
selects, and finally produce error statistics and plots. They differ in where truth comes
from and which estimators they can reach.

### `Android`

The kinematic driver. Ground truth is a *time series*, read either from the Google
Decimeter Challenge CSV via `GroundTruth` or from a GSA-format file via `GroundTruth_GSA`
(the `isGSA` flag selects between them), and epochs are advanced in lockstep with the truth
track — an epoch with no matching truth timestamp is skipped. Because truth varies per
epoch, this driver is the one that can supply a per-epoch true-position list to the TDCP
filters, and it is the only one that runs the GNSS/INS fusion mode.

Estimator coverage: LS/WLS snapshot modes, the basic EKF family, the DBP filters
(`EKFDoppler` / `AKFDoppler`), the TDCP snapshot and filtered modes, `PPP_FLOAT` via
`EKF_PPP`, and `ALL_ANALYSIS` which fans out to several sub-modes in one run.

### `Android_Static`

The static driver. Truth is a single `double[] trueEcef` passed in by the caller, so the
epoch loop has no ground-truth cursor to keep in sync and every logged epoch is processed.
It prints a reproducibility header (`printRunHeader`, including the current git hash via
`getGitHash`) at the top of its output file. Two estimator paths exist only here:

- `EKF_TDCP_ALL_ESTIMATORS`, which runs `EKF_TDCP_ambFix_allEst` to compare LAMBDA
  estimator variants side by side.
- `PPP_FLOAT`, which in this driver constructs `EKF_PPP3` rather than the `EKF_PPP` used by
  `Android` (see [ppp-engine.md](ppp-engine.md)).

It also substitutes `AKFDoppler_Static` for `AKFDoppler` when the mode is `AKF_DBP`, since
the static case has no meaningful Doppler-driven state propagation.

### `Android_Static_Rinex`

A bridge rather than a variant. When Android data has already been converted to RINEX
(newer phones and third-party loggers emit RINEX directly), this driver reads it with
`Rinex.fileParser.ObservationRNX` / `NavigationRNX`, builds `Rinex.models.Satellite` objects
through `IGS.SingleFreq`, and estimates with the `Rinex.estimation` classes. Nothing from
`Android.fileParser`, `Android.models`, or `Android.estimation` is used; only
`EstimatorMode`, `Measurement`, and `State` are shared. Its reachable modes are the
`RINEX_*` family: `RINEX_EKF_CODE`, `RINEX_PPP_LOWCOST`, and `RINEX_PPP_DF`. See
[rinex-igs-pipeline.md](rinex-igs-pipeline.md) for the parsers and estimators it delegates
to.

### `SingleFreq`

Despite the name, this is not a driver — it is the per-epoch satellite assembly step called
from inside the `Android` and `Android_Static` epoch loops. For each observation code and
satellite it resolves the transmission time, obtains satellite position/velocity, and
produces a `Satellite` carrying pseudorange, phase and range-rate that have already been
corrected for satellite clock offset and drift.

It has two branches controlled by the `useIGS` flag:

- **Broadcast branch** (`useIGS == false`): satellite states and the ionosphere/troposphere/
  inter-signal-range-bias corrections all come straight from the derived CSV, keyed by
  receive time with a ±1 ms tolerance to absorb rounding.
- **Precise branch** (`useIGS == true`): satellite position comes from polynomial
  interpolation of a precise orbit (`Orbit.getPV`, 10th order), clocks from `Clock`,
  and the antenna phase-centre offset plus carrier phase wind-up from `Antenna`. Relativistic
  clock correction is applied and the transmission time is recomputed with it. Wind-up is
  accumulated across epochs in a static `phase_windup_map`; when `Antenna` rejects a
  satellite (for example during eclipse), its map entry is removed so wind-up continuity
  restarts cleanly on reacquisition.

> [!WARNING]
> That wind-up map is static, so a single JVM run should process one dataset at a time.

## Constants and configuration

`EstimatorMode` is the enum that selects the estimator. Each constant carries a
`legacyCode` integer (the pre-enum representation, still accepted through
`fromLegacyCode`) and the enum exposes family predicates — `isLLSMode`, `isBasicEKFMode`,
`isDBPMode`, `isTDCPMode`, `isPPPMode`, `isAnalysisMode` — that the drivers branch on
instead of testing individual constants. `getSubModes` expands the composite modes
(`LLS_WLS_BOTH`, `EKF_BASIC_MULTI`, `ALL_ANALYSIS`) into the list of concrete modes to run
in one pass. The full mode-to-class dispatch table is documented in
[architecture.md](architecture.md).

`GnssDataConfig` holds the a-priori variance of unit weight for each measurement type
(pseudorange, Doppler, TDCP, carrier phase, GIM TEC) and the ENU random-walk process noise
vectors `qENU_posRandWalk` / `qENU_velRandWalk`. These are the primary tuning knobs. The
drivers echo these values into the run header so an output file records the tuning it was
produced with.

> [!WARNING]
> The file currently carries commented-out alternatives for kinematic tuning alongside the
> active near-zero static values — check which pair is actually uncommented before
> interpreting a run's behaviour.

The remaining constants classes are lookup data rather than tuning: `Constellation` maps a
constellation letter and band index to a carrier frequency, `ClockAllanVar` converts TCXO
Allan-variance coefficients `h0` / `h_2` into the phase and frequency PSDs (`sf`, `sg`) that
`KFconfig` needs, and `ImuDataSheets` holds per-device IMU noise figures (nested classes
`BMI160`, `Pixel4`, `Pixel4_2`) for INS fusion. `Measurement` and `State` are small enums
used almost exclusively as map keys in the analysis plumbing, so that residuals, redundancy
and covariance can be collected per measurement type and per state block.

## Non-Kalman estimation

These three classes are stateless in the filtering sense — each call solves one epoch (or
one epoch pair) independently. All three follow the same internal shape: build the design
matrix and weight matrix, iterate the linearised solution, optionally run quality control,
and stash residuals, `Cyy`, `Cxx_hat` and the posterior variance of unit weight in static
fields that the caller retrieves through getters.

> [!WARNING]
> That static-getter pattern means these classes are not thread-safe — results must be read
> before the next call, or a concurrent call will overwrite them first.

### `LinearLeastSquare`

The workhorse snapshot estimator. `getEstPos` solves for position and one receiver clock
offset per observation code from pseudoranges; `getEstVel` solves for velocity and clock
drift from Doppler, given a reference position. The `isWLS` flag switches between equal
weighting and elevation/C-N0-based weighting supplied by `utility.Weight`. Optional
`doTest` enables statistical testing and `outlierAnalyze` drives blunder detection through
`qualityControl`, which iteratively removes the largest normalised residual.

Its role goes well beyond the LS estimator modes: nearly every Kalman filter in the codebase
calls it to bootstrap. `EKF`, `EKFDoppler`, `AKFDoppler` and the TDCP filters all initialise
their state from a first-epoch WLS solution, and `EKFDoppler` / `AKFDoppler` call
`getEstVel` every epoch as their state-propagation input. Changing the LS behaviour
therefore changes filter behaviour too.

### `LLS_TDCP`

Estimates velocity from the difference of carrier-phase observations between two
consecutive epochs. It matches satellites common to both epochs by observation code plus
SVID, forms a `TDCP` record per common satellite (delta range, unit line-of-sight, and a
correction for satellite motion between epochs), and solves for the three velocity
components plus a clock-drift term. Because it differences phase rather than using it
absolutely, the ambiguity cancels — provided no cycle slip occurred. Reached through
`EstimatorMode.LLS_TDCP_VEL`.

### `LLS_TDCP_ambFix`

Extends the `LLS_TDCP` idea with cycle-slip detection and repair, and is reached through
`EstimatorMode.TDCP_CSDR`. Instead of `TDCP` records it builds `CycleSlipDetect` records,
which hold both the Doppler-derived and phase-derived delta ranges for the same satellite;
their disagreement is the cycle-slip observable. Detected slips get an integer ambiguity
parameter that is resolved with `helper.lambda.Lambda` and folded back into the solution.
The class tracks `ambDetectedCount` / `ambRepairedCount` (both totals and per-epoch maps) so
that detection and repair rates can be plotted.

> [!NOTE]
> This class uses the older `helper.lambda` package, while the filtered equivalent
> `EKF_TDCP_ambFix2` uses the newer `helper.lambdaNew`. See
> [lambda-ambiguity-resolution.md](lambda-ambiguity-resolution.md) for the difference.

## File parsers and data models

Covered in depth in [file-formats-and-parsers.md](file-formats-and-parsers.md); in brief:

- `GNSS_Log` reads the `#`-delimited Android `gnss_log.txt`, splitting `Raw` rows into a
  `TreeMap` keyed by receive-time-in-milliseconds then by observation code, and `UncalAccel`
  / `UncalGyro` / `UncalMag` rows into a time-sorted `IMUsensor` list. It also aligns IMU
  timestamps to GPS time using the logged boot-time reference. Both results are exposed
  through static getters.
- `DerivedCSV` reads the derived/`device_gnss` CSV into a nested map (receive time →
  observation code → SVID → `Derived`), translating Google's signal-type strings into
  RINEX-style observation codes.
- `GroundTruth` and `GroundTruth_GSA` each return an `ArrayList<double[]>` of truth records;
  the former in geodetic coordinates from a CSV, the latter already in ECEF from a
  whitespace-delimited file.
- `GNSSLog` models a single raw measurement row and derives the observation code, receive
  time and week number from it. `Satellite` extends `GNSSLog` and is the object every
  estimator actually consumes — it adds satellite ECEF/velocity, corrected pseudorange,
  phase and range rate, elevation/azimuth, ionospheric and tropospheric error terms, and
  per-satellite flags such as `cycleSlipDetected` and `isValidPhase`.
- `Derived`, `TDCP`, `CycleSlipDetect` and `IMUsensor` are plain carriers for, respectively,
  one derived-CSV row, one time-differenced phase observation, one cycle-slip candidate with
  its float and fixed ambiguity, and one inertial/magnetic sample.

## INS and sensor fusion support

Three classes prepare the inputs that `INSfusion` (see [kalman-filters.md](kalman-filters.md))
consumes. They are only exercised by `EstimatorMode.GNSS_INS`, which is reachable from
`Android` only.

`IMUconfigure.configure` takes the flat, time-sorted `IMUsensor` list from `GNSS_Log` and
produces a `TreeMap<Long, HashMap<AndroidSensor, IMUsensor>>` on a uniform grid at the
requested sample rate, phase-aligned to the first GNSS epoch. Android sensors are not
sampled synchronously, so each stream is downsampled independently and linearly
interpolated (`utility.Interpolator`) onto the common grid; the first and last grid points
are dropped because they may be extrapolations. The driver asserts afterwards that every
grid point carries all three sensor types.

`StateInitialization.initialize` produces the initial body-to-ENU direction cosine matrix.
It averages accelerometer and magnetometer samples over the first 15 grid points — assuming
the device is stationary during that window — and builds an orthonormal rotation from the
gravity and magnetic vectors (`getRotationMatrix`, which throws if the device is in
near-free-fall or too close to the magnetic pole). Because that rotation is referenced to
*magnetic* north, it then queries `GeomagneticField` for the declination at the initial
position and time and rotates the DCM into geographic north.

`GeomagneticField` is a self-contained WMM-2020 implementation (spherical-harmonic
expansion with a `LegendreTable` helper) providing declination, inclination, and field
strength. It exists solely so the initial attitude can be referenced to true north without
an external dependency.

## Extending this

**Adding a new estimator mode.** Add a constant to `EstimatorMode` with an unused
`legacyCode` (the 1–99 band is Android, 100+ RINEX, 200+ IGS), slot it into the appropriate
`isXxxMode()` predicate if it belongs to an existing family, then add a branch in the
relevant driver's epoch-loop that constructs your estimator and writes into `estPosMap` /
`estVelMap` under a display key. Reusing an existing predicate is preferable to adding an
equality test, because the drivers' analysis blocks already key off the predicates.

**Adding a new entry-point scenario.** If the new scenario differs only in configuration
(new dataset, different signals, different truth file), add a `RUN` constant and a `case`
block in `MainApp` rather than a new driver class. Write a new driver only when the *shape*
of the epoch loop changes — a different truth cadence, a different observation source, or a
different `Satellite` model — which is precisely what separates the three existing drivers.

> [!NOTE]
> Output paths are currently hard-coded inside each driver's `posEstimate`, so a new scenario
> means editing that string too.

**Adding a new tuning parameter.** Put it in `GnssDataConfig` (measurement or process noise)
or `ImuDataSheets` (a new device), and echo it in `Android_Static.printRunHeader` so runs
stay self-documenting.

**Extending `SingleFreq`.** New correction terms belong in the `useIGS` branch alongside
wind-up and PCO, and should be applied to `corrPR` / `corrPhase` / `corrPRrate` before the
`Satellite` is constructed, so that every downstream estimator sees them uniformly.

## See also

- [architecture.md](architecture.md) — system overview and the `EstimatorMode` dispatch table
- [kalman-filters.md](kalman-filters.md) — the recursive filters these drivers invoke
- [ppp-engine.md](ppp-engine.md) — `EKF_PPP`, `EKF_PPP3` and the PPP state layout
- [rinex-igs-pipeline.md](rinex-igs-pipeline.md) — the stack `Android_Static_Rinex` delegates to
- [file-formats-and-parsers.md](file-formats-and-parsers.md) — parser and data-model details
- [corrections-and-models.md](corrections-and-models.md) — ionosphere, troposphere, tides, antenna models
- [lambda-ambiguity-resolution.md](lambda-ambiguity-resolution.md) — the LAMBDA packages used by `LLS_TDCP_ambFix`
- [utilities.md](utilities.md) — `Weight`, `Matrix`, `LatLonUtil`, `Analyzer`, `GraphPlotter`
- [thesis/03-pillar1-dbp.md](thesis/03-pillar1-dbp.md), [thesis/05-cycle-slip-detection.md](thesis/05-cycle-slip-detection.md) — theory and validation background
