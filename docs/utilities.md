# Utilities

`com.gnssAug.utility` is the codebase's shared foundation: the linear algebra helpers that sit on
top of EJML, the geodesy and reference-frame conversions, GPS time handling, measurement weighting
and covariance construction, and the diagnostics tooling — plotting, residual analysis, and CSV/JSONL
export — that turns a filter run into something inspectable. It is a flat, dependency-light package
with no internal hierarchy: almost every class is a bag of static methods, and almost every other
package in the project imports from it. Two classes break that pattern and are large enough to
deserve their own treatment: `GraphPlotter`, which is the project's entire visualisation layer, and
`Analyzer`, which computes the truth-referenced error series those plots consume.

There is no state shared between these classes and no initialisation order to respect. The one
external assumption is EJML's `SimpleMatrix`, which is the codebase's standard dense matrix type and
appears in most of the signatures here.

## Class index

| Class | Responsibility |
|---|---|
| `Matrix` | Projection/pseudo-inverse operators, covariance diagnostics, EJML↔array bridging |
| `Vector` | 3-vector arithmetic on raw `double[]` — cross, dot, norm, normalise |
| `MathUtil` | Euclidean distance, factorials/combinations, RMS and MAE aggregates |
| `LatLonUtil` | WGS-84 constants and all ECEF/geodetic/ENU/NED conversions and Earth radii |
| `Time` | GPS week-and-seconds ↔ `Calendar`/`ZonedDateTime`/Unix conversions |
| `Weight` | Elevation- and C/N₀-dependent measurement variance and covariance matrices |
| `SatUtil` | Line-of-sight unit vectors and satellite-list bookkeeping for estimators |
| `Interpolator` | Lagrange (with analytic derivative) and linear interpolation |
| `Closest` | Nearest-value index lookup, binary-search and linear variants |
| `Combination` | Enumerates all size-*r* combinations of an integer array |
| `ArrayIndexComparator` | Sorts an index array by the descending values of a backing `double[]` |
| `StringUtil` | Fixed-width column splitting and string↔array parsing |
| `SymbolToken` | Splits Fortran-formatted numeric fields that ran together |
| `GraphPlotter` | The project's JFreeChart layer — every diagnostic plot it produces |
| `Analyzer` | Builds truth-referenced error and residual series to feed the plots |
| `MakeCSV` | Structured CSV/JSONL export of filter states, residuals, and quality metrics |
| `Trajectory` | Writes a wide multi-estimator position/velocity trajectory CSV |

## Linear algebra — `Matrix` and `Vector`

`Vector` is deliberately primitive: static operations on raw `double[]`, mostly fixed at length 3,
used wherever allocating a matrix object would be wasteful — line-of-sight geometry, satellite body
frames, tide displacement vectors. `scale` mutates its argument in place while `add`/`subtract`
allocate, which is an inconsistency worth remembering before you pass a shared array to it.

`Matrix` is the more interesting of the two. Alongside array/`SimpleMatrix` bridging
(`matrix2Array`, `matrix2ArrayVec`, `ArrayList2Vector`) and small fixed-size 3×3 helpers, it provides
the least-squares operators the estimators are actually written in terms of:

- `getPseudoInv(A, W)` — the weighted pseudo-inverse `(AᵀWA)⁻¹AᵀW`, the least-squares solution
  operator.
- `getProjection(A, W)` — `A(AᵀWA)⁻¹AᵀW`, the weighted projector onto the column space of the design
  matrix.
- `getPerpendicularProjection(A, W)` — `I − P_A`, the projector onto the residual space. Post-fit
  residuals, redundancy numbers, and reliability measures are all expressed through this operator
  rather than computed ad hoc.
- `getNorm(A, B)` — the quadratic form `AᵀB⁻¹A`, i.e. a squared Mahalanobis distance; used for
  statistical testing against a chi-square threshold.
- `getSkewSymMat(a)` — the cross-product matrix, used in the INS mechanisation.

Two covariance diagnostics round it out. `computeCorrelationMatrix` normalises a covariance matrix
by its diagonal standard deviations, which is how you see at a glance whether ambiguity states are
badly coupled. `computeDecorrelationNumber` returns the ratio of largest to smallest eigenvalue —
the condition number of the covariance — which is the standard way of quantifying how elongated an
ambiguity search space is before decorrelation, and therefore how much the LAMBDA Z-transformation
is buying you. See [lambda-ambiguity-resolution.md](lambda-ambiguity-resolution.md).

`MathUtil` holds the remaining scalar odds and ends: `getEuclidean` (fixed at 3 components,
regardless of input length — do not use it for general vectors), `getFact`/`getCombCount`, and the
`RMS`/`MAE` aggregates used when summarising an error series.

## Geodesy — `LatLonUtil`

The single most-imported class in the package. It owns the WGS-84 constants (`a`, `b`, `f`, `omega_ie`,
`mu`) and every coordinate transform the codebase uses.

**Conversions.** `ecef2lla` is an iterative solution — three fixed Bowring-style iterations on the
reduced latitude — with explicit longitude quadrant handling. `lla2ecef` is the closed-form inverse.
Both take a boolean controlling degrees versus radians, and the defaults differ between them
(`ecef2lla` returns degrees by default; `lla2ecef` requires you to say). Mixing these up is the most
likely source of a silent order-of-magnitude error when using this class, so read the flag at every
call site.

**Local frames.** `enu2ecef` / `ecef2enu` / `ned2ecef` / `ecef2ned` all take a reference ECEF position
and an `isPos` flag distinguishing a *point* (translate by the reference) from a *vector* (rotate
only). Getting `isPos` wrong on a displacement vector adds an Earth radius to your answer. The NED
variants are thin wrappers that flip axes via `enu_ned_convert` — the direction cosine matrix is its
own inverse for that swap, which is why one method handles both directions. `getEnu2EcefRotMat` and
`getEcef2EnuRotMat` expose the rotation alone, which is what the error-analysis and plotting code
uses to express position differences in local East/North/Up.

**Distances and radii.** `getHaversineDistance` (spherical, fast) and `getVincentyDistance`
(ellipsoidal, iterative, accurate) cover horizontal distance. `getNormalEarthRadius` and
`getMeridianEarthRadius` give the prime-vertical and meridian radii of curvature, and `getGravity`
implements the Somigliana normal gravity formula with an altitude correction — all three exist for
the INS mechanisation rather than for GNSS positioning.

**Odds and ends.** `gd2gc` converts geodetic to geocentric latitude, which matters for the
ionosphere thin-shell geometry (the IONEX path passes geocentric latitude, not geodetic).
`lonAddDD` adds a longitude increment with wraparound to (−180, 180], used by the IONEX time
interpolation.

One naming trap to be aware of: the private field `e` is the first eccentricity while the public
`e2` is the *second* eccentricity (not `e` squared), despite the name.

## Time — `Time`

GPS time handling is where GNSS codebases usually accumulate bugs, and this class is explicit about
why. Java's `Calendar` and `Date` do not model leap seconds, so `getGPSTime(...)` computes the
offset from the GPS epoch (1980-01-06) using millisecond arithmetic that is *deliberately* not
leap-second aware — the input is already on the GPS timescale, so no leap-second handling is wanted.
The source comments say so; do not "fix" this by introducing a UTC conversion.

The genuinely leap-second-aware path is separate: `convertToZonedDateTime(week, secondsOfWeek)`
adds the total GPS seconds to the epoch and then subtracts a `LEAP_SECONDS` constant to land on UTC.
That constant is hard-coded and must be updated if a new leap second is ever introduced — it is the
one piece of this class that goes stale with time.

Beyond that, the class is a set of conversions covering the formats the parsers and models actually
produce: `getGPSTime(String[])` for RINEX/SP3 timestamp tokens, `getGPSTimeFromYDOY` for the
year:day-of-year:seconds form used in SINEX-BIAS validity windows, `UnixTimeToGPSTimeMilli` for
Android's millisecond timestamps, and `convertUsingToInstant` to bridge a `Calendar` to
`ZonedDateTime` for the Orekit-based models.

`getDate(gpsTime, weekNo, longitude)` and `convertToUTC(cal, longitude)` are a matched pair
implementing *local* GPS time — GPS time shifted by longitude — which the Klobuchar ionosphere model
needs. As the source comment notes, there is no such thing as a local GPS time in reality; this is a
representation convenience for that one model, and should not leak into general timekeeping.

## Weighting — `Weight`

The stochastic model. `computeCoVariance(CNo, elevation)` is the core formula:

```
var = 10^(−CNo/10) / sin²(elevation)
```

combining a carrier-to-noise-density term with the standard elevation-dependent term, so a
low-elevation weak-signal satellite is downweighted on both counts. Everything else in the class
builds on it.

`computeCovInvMat` / `computeCovInvMat2` produce a diagonal inverse-covariance (weight) matrix over
a satellite list — the duplication exists because `Rinex.models.Satellite` and
`Android.models.Satellite` are separate types with different accessor names, which is a recurring
shape in this codebase. `getNormCyy` / `igs_getNormCyy` are the same pair for covariance rather than
weight: they normalise the weights by their maximum so the largest weight is 1, then scale by an a
priori variance of unit weight to produce a covariance matrix in physical units. That normalisation
step is what makes the a priori variance of unit weight interpretable and comparable across epochs,
which matters for the variance component estimation in the adaptive filters — see
[kalman-filters.md](kalman-filters.md). `normalize(W)` is an alternative RMS-based normalisation of
an existing weight matrix.

An elevation-squared-only alternative is present but commented out in both variants; that is the
switch to flip if you want to isolate the effect of the C/N₀ term.

## Satellite bookkeeping — `SatUtil`

Two jobs. The first is geometry: `getUnitLOS` returns the unit line-of-sight vector from user to
satellite, either for a single pair or as an `n×3` matrix over a satellite list — this is the
geometry block of every design matrix in the project. As with `Weight`, there are parallel
`getUnitLOS` and `igs_getUnitLOS` methods for the two `Satellite` types.

The second is state-vector bookkeeping. `findObsvCodeArray`, `findObsvCodeSet`, and `findSSIset`
extract the ordered set of observation codes or constellations present in a satellite list —
which determines how many receiver clock-bias states the filter needs. `resetVar` handles the
consequence: when a constellation drops out of view and later returns, the state vector's clock
biases no longer line up with the filter's fixed layout. `resetVar` re-expands the state and its
covariance back to the full layout, placing each surviving clock state at its correct index and
filling the reintroduced entries with a large variance (1e10) so the filter treats them as
uninformative rather than trusting a stale value. `createCopy` deep-copies a satellite list, used
where an estimator needs to test a modified list without disturbing the original.

## Numerics — `Interpolator`, `Closest`, `Combination`, `ArrayIndexComparator`

`Interpolator.lagrange(X, Y, x, getDeriv)` is the workhorse behind SP3 orbit interpolation. Its
notable feature is that it returns `[value, derivative]`: the analytic derivative of the Lagrange
polynomial is computed alongside the value, which is how satellite velocity is obtained from
position-only SP3 products without finite differencing. `linear` has `double[]` and `long[]`
overloads because time is represented both ways in different parts of the codebase.

`Closest.findClosest_BS` uses `Collections.binarySearch` and correctly handles the negative
insertion-point return, including both boundary cases; `findClosest` is the O(n) equivalent for
unsorted `Double` lists. `Combination.getCombination(arr, n, r)` enumerates all size-*r* subsets,
used by the subset-testing and partial-ambiguity-resolution logic. `ArrayIndexComparator` sorts an
index array by the backing values — note it sorts **descending** (the comparator returns `+1` when
the first value is smaller), which is what you want for ranking satellites by elevation or C/N₀ but
is the opposite of the usual convention.

## String handling — `StringUtil` and `SymbolToken`

Both exist because GNSS text formats are fixed-column Fortran output, not delimited data.

`StringUtil.splitter(str, lens...)` splits a line into trimmed fields of the given widths. The
`flag` overload controls what happens when the widths do not sum to the line length: with `flag`
true it reports an error and returns `null`; with `flag` false (the default) it treats the remainder
as an extra trailing field. That lenient mode is what makes it usable on records whose last column
has variable width — serial numbers in SINEX, comments in CLK headers — and it is used throughout
`Orbit`, `Clock`, `SINEX`, `IONEX`, and `DCB_Bias`.

`SymbolToken.split(token)` addresses a narrower problem: in RINEX navigation records, adjacent
fixed-width scientific-notation numbers can run together with no separating space when a negative
sign consumes it, so whitespace splitting yields a single unparseable token. It splits on the
`digit` followed by `-` boundary to recover the individual numbers.

The remaining `StringUtil` methods (`str2arr_D`, `str2mArr_D`, `str2arr_L`) parse Java's own
`Arrays.toString`-style bracketed output back into arrays, with an optional scale factor — used for
reading configuration values and previously-logged vectors.

## Plotting — `GraphPlotter`

`GraphPlotter extends ApplicationFrame` and is the project's complete visualisation layer, built on
JFreeChart. **Constructing one opens a Swing window**; there is no headless mode, and the
`ChartUtils.saveChartAsJPEG` calls are present but commented out throughout, so plots are viewed
interactively rather than written to disk. Use `MakeCSV` if you want persisted output.

The class is large — roughly two dozen constructors and forty-odd static `graph*` methods — because
it uses constructor overloading as its dispatch mechanism: the shape of the data you pass determines
the plot you get. There is no common interface, so the practical way to find the right entry point
is to match your data structure against the parameter types. The plots fall into a handful of
families:

- **Residual and innovation diagnostics.** Per-satellite residual and innovation time series from
  `SatResidual` maps, optionally overlaid with outlier markers and normal-distribution envelopes
  (`graphSatRes`, and the `HashMap<String, ArrayList<SatResidual>>` constructor). Also the
  statistical quality series: posterior variance of unit weight (`graphPostUnitW`) and redundancy
  numbers (`graphRedundancy`, `graphRedundancyPPP`), both computed per measurement type.
- **Position and velocity error.** East/North/Up error against a reference (`graphENU`, with an
  overload adding an error-bound band), trajectory overlays across multiple estimators
  (`graphTrajectory`), GNSS/INS comparison (`graphGnssIns`), true-error and outlier scatter
  (`graphTrueError`, `graphOutlier`), and dilution of precision alongside satellite count
  (`graphDOP`, `graphSatCount`).
- **Ambiguity and cycle slip.** Detected-slip time series per satellite (`graphCycleSlip`,
  `graphCycleSlipAllEst`), fixed-ambiguity counts (`graphAmbiguityCount`), and success-rate versus
  fix-rate comparisons across LAMBDA estimator types (`graphSRFR`).
- **PPP state histories.** `createPPPplots(...)` has one overload per PPP filter variant
  (`Rinex.estimation.EKF_PPP`, `EKF_PPP_DF`, `EKF_PPP_LowCostRx`, and the Android `EKF_PPP3`) and
  produces the standard set in one call: receiver clock offsets and drifts per constellation, zenith
  wet delay with its σ envelope, per-satellite ambiguity and ionosphere states, and position σ. This
  is the entry point you want after a PPP run, not the individual constructors.
- **Raw measurement and sensor inspection.** IMU triads (`graphIMU`), Android receiver clock and
  timing fields (`graphAndroidRawGNSStimeParams`), and time-differenced range checks
  (`graphDeltaRange`).

Two export methods live here rather than in `MakeCSV` because they are tied to plot data:
`exportADRStateToCSV` dumps Android ADR state flags per satellite, and `writeEpochJSONL` writes one
JSON object per epoch containing the full filter snapshot — position, ENU error, per-measurement
residuals and innovations, posterior variance of unit weight, redundancy, clock states, zenith wet
delay, DOP, σ values, per-satellite ambiguity and ionosphere states, GIM residuals, and cycle-slip
detections. That JSONL file is the most complete machine-readable record of a run and the right
starting point for external analysis in Python or a notebook.

## Analysis and export — `Analyzer`, `MakeCSV`, `Trajectory`

`Analyzer` sits between the estimators and the plots. `processAndroid` and `processIGS` are the two
entry points — one per satellite type, mirroring the split seen in `Weight` and `SatUtil`. Each
takes the per-epoch satellite map, a ground-truth position and velocity series, the estimated
position and velocity from each estimator, and a residual map, then computes truth-referenced error
series: geometric range and range-rate from the true position, per-satellite pseudorange and Doppler
residuals against truth, and the resulting statistics. Its output populates the maps that
`GraphPlotter` then renders, and it drives the plotting calls directly.

`getVel` and `getOriginalVel` derive velocity by differencing successive positions (`getVel` also
produces acceleration). Note that they difference position without dividing by the time step, so the
result is a per-epoch displacement — correct as a velocity only at a 1 Hz sampling rate.

`MakeCSV` is the persisted-output counterpart to `GraphPlotter`, and unlike the plotting layer it
writes files you can keep:

- `exportPPPresultsToCSV` / `exportAndroidPPPToCSV` — the full PPP result set, one row per epoch.
- `exportSatResToCSV` — per-satellite residuals organised by measurement type and estimator.
- `exportQualityMetricsToCSV` / `exportAndroidQualityMetricsToCSV` — elevation, C/N₀, DOP,
  redundancy, and outlier flags per satellite per epoch.
- `exportENUToCSV` — local-frame error series for multiple estimators side by side.
- `exportDenseCSV` / `exportSparseCSV` / `exportScalarCSV` — three generic shapes for
  `TreeMap`-keyed time series: fixed-width `double[]` rows, sparse string-keyed maps (columns
  discovered from the data), and a single scalar per epoch. Reach for these when adding a new
  diagnostic rather than writing another bespoke exporter.

`Trajectory.createCSV` writes one wide CSV with a column triple per estimator for both position and
velocity, so multiple estimators' trajectories can be compared in a spreadsheet. It uses `-999` as
its missing-data sentinel and emits empty cells for those epochs.

## Extending this

- **Adding a new diagnostic output.** Prefer `MakeCSV.exportDenseCSV` / `exportSparseCSV` /
  `exportScalarCSV` over a new bespoke exporter — they cover most `TreeMap`-keyed time series
  already. If the diagnostic belongs to a filter run as a whole, add a field to
  `GraphPlotter.writeEpochJSONL` instead, so it travels with everything else from that epoch.
- **Adding a plot.** Follow the existing pattern: a static `graph*` method that takes the data
  structure and constructs the frame, rather than another constructor overload — the overload set is
  already hard to navigate. If you want the chart written to disk, uncomment (or add) the
  `ChartUtils.saveChartAsJPEG` call; it is present but disabled throughout.
- **Changing the stochastic model.** `Weight.computeCoVariance` is the single point of definition;
  everything else in the class composes it. Change it there rather than at call sites.
- **Adding a coordinate transform.** Put it in `LatLonUtil` alongside the existing frame
  conversions, and match the established signature conventions: reference position as `double[]`,
  an `isPos` flag for point-versus-vector semantics, and an explicit degrees/radians boolean where
  ambiguity is possible.
- **The duplicated-`Satellite` problem.** `Weight`, `SatUtil`, `Analyzer`, and `MakeCSV` all carry
  parallel methods for `Rinex.models.Satellite` and `Android.models.Satellite`. If you add a method
  that takes a satellite list, expect to add both variants — or, better, extract the accessors these
  methods actually need into a shared interface and take that instead.
- **Adding a numeric helper.** `MathUtil` and `Vector` are the right homes for scalar and 3-vector
  operations respectively; anything involving `SimpleMatrix` belongs in `Matrix`. Note that several
  existing helpers are hard-coded to three components (`MathUtil.getEuclidean`, `Vector.dotProd`)
  regardless of input length — do not assume they generalise.

## See also

- [architecture.md](architecture.md) — how `utility/` relates to the rest of the packages
- [corrections-and-models.md](corrections-and-models.md) — the physical models built on
  `LatLonUtil`, `Time`, and `Vector`
- [file-formats-and-parsers.md](file-formats-and-parsers.md) — every parser depends on
  `StringUtil`, `SymbolToken`, `Time`, and `Interpolator`
- [kalman-filters.md](kalman-filters.md) — how `Weight` and `Matrix` feed the filter's stochastic
  model and variance component estimation
- [lambda-ambiguity-resolution.md](lambda-ambiguity-resolution.md) — where
  `Matrix.computeDecorrelationNumber` and `Combination` are used
- [ppp-engine.md](ppp-engine.md) — the filters whose state histories `GraphPlotter.createPPPplots`
  and `MakeCSV` render and export
