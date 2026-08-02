# Corrections and Physical Models

The `com.gnssAug.helper` package holds the physical and geometric models that turn a raw
GNSS observable into something a positioning filter can use: where the satellite was when it
transmitted, where it sits in the sky relative to the receiver, and how much the ionosphere,
troposphere, and a deforming Earth have distorted the measurement in between. These classes
are deliberately stateless-or-nearly-so and know nothing about estimators — they are called by
the pipeline entry points (`IGS.IGS`, `Android.Android`, `Android.Android_Static_Rinex`) and by
the PPP filters, which decide whether a given correction is subtracted from the measurement or
carried as an estimated state. Several models delegate to [Orekit](https://www.orekit.org/) for
ephemerides, geoid undulation, and IERS conventions, so Orekit's `DataProvidersManager` must be
initialised (via `MainApp.buildGeoid()`) before most of them will work.

## Class index

| Class | Responsibility |
|---|---|
| `helper.ComputeSatPos` | Propagates broadcast ephemeris to satellite ECEF position, velocity, and clock offset |
| `helper.ComputeEleAzm` | Converts a receiver/satellite ECEF pair to elevation and azimuth in the local ENU frame |
| `helper.ComputeIonoCorr` | Klobuchar broadcast ionosphere model, driven by RINEX NAV header coefficients |
| `helper.ComputeIPP` | Thin-shell ionospheric pierce point geometry, shared by the GIM interpolator |
| `helper.ComputeTropoCorr` | Saastamoinen zenith delay from tabulated meteorology plus Niell elevation mapping |
| `helper.ComputeSolidEarthTide` | IERS 2010 solid Earth tide station displacement via Orekit |
| `helper.Rotation` | Direction-cosine-matrix / Euler-angle conversions and DCM re-orthonormalisation |

Two further corrections live outside this package because they are inseparable from the
product files they are parsed from — satellite antenna phase-center offsets and carrier-phase
wind-up (`Rinex.fileParser.Antenna`), and code/phase hardware biases
(`Rinex.fileParser.DCB_Bias`, `Rinex.fileParser.OSB_Bias`). They are summarised below and
documented in full in [file-formats-and-parsers.md](file-formats-and-parsers.md).

## Satellite position and clock — `ComputeSatPos`

`computeSatPos(NavigationMsg, tSV, tRX, ISC)` is the broadcast-ephemeris path: the standard
GPS interface-control-document Keplerian propagation, with a fixed 10-iteration Newton solve
for eccentric anomaly, WGS-84 values for the gravitational parameter and Earth rotation rate,
and the `F·e·√A·sin(E)` relativistic clock term. It returns an `Object[]` bundling four things
the callers need together: ECEF position with satellite clock offset, the satellite velocity
vector, a derived clock drift, the corrected transmission time, and the position rotated into
an inertial (ECI) frame at the reception epoch.

Two details are worth knowing before you touch this class:

- **Velocity is analytic, not differenced.** The velocity vector comes from differentiating the
  ICD orbital-frame position and the ICD→ECEF rotation matrix, not from finite differences of
  successive positions. This matters for Doppler-based estimation, which needs a velocity that
  is consistent with the position at the same instant.
- **The ECI output exists for Sagnac handling.** `xk_ECEF` is the position at transmission in
  the ECEF frame of that instant; `ECI` re-rotates it by `OMEGA_E_DOT · (tRX − t)` so that
  geometric range is computed in a single consistent frame. Downstream code (`SatUtil.getUnitLOS`,
  the EKF design matrices) reads `getSatEci()`, not the raw ECEF.

When precise MGEX products are enabled, this class is bypassed entirely: `Rinex.fileParser.Orbit`
Lagrange-interpolates SP3 positions and `Rinex.fileParser.Clock` linearly interpolates CLK clock
biases instead. `IGS.SingleFreq` is the switch point between the two paths.

## Elevation and azimuth — `ComputeEleAzm`

`computeEleAzm(userECEF, satECEF)` rotates the user→satellite range vector into the local ENU
frame using the standard geodetic rotation (Soler & Eisemann formulation, referenced in the
source) and returns `[elevation, azimuth]` in radians. It is the cheapest and most widely called
model here: elevation masking, C/N₀-and-elevation weighting (`utility.Weight`), the ionosphere
obliquity factor, the troposphere mapping function, and adaptive-filter binning all depend on it.

The azimuth is reconstructed from `atan(E/N)` with explicit quadrant fixes rather than `atan2`;
`Rotation.check` exists as a debug assertion that the two agree.

## Ionosphere

The codebase supports two mutually exclusive ionosphere paths, selected per run by whether an
IONEX product was loaded.

**Klobuchar broadcast model — `ComputeIonoCorr`.** Takes elevation, azimuth, user geodetic
position, receiver time, the eight `IonoCoeff` alpha/beta terms from the RINEX NAV header, the
signal carrier frequency, and a `Calendar` epoch. It computes the earth-central angle, the
pierce-point latitude/longitude in semicircles, geomagnetic latitude, local time, and obliquity
factor, then evaluates the cosine-series delay and scales it from L1 to the requested frequency
by `f_L1²/f²`. The same inputs are fed to Orekit's `KlobucharIonoModel` and the two results are
compared, with a warning printed to `stderr` on disagreement above 0.1 m — a cheap built-in
regression check rather than a fallback.

**GIM path — `ComputeIPP` + `Rinex.fileParser.IONEX`.** `ComputeIPP.computeIPP(elev, azm, lat,
lon, Re, h)` is the standalone thin-shell pierce-point geometry, parameterised by the shell
radius and height so it can be driven directly from the IONEX header rather than hard-coded
constants. It includes the polar special cases (|lat| > 70°) where the pierce point crosses the
pole and the longitude branch flips. `IONEX.computeIonoCorr` calls it, bilinearly interpolates
VTEC in latitude/longitude and linearly in time (with the standard co-rotating-longitude
correction between the two bracketing maps), applies the obliquity factor, and converts slant
TEC to metres of delay at the requested frequency.

**Where the choice is made.** `IGS.filterSat` prefers the GIM if an `IONEX` object was passed and
falls back to Klobuchar otherwise. In PPP mode the computed delay is *stored* on the `Satellite`
object (`setIonoErr`) but not subtracted — the filter estimates ionospheric delay as a state and
uses the stored value as an a priori pseudo-observation. In non-PPP modes it is subtracted from
the pseudorange and added to the carrier phase, reflecting the opposite sign of ionospheric group
delay versus phase advance.

## Troposphere — `ComputeTropoCorr`

Unlike the other models in this package, this one is instantiated per epoch:
`new ComputeTropoCorr(llh, time, geoid)` fixes the site latitude, day of year, and orthometric
height, then precomputes the zenith delays once; `getSlantDelay(elevation)` is then called per
satellite. Construction is where the cost is, so hoist it out of your satellite loop.

The zenith delay follows the UNB3m structure: average and seasonal-amplitude tables of barometric
pressure, temperature, relative humidity, temperature lapse rate, and water-vapour pressure height
factor, indexed by five latitude bands (15°, 30°, 45°, 60°, 75°) and interpolated linearly in
latitude, then modulated by a `cos(2π(D − Dmin)/365.25)` annual term with the phase reference
flipped between hemispheres. Those meteorological values feed a Saastamoinen-family formulation
that splits the result into hydrostatic (dry) and wet zenith delays.

Height enters twice and it is easy to conflate them. The constructor converts the ellipsoidal
height in `llh[2]` to orthometric height using Orekit's `Geoid.getUndulation` — that orthometric
height `H` is what the Saastamoinen and Niell height-correction formulas expect. Passing an
ellipsoidal height directly would introduce a geoid-undulation-sized error in the delay.

Mapping to the line of sight uses the Niell functions, implemented here as a Marini continued
fraction normalised to unity at zenith, with the standard height correction applied to the dry
mapping only. Orekit's `NiellMappingFunctionModel` is evaluated alongside as a cross-check (the
comparison assertions are currently commented out; a source note records that the dry `c`
coefficient differs from Orekit's while the rest agree).

`getSlantDelay` returns `[slantDelay, wetMappingFunction]`. The second element is the reason
this class returns an array rather than a scalar: PPP filters apply the full slant delay to the
measurement and simultaneously retain the wet mapping function as the partial derivative of the
estimated zenith wet delay state. `getM_wet(E)` exposes the same quantity on its own.

## Solid Earth tide — `ComputeSolidEarthTide`

Wraps Orekit's `TidalDisplacement` under IERS 2010 conventions, with Sun and Moon ephemerides
from `CelestialBodyFactory` and ITRF built with Earth-orientation parameters.
`calculateTimeVaryingTides(stationCoords, isLLA, dateTime)` returns the ECEF displacement vector
to **add** to the station position; the `isLLA` flag lets callers pass either geodetic or ECEF
coordinates.

The one configuration choice that materially changes the result is the final `false` argument to
the `TidalDisplacement` constructor (`removePermanentDeformation`). With it false, Orekit's
output already embeds the permanent tide, so no separate permanent-tide term may be added
elsewhere in the pipeline — doing so double-counts and produces a static vertical bias.
`getMeanTideCorrection` exists for the case where you deliberately want to move between the
tide-free and mean-tide conventions; it is a separate, explicit call and is not composed
automatically with `calculateTimeVaryingTides`.

Callers (`Rinex.estimation.EKF_PPP`, `Android.estimation.KalmanFilter.EKF_PPP`) project the
displacement onto each satellite's unit line-of-sight vector and add the resulting range
correction to both pseudorange and carrier phase, since a station displacement affects both
identically.

Ocean tide loading and pole tides are not modelled.

Note that `ComputeSolidEarthTide implements StationDisplacement` but its `displacement` override
returns `null` — the interface is vestigial and the class is used purely through its static
methods. Do not pass an instance to Orekit code expecting a working `StationDisplacement`.

## Rotation utilities — `Rotation`

Small, dependency-light helpers for attitude representation, used mainly by the INS/GNSS fusion
path rather than by the GNSS-only estimators:

- `dcm2euler` / `euler2dcm` — conversions in a yaw-pitch-roll ordering.
- `reorthonormDcm` — one iteration of the classic Gram-Schmidt-style drift correction that
  removes accumulated non-orthogonality from a DCM propagated by integrating gyro rates. If you
  are integrating attitude over long runs, this is what stops the matrix from slowly ceasing to
  be a rotation.
- `check` — a debug assertion comparing a manual quadrant-corrected `atan` against `atan2`.

## Antenna phase center and phase wind-up

`Rinex.fileParser.Antenna` sits with the parsers because it is fed by ANTEX data, but it
performs two corrections that belong conceptually here:

- **Phase center offset.** The satellite position from SP3 (or from broadcast ephemeris) refers
  to the antenna mechanical center. `getSatPC_windup` builds the nominal yaw-steering body frame
  (Z toward nadir, Y along the solar panel axis from `nadir × toSun`, X completing the right-hand
  set), rotates the ANTEX PCO vector into ECEF, and returns the corrected phase center position.
- **Carrier-phase wind-up.** The relative rotation between transmitter and receiver antennas
  changes the observed carrier phase; the correction (Wu et al. 1993) is unwrapped against the
  previous epoch's value to preserve continuity, which is why the caller must pass the stored
  previous wind-up.

**Eclipse exclusion.** During eclipse the satellite's yaw attitude departs from the nominal
model, so the body frame — and therefore both the PCO rotation and the wind-up — becomes
unreliable. `Antenna` detects umbra with a cylindrical shadow model and also catches the
noon/midnight yaw singularity, then returns `null` from `getSatPC_windup` for the eclipse
duration plus a 30-minute recovery window. **Callers must treat `null` as "drop this satellite
this epoch and clear its wind-up history"** — reusing a stale wind-up value across an eclipse gap
produces a phase discontinuity that cycle-slip detection will happily report as a real slip.

## Hardware biases

Broadcast and precise clock products are referenced to specific signal combinations, so raw
pseudoranges on other signals need a bias correction before they are consistent with the clock:

- **DCB (`DCB_Bias`)** — differential code biases, converted into per-signal inter-signal
  corrections and folded into the satellite clock inside `Clock.getBiasAndDrift(..., applyDCB)`.
  Per-constellation logic handles the GPS TGD/ISC relationship, Galileo BGD, and BeiDou's
  B1I/B3I-referenced convention.
- **OSB (`OSB_Bias`)** — observable-specific biases, covering code *and* phase, with validity
  windows and unit handling (`ns` versus `cyc`). Phase biases are what make integer ambiguity
  resolution possible at all in PPP-AR; without them the ambiguities are not integer-valued.

See [lambda-ambiguity-resolution.md](lambda-ambiguity-resolution.md) and
[ppp-engine.md](ppp-engine.md) for how this distinction propagates into the estimators.

## Where each correction is applied

| Correction | Applied in | Applied how |
|---|---|---|
| Satellite position/clock | `IGS.SingleFreq`, `Android.SingleFreq` | Broadcast (`ComputeSatPos`) or precise (`Orbit`/`Clock`) |
| Elevation/azimuth | `IGS.IGS`, `Android.Android` per epoch | Stored on each `Satellite`; drives masks and weights |
| Ionosphere | `IGS.filterSat`, `Android.Android` | Subtracted in non-PPP; stored as a priori state in PPP |
| Troposphere | `IGS.filterSat`, `Android.Android` | Full slant delay subtracted; wet mapping retained as a partial |
| Solid Earth tide | `Rinex.estimation.EKF_PPP`, `Android...EKF_PPP` | Projected on line of sight, added to code and phase |
| PCO + phase wind-up | `IGS.SingleFreq`, `Android.SingleFreq` | Corrected satellite position and phase term; `null` on eclipse |
| DCB / OSB | `Clock.getBiasAndDrift`, PPP filters | Folded into satellite clock / applied per observable |

## Extending this

- **Adding a new correction model.** Put it in `com.gnssAug.helper` as a class named
  `Compute<Effect>`. Follow the existing split: a static method if the model is a pure function of
  its arguments (`ComputeEleAzm`, `ComputeIPP`), an instantiable class if there is per-epoch or
  per-site setup worth amortising (`ComputeTropoCorr`). Return arrays rather than allocating a
  wrapper type — the rest of the package does, and the estimators unpack them positionally.
- **Applying it in the pipeline.** Non-PPP corrections belong in `IGS.filterSat` (and its Android
  counterpart) where elevation/azimuth are already available and the sign convention for code
  versus phase is applied in one place. Station-displacement corrections belong in the EKF loop
  where a reference position exists to project onto the line of sight.
- **Swapping the troposphere model.** Keep the `getSlantDelay` contract returning
  `[slantDelay, wetMappingFunction]`; the PPP filters depend on getting the wet mapping function
  back as an estimation partial, not just a delay.
- **Validating against Orekit.** `ComputeIonoCorr` and `ComputeTropoCorr` both evaluate an Orekit
  equivalent alongside the hand-rolled implementation and warn on divergence. That is a useful
  pattern to copy when adding a model that Orekit also implements — but keep it behind a flag or
  remove it once stable, since it doubles the per-satellite cost.
- **Adding a constellation.** Frequency tables live in `Android.constants.Constellation`, and
  bias handling in `Clock.getBiasAndDrift` and `DCB_Bias.getISC` is an explicit
  per-constellation `if` chain — both need a new branch.

## See also

- [ppp-engine.md](ppp-engine.md) — how the PPP filter consumes these corrections and which ones
  become estimated states instead of applied corrections
- [file-formats-and-parsers.md](file-formats-and-parsers.md) — the products these models read
  their inputs from (SP3, CLK, ANTEX, IONEX, DCB/OSB)
- [rinex-igs-pipeline.md](rinex-igs-pipeline.md) — the RINEX/IGS entry point and the order in
  which corrections are applied
- [android-pipeline.md](android-pipeline.md) — the smartphone entry point and its correction subset
- [utilities.md](utilities.md) — `LatLonUtil`, `Time`, `Weight`, and the other helpers these
  models build on
- [architecture.md](architecture.md) — how `helper/` relates to the rest of the codebase
- [thesis/08-corrections-and-models.md](thesis/08-corrections-and-models.md) — the research-track
  treatment, including the modelling rationale behind these choices
