# Physical Correction Models

The estimators in Pillars 1–3 all consume a common set of geophysical/geometric corrections, implemented in `com.gnssAug.helper`. This document is the reference for what each one does and where it's used.

## Satellite position & clock — `ComputeSatPos.java`

`computeSatPos(NavigationMsg Sat, double tSV, double tRX, double ISC)` implements the standard GPS ICD broadcast-ephemeris-to-ECEF algorithm: Kepler orbit propagation with Newton iteration for eccentric anomaly, WGS-84 gravitational parameter `Mu` and Earth rotation rate `OMEGA_E_DOT`, plus the relativistic clock correction term. This is the fallback path when precise MGEX products aren't used (`useIGS = false`); when they are, `IGSOrbit`/`IGSClock` interpolation (see [docs/09](09-file-formats-and-parsers.md)) takes over.

## Elevation & azimuth — `ComputeEleAzm.java`

`computeEleAzm(double[] userECEF, double[] satECEF)` rotates the user-to-satellite line-of-sight vector into a local ENU frame (per NOAA's Soler & Eisemann formulation) and returns `[elevation, azimuth]` in radians. This feeds directly into elevation masking, the Klobuchar and GIM ionosphere models, tropospheric mapping functions, and the elevation/C-N₀ binning used by the [Adaptive Kalman Filter](04-pillar1-akf-vce.md).

## Ionosphere

Two independent paths, selected per run:

- **Klobuchar broadcast model** — `ComputeIonoCorr.java`: `computeIonoCorr(...)` computes the ionospheric pierce point (IPP) latitude/longitude, geomagnetic latitude, local time, and obliquity factor from elevation/azimuth/user position/time and the broadcast `IonoCoeff` coefficients (parsed from the RINEX nav header). Used when no GIM product is available.
- **IONEX Global Ionosphere Map (GIM)** — `Rinex.fileParser.IONEX.java` parses a 3D VTEC grid (lat/lon spacing, epoch interval, shell radius/height), and `ComputeIPP.java` provides the standalone thin-shell pierce-point geometry (`computeIPP(ElevAng_rad, AzmAng_rad, userLat_deg, userLong_deg, Re, h)`) used to interpolate VTEC at the IPP for a given satellite/epoch. This is the path used for the [PPP engine's](07-ppp-engine.md) GIM pseudo-observation constraint (`I_TECU^{s,T}` state, σ² = 9.0 TECU² a priori).

## Troposphere — `ComputeTropoCorr.java`

A Saastamoinen-family model using tabulated average and seasonal meteorological parameters (pressure, temperature, relative humidity, temperature lapse rate, water-vapour pressure height factor) indexed by latitude band — the same structure as the UNB3/UNB3m models referenced in the thesis (Ch. 7.2.1.1) — combined with Orekit's `NiellMappingFunctionModel` for elevation-dependent dry/wet mapping. In the PPP engine, the dry delay is applied directly and the wet delay's zenith component (`Z_w`) is carried as an estimated state (see [docs/07](07-ppp-engine.md#state-vector)).

## Solid Earth tide — `ComputeSolidEarthTide.java`

`implements StationDisplacement`; wraps Orekit's `TidalDisplacement` (IERS 2010 conventions, Sun/Moon ephemerides via `CelestialBodyFactory`) to compute the displacement vector to *add* to a station's ECEF position. `removePermanentDeformation = false` is used deliberately — Orekit's output in this mode already embeds the permanent tide, so no separate permanent-tide term should be added elsewhere in the filter (a ~2–5 cm static vertical bias from exactly this double-counting bug was found and fixed — see the "Fix permanent tide double-counting" commit). Ocean tide loading and pole tides are treated as negligible against the smartphone measurement noise floor and are not modelled.

## Phase wind-up & satellite eclipse exclusion — `Rinex.fileParser.Antenna.java`

Carrier-phase wind-up correction (`getSatPC_windup`) accounts for the relative rotation between transmitter and receiver antennas. During a satellite eclipse, the satellite's attitude control can behave unpredictably (yaw maneuvers), which corrupts wind-up computation — the antenna model implements **cylindrical umbra eclipse geometry with a 30-minute post-eclipse recovery window**, returning `null` during exclusion so callers (`IGS.SingleFreq`, `Android.SingleFreq`) can clear phase-continuity state for that satellite rather than feed it a corrupted correction.

## Antenna phase centers (PCO/PCV)

`Rinex.fileParser.Antenna.java` parses ANTEX-format calibration data (currently `igs20.atx`) plus a supplementary `antenna.csv` table into `IGS.models.IGSAntenna` objects — zenith/azimuth-dependent phase center offset and variation grids, with a receiver-PCO fallback for antenna models not present in the ANTEX file (see the "Switch to igs20.atx, add receiver PCO fallback" commit).

## Where these plug into the pillars

| Correction | Consumed by |
|---|---|
| Sat pos/clock, ele/azm | Every estimator ([docs/03](03-pillar1-dbp.md)–[07](07-ppp-engine.md)) — geometry design matrices, elevation masking |
| Klobuchar / GIM ionosphere | [PPP engine](07-ppp-engine.md) (`I_TECU` state + GIM pseudo-observation); single-frequency Android/IGS code positioning |
| Troposphere (Saastamoinen + Niell) | [PPP engine](07-ppp-engine.md) (`Z_w` state) |
| Solid Earth tide | [PPP engine](07-ppp-engine.md) station-coordinate correction, all EKF_PPP variants |
| Phase wind-up + eclipse exclusion | [PPP engine](07-ppp-engine.md) and [TDCP/cycle-slip filters](05-cycle-slip-detection.md) — corrupted wind-up during eclipse would otherwise masquerade as a cycle slip |
| Antenna PCO/PCV | Geometric range computation in `IGS.IGS` and `Rinex.estimation.*` |

## Code index

| Class | Role |
|---|---|
| `helper.ComputeSatPos` | Broadcast ephemeris → ECEF position/clock |
| `helper.ComputeEleAzm` | ECEF → elevation/azimuth |
| `helper.ComputeIonoCorr` | Klobuchar broadcast ionosphere model |
| `helper.ComputeIPP` | Thin-shell ionospheric pierce point geometry |
| `helper.ComputeTropoCorr` | Saastamoinen troposphere + Niell mapping |
| `helper.ComputeSolidEarthTide` | Orekit-based IERS 2010 solid Earth tide |
| `Rinex.fileParser.Antenna` | ANTEX PCO/PCV parsing, phase wind-up, eclipse exclusion |
| `Rinex.fileParser.IONEX` | IONEX GIM VTEC grid parsing |

## Next

[docs/09 — File Formats & Parsers](09-file-formats-and-parsers.md) covers the raw data formats these corrections are built from.
