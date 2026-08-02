# File Formats & Parsers

This codebase reads two broad categories of input: raw/consumer GNSS data (Android logs, Google Decimeter Challenge files) and geodetic/precise-product formats (RINEX and MGEX products). This document lists what's supported and where.

## Geodetic / precise-product formats (`com.gnssAug.Rinex.fileParser`, `com.gnssAug.IGS.models`)

| Format | Purpose | Parser | Model |
|---|---|---|---|
| **RINEX 3 Observation** | Code, phase, Doppler, C/N₀ per satellite/epoch | `ObservationRNX.rinex_obsv_process(path, useSNX, sinex_path, obsvCode, phaseReq, dopplerReq)` | `Rinex.models.ObservationMsg` |
| **RINEX 3 Navigation** | Broadcast ephemeris, Klobuchar iono coefficients, time corrections | `NavigationRNX.rinex_nav_process(path, getIonoOnly)` | `Rinex.models.NavigationMsg`, `IonoCoeff`, `TimeCorrection` |
| **SP3** | MGEX precise orbits (per-constellation/PRN ECEF position) | `Orbit.java` | `IGS.models.IGSOrbit` |
| **CLK** | MGEX precise clocks (per-constellation/PRN clock bias) | `Clock.java` | `IGS.models.IGSClock` |
| **ANTEX** | Satellite/receiver antenna PCO/PCV (currently `igs20.atx`) | `Antenna.java` (+ supplementary `antenna.csv`) | `IGS.models.IGSAntenna` |
| **IONEX** | Global Ionosphere Maps (GIM), 3D VTEC grid | `IONEX.java` | `vtecMap` (lat/lon/epoch-indexed) |
| **SINEX** | Station coordinates/ARP, phase-center offsets, monument eccentricity | `SINEX.sinex_process(path, siteCode, obsvCode)` | — |
| **DCB (SINEX BIAS `.BSX`/`.BIA`)** | Differential Code Biases — makes broadcast/legacy pseudoranges compatible with precise clocks | `DCB_Bias.java` | keyed by constellation/PRN/observable pair |
| **OSB (SINEX BIAS `.BIA`)** | Observable-Specific Biases — code *and* phase bias corrections for PPP-AR | `OSB_Bias.java` | `BiasEntry` (validity window, value, std, unit) |

See [docs/07](07-ppp-engine.md#standard-ppp-vs-ppp-ar) for why DCB vs. OSB correction determines whether PPP-AR is even possible.

### Product naming and retrieval

Precise products follow the IGS long-filename MGEX convention `AGENCY_TYPE_YYYYDDD0000_01D_RES_TYPE.EXT`, resolved per year/DOY by `MainApp.pDir()`. `tools/download_igs_products.py` automates fetching a full set (orbits, clocks, OSB/DCB, IONEX, SINEX, broadcast nav) for a given RINEX obs file from IGN/BKG/CDDIS mirrors, across a choice of analysis centres (WUM default, COD, GFZ, JAX).

## Consumer / Android formats (`com.gnssAug.Android.fileParser`)

| Format | Purpose | Parser |
|---|---|---|
| **Android `GnssLogger` raw log** | Raw pseudorange/carrier-phase/Doppler/ADR-state measurements + IMU samples, split on `Raw`/sensor markers | `GNSS_Log.process(...)` → per-epoch `GNSSLog` map + `IMUsensor` list |
| **Google Decimeter Challenge `derived.csv`** | Pre-corrected measurement-derived CSV; signal-name → observable-code mapping (`GPS_L1→G1C`, `GPS_L5→G5X`, `GAL_E1→E1C`, `GAL_E5A→E5X`, `GLO_G1→R1C`, `BDS_B1I→C2I`, `QZS_J1→J1C`, `QZS_J5→J5X`) plus a QZSS PRN remap table | `DerivedCSV.java` |
| **Decimeter Challenge ground truth** | Per-epoch `[GPStime, weekNo, lat, lon, alt, vel]` | `GroundTruth.java` |
| **Alternate ECEF ground truth** | `[GPStime, weekNo, x, y, z]`, whitespace-delimited (NMEA-GSA/reference-receiver style) | `GroundTruth_GSA.java` |

The **ADR (Accumulated Delta Range) state** field inside the raw Android log is the hardware-layer signal used by the [cycle-slip detection framework](05-cycle-slip-detection.md#hardware-layer-screening-android-adr-states) — it's parsed here and consumed there.

## Which entry point reads which format

| Entry point | Observation format | Precise products |
|---|---|---|
| `Android.Android` / `Android_Static` | Raw `GnssLogger` log (+ optional `derived.csv`) | Optional MGEX (SP3/CLK/IONEX/DCB/OSB) |
| `Android.Android_Static_Rinex` | RINEX 3 (from an Android RINEX export) | MGEX |
| `IGS.IGS` | RINEX 3 | MGEX + SINEX |

## Code index

| Class | Format |
|---|---|
| `Rinex.fileParser.ObservationRNX` | RINEX 3 OBS |
| `Rinex.fileParser.NavigationRNX` | RINEX 3 NAV |
| `Rinex.fileParser.Orbit` | SP3 |
| `Rinex.fileParser.Clock` | CLK |
| `Rinex.fileParser.Antenna` | ANTEX |
| `Rinex.fileParser.IONEX` | IONEX |
| `Rinex.fileParser.SINEX` | SINEX |
| `Rinex.fileParser.DCB_Bias` | SINEX BIAS (DCB) |
| `Rinex.fileParser.OSB_Bias` | SINEX BIAS (OSB) |
| `Android.fileParser.GNSS_Log` | Android raw `GnssLogger` |
| `Android.fileParser.DerivedCSV` | GDC `derived.csv` |
| `Android.fileParser.GroundTruth` / `GroundTruth_GSA` | GDC / alternate ground truth |

## Next

[docs/10 — Results Summary](10-results-and-recommendations.md) collects the headline numbers from every pillar in one place.
