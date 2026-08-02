# gnssAug Documentation

Detailed technical documentation for the gnssAug GNSS positioning engine, tracing each module back to the PhD thesis it implements: *"Improving Precise Smartphone GNSS with Robust Dynamics, Adaptive Stochastics, and Cycle Slip Repair"* (N. Agarwal, University of Calgary, 2026). See [CITATION.md](../CITATION.md) and the top-level [README.md](../README.md) for the project overview.

## Reading order

| # | Document | Covers |
|---|---|---|
| 1 | [Motivation & Problem Statement](01-motivation-and-problem.md) | Why this exists, the three failure points in smartphone GNSS, the three-pillar solution, performance metric definitions, publications |
| 2 | [Architecture](02-architecture.md) | Package layout, `EstimatorMode` dispatch table, entry points, processing pipeline |
| 3 | [Pillar 1a — Doppler-Based Prediction](03-pillar1-dbp.md) | Complementary-filter theory, automated process-noise (`Q`) estimation, FDE, results |
| 4 | [Pillar 1b — Adaptive Kalman Filter (VCE)](04-pillar1-akf-vce.md) | Variance Component Estimation, elevation/C-N₀ binning, real-time `R` adaptation, results |
| 5 | [Pillar 2a — Cycle Slip Detection](05-cycle-slip-detection.md) | Android ADR hardware states, geometry-free/geometry-based detection, the 9-step hierarchical pipeline |
| 6 | [Pillar 2b — LAMBDA Estimators & Stochastic Repair](06-lambda-estimators.md) | ILS, IA-FFRT, BIE, PAR, PAR-FFRT — Monte Carlo variance estimation, comparative results, code map |
| 7 | [Pillar 3 — UU-PPP Engine](07-ppp-engine.md) | Observation model, state vector, rank deficiencies, CSDR integration, IGS + smartphone validation results |
| 8 | [Physical Correction Models](08-corrections-and-models.md) | Satellite position/clock, ionosphere, troposphere, solid Earth tide, phase wind-up, antenna PCO/PCV |
| 9 | [File Formats & Parsers](09-file-formats-and-parsers.md) | RINEX, SP3, CLK, ANTEX, IONEX, SINEX, DCB/OSB, Android raw logs, Decimeter Challenge CSVs |
| 10 | [Results Summary & Recommendations](10-results-and-recommendations.md) | Every headline number in one place, implementation recommendations, device suitability checklist, future work |

## Quick map: thesis chapter → doc → code

| Thesis chapter | Doc | Primary package |
|---|---|---|
| Ch. 1–2 (Intro, literature) | [01](01-motivation-and-problem.md) | — |
| Ch. 3 (Robust Dynamics / DBP) | [03](03-pillar1-dbp.md) | `Android.estimation.KalmanFilter.EKFDoppler` |
| Ch. 4 (Adaptive Variance / AKF) | [04](04-pillar1-akf-vce.md) | `Android.estimation.KalmanFilter.AKFDoppler` |
| Ch. 5 (Cycle Slip Detection) | [05](05-cycle-slip-detection.md) | `Android.estimation.KalmanFilter.EKF_TDCP_ambFix2` |
| Ch. 6 (Stochastic Repair Engine) | [06](06-lambda-estimators.md) | `helper.lambdaNew.*` |
| Ch. 7 (PPP Validation) | [07](07-ppp-engine.md) | `Rinex.estimation.EKF_PPP` |
| Ch. 8 (Conclusion) | [10](10-results-and-recommendations.md) | — |
