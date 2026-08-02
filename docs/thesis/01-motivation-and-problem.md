# Motivation and Problem Statement

> Source: PhD thesis *"Improving Precise Smartphone GNSS with Robust Dynamics, Adaptive Stochastics, and Cycle Slip Repair"* (N. Agarwal, University of Calgary, Department of Geomatics Engineering, 2026), Chapter 1. DOI: [10.11575/PRISM/51108](https://doi.org/10.11575/PRISM/51108). See [CITATION.md](../../CITATION.md).

## Why this exists

GNSS is the primary source of global autonomous positioning, and the vast majority of GNSS receivers today are inside smartphones — **1.5 billion devices are produced every year**. Since 2016, Google has exposed raw GNSS measurements on most Android phones, opening the door to algorithms that go beyond the metre-level fixes smartphones traditionally produce. Recent literature suggests decimeter- or even centimeter-level positioning is possible on this hardware, which would unlock lane-level driving navigation, ubiquitous geo-surveying, precision agriculture, and precise augmented-reality applications — bringing survey-grade capability to commodity devices.

Reaching that precision means shifting reliance from noisy pseudorange (code) measurements to carrier-phase measurements, which are orders of magnitude more precise but ambiguous and fragile: when phase tracking is interrupted, a **cycle slip** occurs, introducing an integer discontinuity into the measurement.

> **For geodetic receivers, cycle slips are occasional anomalies. For smartphones, they are a chronic condition.**

Linearly polarized antennas, low-cost oscillators, and erratic user dynamics give smartphones frequent loss of lock and low carrier-to-noise density (C/N₀). The standard remedy — resetting the ambiguity term whenever a slip is detected — keeps the filter stable but is catastrophic for precision: every reset discards the accumulated geometric information and forces the filter back to square one.

**Central hypothesis of the thesis:** precise smartphone positioning cannot be achieved by simply discarding corrupted phase arcs. Instead, cycle slips must be *repaired* — the slip magnitude estimated as a stochastic value with an associated variance to bridge the discontinuity — while the filter's dynamics and stochastic models are made robust enough to keep that repair numerically stable in the first place.

## Three independent failure points

The thesis identifies three failure points in conventional smartphone positioning engines, and this codebase is organized around fixing each one:

### 1. Heuristic dynamics modelling

Conventional Kalman filters need *a priori* knowledge of the system dynamics (the process noise covariance, `Q`). Professional receivers are typically static or vehicle-mounted with predictable motion; smartphones move erratically — pocketed, handheld, running, driving — and standard heuristically-tuned models (e.g. constant velocity) fail whenever the user changes pace or direction. An incorrect dynamics model makes the filter diverge and cripples its ability to detect measurement outliers, because outlier tests depend on `Q` being right. There was no method that could adapt the filter's dynamic constraints to real smartphone motion instantaneously, without manual tuning.

**Codebase:** [Pillar 1 — Doppler-Based Prediction](03-pillar1-dbp.md), implemented in `Android.estimation.KalmanFilter.EKFDoppler` / `EKF`.

### 2. Static stochastic models

Precise positioning also depends on correctly weighting each measurement (the measurement noise covariance, `R`). Standard models assign weights from fixed elevation- or C/N₀-based look-up tables built for ideal open-sky conditions. A smartphone's environment changes rapidly — open sky to urban canyon — and a static weighting model either over-trusts noisy data or under-trusts good data. There was no way to estimate measurement variance in real time, independent of a pre-set model.

**Codebase:** [Pillar 1 — Adaptive Kalman Filter (VCE)](04-pillar1-akf-vce.md), implemented in `Android.estimation.KalmanFilter.AKFDoppler`.

### 3. The "reset" paradigm in carrier-phase positioning

To reach high precision, carrier phase has to be used, but low-quality smartphone antennas mean frequent cycle slips. Existing algorithms handle these by resetting the ambiguity term — stable, but it destroys precision by throwing away accumulated convergence. Detection methods existed, but there was a significant gap in *repair*: methods that estimate the float slip magnitude with its associated variance to maintain phase continuity.

**Codebase:** [Pillar 2 — Cycle Slip Detection](05-cycle-slip-detection.md) and [Stochastic Repair Engine / LAMBDA estimators](06-lambda-estimators.md), implemented in `Android.estimation.KalmanFilter.EKF_TDCP_ambFix2` and `helper.lambdaNew`.

## The three pillars

| # | Pillar | Contribution | Where |
|---|---|---|---|
| 1 | **Robust Dynamics Modelling** | Doppler-Based Prediction (DBP): use Doppler-derived velocity to drive state propagation, which automates `Q` estimation and improves both accuracy and Fault Detection & Exclusion (FDE) | [docs/03](03-pillar1-dbp.md) |
| 2 | **Adaptive Variance Estimation** | An Adaptive Kalman Filter (AKF) based on Variance Component Estimation (VCE) that learns the true measurement noise `R` in real time from the innovation sequence, binned by elevation/C/N₀ | [docs/04](04-pillar1-akf-vce.md) |
| 3 | **Cycle Slip Detection and Repair (CSDR)** | A hierarchical framework that treats cycle slips as an integer estimation problem rather than a reset trigger, comparing BIE, PAR, PAR-FFRT and other LAMBDA estimators, validated inside a custom Undifferenced Uncombined (UU) PPP engine | [docs/05](05-cycle-slip-detection.md), [docs/06](06-lambda-estimators.md), [docs/07](07-ppp-engine.md) |

A deliberate scope decision underlies pillar 3: this work targets **PPP-Float**, not PPP-RTK. Given how noisy and unpredictable smartphone measurements are, the primary challenge is simply maintaining a stable, continuous signal lock — attempting the rigorous integer resolution required for PPP-RTK is premature before that problem is solved. The thesis "deliberately focuses on establishing a reliable float ambiguity solution as a necessary prerequisite, laying the groundwork upon which future smartphone PPP-RTK implementations can be built." (The codebase's `EstimatorMode.IGS_PPP_AR` and the OSB-corrected code path exist and are architecturally supported — see [docs/07](07-ppp-engine.md#rank-deficiencies) — but full PPP-AR validation on smartphone data was out of scope.)

## Performance metric definitions (as used throughout)

- **Accuracy** — closeness of an estimated position to ground truth, quantified as RMSE against a post-processed high-precision reference trajectory.
- **Precision** — reproducibility/repeatability of the estimate, independent of closeness to truth; characterized by standard deviation.
- **Convergence** — the time for the PPP position solution to transition from metre-level (pseudorange-dominated) accuracy to decimeter-level (carrier-phase-dominated) accuracy, i.e. the time for RMSE to settle to a steady minimum.
- **Robustness** — the filter's resistance to divergence during abrupt maneuvers (sharp turns, rapid acceleration) where constant-velocity/acceleration models fail, combined with its ability to identify and de-weight erroneous data (e.g. severe multipath) without losing the true signal.

## Publications this codebase implements

1. N. Agarwal and K. O'Keefe, "Use of GNSS Doppler for Prediction in Kalman Filtering for Smartphone Positioning," *IEEE Journal of Indoor and Seamless Positioning and Navigation*, vol. 1, pp. 1–10, Nov. 2023. doi: [10.1109/JISPIN.2023.3337188](https://doi.org/10.1109/JISPIN.2023.3337188)
2. N. Agarwal, K. O'Keefe, and R. Klukas, "Alternative Approach to Integrate GNSS Doppler in Kalman Filter for Smartphone Positioning," *2023 13th IPIN*, Nuremberg, Sept. 2023, pp. 1–6. doi: [10.1109/IPIN57070.2023.10332485](https://doi.org/10.1109/IPIN57070.2023.10332485)
3. N. Agarwal and K. O'Keefe, "Application of Adaptive Kalman Filtering on Smartphone Positioning," *ION GNSS+ 2024*, Baltimore MD, Sept. 2024, pp. 2576–2588. doi: [10.33012/2024.19884](https://doi.org/10.33012/2024.19884)
4. N. Agarwal and K. O'Keefe, "Hybrid Cycle Slip Detection Method for Smartphone Global Navigation Satellite System," *Engineering Proceedings* (ENC 2024), vol. 88, no. 1, Art. 1, 2025. doi: [10.3390/engproc2025088010](https://doi.org/10.3390/engproc2025088010)
5. N. Agarwal and K. O'Keefe, "Cycle Slip repair for single-frequency smartphone GNSS using the Best Integer Equivariant estimator," *2025 IEEE/ION PLANS*, Monterey CA, Apr. 2025, pp. 1536–1543. doi: [10.1109/PLANS61210.2025.11028424](https://doi.org/10.1109/PLANS61210.2025.11028424)
6. N. Agarwal and K. O'Keefe, "Investigating Cycle Slip Repair for Single and Multi-Frequency Smartphone GNSS," *ION GNSS+ 2025*, Denver CO, Sept. 2025.

Full thesis: Agarwal, N. (2026). *Improving precise smartphone GNSS with robust dynamics, adaptive stochastics, and cycle slip repair* (Doctoral thesis, University of Calgary). DOI: [10.11575/PRISM/51108](https://doi.org/10.11575/PRISM/51108). https://hdl.handle.net/1880/124167

## See also

- [docs/02 — Architecture](02-architecture.md) for how these pillars map onto the codebase's packages and entry points.
- [docs/10 — Results Summary](10-results-and-recommendations.md) for every headline number in one place, plus the thesis's own implementation recommendations.
