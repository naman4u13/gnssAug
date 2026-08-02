# Results Summary & Recommendations

> Source: thesis Chapter 8, "Conclusion and Future Scope," plus the headline results from Chapters 3, 4, and 7.

## Headline results, one per pillar

| Pillar | Result |
|---|---|
| **1a — DBP** ([docs/03](03-pillar1-dbp.md)) | Under mismatched process noise, DBP reduces Horizontal RMSE by **57%** vs. standard Doppler-aided filters (17.89 m → 7.75 m) and eliminates **1,800+** false outlier detections — while matching VRWD's accuracy when `Q` happens to be well tuned. |
| **1b — AKF/VCE** ([docs/04](04-pillar1-akf-vce.md)) | Horizontal accuracy gains of **35.4% (static)** and **27.3% (vehicular)** vs. the best conventional baseline; only 67 of 66,185 measurements needed flagging as outliers, vs. 650+ for fixed-`R` filters — geometry is preserved instead of thrown away. |
| **2 — CSDR (detection + repair)** ([docs/05](05-cycle-slip-detection.md), [docs/06](06-lambda-estimators.md)) | Hybrid detection found **1,783** unique Galileo slips vs. 1,108 (chipset-only) or 1,243 (model-only) on real data. PAR-FFRT repair: **~79% repair rate** at **<1% failure rate** on L1+L5. |
| **3 — UU-PPP validation** ([docs/07](07-ppp-engine.md)) | On IGS data with 1,045 artificial slips: repaired accuracy statistically equals the clean baseline (Up RMSE 15.9 cm vs. 16.0 cm), vs. 52.9 cm for standard reset. On Pixel 4 real data: Galileo Vertical RMSE improved **45%** (150.1→82.2 cm); 2-minute cold-start convergence to sub-metre accuracy improved from **65% to 88%** probability. |

## The three pillars, in the thesis's own closing words

**Robust Dynamics (DBP):** *"Integrating Doppler measurements into the prediction step significantly enhances positioning robustness, particularly in high-dynamic scenarios. This approach decouples the dynamic model from erratic user behavior, ensuring state propagation is driven by the deterministic physics of the signal rather than in assumptions about motion."*

**Adaptive Stochastics (AKF/VCE):** *"The stochastic properties of smartphone GNSS data are volatile and environment-dependent. Static weighting models are obsolete in this context. Adaptive methods like VCE are essential to bridge the gap between the assumed error models of the Kalman Filter and the messy reality of urban signal propagation."*

**Hierarchical CSDR:** *"The ability to repair cycle slips transforms the smartphone positioning problem. It converts a fragmented series of short phase arcs into a continuous, connected observable, allowing the PPP filter to accumulate convergence over time despite frequent signal interruptions."*

## Implementation recommendations (Ch. 8.4, verbatim intent)

For anyone building on or adapting this codebase:

1. **Require dual-frequency (L1/L5 or E1/E5a).** Single-frequency cycle-slip repair is unreliable — the physics of ionospheric refraction mask the integer nature of the slip on one frequency alone.
2. **Decouple the code and phase clocks.** A single shared receiver clock state will diverge in the presence of the "dual-drift" behaviour documented in mobile chipsets; model `dt_r^Code` and `dt_r^Phase` (and their per-signal DSBs) as separate states.
3. **Use PAR-FFRT for ambiguity/cycle-slip resolution.** Full ILS is too risky at smartphone noise levels (~1-in-20 wrong fixes on L1); a fixed-ratio aperture test is too conservative (PAR alone repaired 0 of 153 L1 slips). PAR-FFRT is the balance point.
4. **Don't trust hardware ADR flags alone.** They're necessary (as a veto) but not sufficient (missed detections caused >28 m vertical error in one test) — always pair with a model-based (TDCP innovation) check.
5. **Propagate repair uncertainty, don't hard-fix.** Treat a repaired ambiguity as a weighted observation with its estimated variance, never as a deterministic constant (variance = 0).

## Device suitability checklist (from Ch. 7.5.5 / 8.3.3)

Before pointing this codebase's CSDR/PPP path at a new device, check for the three failure modes that ruled out most devices tested:

- **Cycle-slip density** — needs to be non-trivial for repair to matter statistically (Pixel 7 Pro: 0.07–0.49%, too clean to show a difference, though repair still worked correctly on the few slips present).
- **Clock steering behaviour** — some devices (OnePlus Nord 2) apply aggressive periodic code-clock steering that looks like a geometric error, not a cycle slip, and defeats integer repair entirely.
- **Dual-frequency availability** — required; see recommendation 1.

## Future scope (Ch. 8.5, as stated)

1. **CSDR's effect on full PPP-AR.** This thesis validates PPP-Float with repaired arcs; whether "stochastically repaired" arcs can support full Integer Ambiguity Resolution for instantaneous centimeter-level smartphone PPP is unexplored. (Architecturally supported already — see `EstimatorMode.IGS_PPP_AR` and the OSB bias path in [docs/07](07-ppp-engine.md#standard-ppp-vs-ppp-ar) — but not validated.)
2. **Kinematic validation.** All CSDR-PPP validation in Chapter 7 is static, to isolate clock/hardware artifacts from motion noise. Whether CSDR improves kinematic PPP convergence is future work.
3. **Tight IMU coupling.** This thesis is GNSS-only. Smartphone IMUs could plausibly bridge signal gaps and constrain the cycle-slip search space further, following the tightly-coupled GNSS/INS precedent from professional-grade receivers (see [docs/01's literature review notes](01-motivation-and-problem.md), Ch. 2.7).

## See also

- [docs/01](01-motivation-and-problem.md) — motivation, problem statement, publication list
- [docs/02](02-architecture.md) — how the pillars map to code
- [docs/03](03-pillar1-dbp.md), [docs/04](04-pillar1-akf-vce.md) — Pillar 1
- [docs/05](05-cycle-slip-detection.md), [docs/06](06-lambda-estimators.md) — Pillar 2
- [docs/07](07-ppp-engine.md) — Pillar 3
- [docs/08](08-corrections-and-models.md), [docs/09](09-file-formats-and-parsers.md) — supporting infrastructure
