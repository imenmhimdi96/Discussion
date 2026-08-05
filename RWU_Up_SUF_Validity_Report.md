# Validity of the Plant-Driven Uptake (Ûp) Assumption under a Mechanistically Simulated Root–Soil System

### A Virtual SWaP Experiment — Depth- and Time-Resolved Analysis of Compensatory Root Water Uptake

**Prepared for:** Prof. Mathieu Javaux, Prof. Jan Vanderborght, Dr. Dagmar van Dusschoten
**Date:** update on each edit

---

## 1. Objective

The Soil Water Profiler (SWaP) methodology of van Dusschoten et al. (2020) estimates the plant-driven component of root water uptake, Ûp(z), from soil water depletion rate data alone, by regressing −dθ/dt against total transpiration (Tact) within a short (≤6 h) time window, under the assumption that soil water redistribution (Sr) can be approximated as linear in time within that window.

This report tests that assumption directly. Instead of applying the SWaP regression to real sensor data, we generate a fully mechanistic, ground-truth dataset using the Couvreur hydraulic-architecture root water uptake model (Couvreur et al., 2012; 2014a), coupled to a 1-D Richards equation solver. Because the true sink term — including its compensatory component — is known exactly in this simulated system, the SWaP extraction method can be validated against ground truth rather than assumed to work.

**Central question:**

> Under what depth and time conditions does the SWaP-estimated Ûp(z) equal the true Standard Uptake Fraction SUF(z), and when does it deviate — and why?

---

## 2. Governing Equations

Mechanistic sink term (Couvreur, Vanderborght, Javaux, 2012):

$$S(z,t) = Up(z,t) + Us(z,t) = SUF(z)\cdot Tact(t) \;+\; SUF(z)\cdot K_{comp}\cdot(H_{sr}(z,t) - H_{s,eq}(t))$$

Extended 1-D Richards equation (Couvreur et al., 2014a):

$$\frac{\partial \theta(z,t)}{\partial t} = Up(z,t) + Us(z,t) + \frac{\partial}{\partial z}\Big[K(\theta)\big(\tfrac{\partial h}{\partial z} - 1\big)\Big] \quad (+\ \text{evaporation at } z=0)$$

van Dusschoten et al. (2020) lump Us with pure soil capillary flux divergence into a single redistribution term Sr, linearized within each fitting window:

$$-\frac{d\theta(z,t)}{dt} = \hat{U}p(z)\cdot Tact(t) + Sr(z,t), \qquad Sr(z,t) \approx p_1(z) + p_2(z)\cdot t$$

The regression slope, Ûp(z), is theoretically identical to SUF(z) whenever Sr is separable from Tact within the window:

$$\hat{U}p(z) \equiv \frac{Up(z,t)}{Tact(t)} \equiv SUF(z)$$

Ground-truth diagnostic (needs no estimation — S(z,t) is known exactly from the mechanistic simulation):

$$RelSink(z,t) = \frac{S(z,t)}{Tact(t)} = SUF(z) + \frac{Us(z,t)}{Tact(t)}$$

---

## 3. Methodology — the Virtual SWaP Experiment

- **Ground truth:** 1-D Richards equation solved with the full Couvreur mechanistic sink (real SUF, Kr, kr, B, Krs from root geometry/conductances), driven by alternating high/low transpiration forcing over 7 days.
- **Virtual sensing:** the resulting θ(z,t) is treated exactly as SWaP sensor output — the real experiment of van Dusschoten et al. (2020) replicated in silico.
- **Extraction:** Ûp(z) estimated per depth, per 6-cycle window per day, by linear regression of −dθ/dt against Tact — both with the original single-predictor form and with selective time-point sampling (last-of-high-plateau, transition/ramp points, one low-plateau point).
- **Validation:** because Us(z,t) is known exactly in simulation, the SWaP linearity assumption can be tested directly, not only justified by the general diffusion-timescale argument.

This is a stronger test than van Dusschoten et al.'s own in-silico validation, which imposed a purely passive, non-compensating Ûp(z) shape as input. Here, real, mechanistically generated compensation is already active in the ground truth.

---

## 4. Depth-Resolved Analysis

RelSink(z,t) (colour-coded by minute), Ûp(z) (dashed), and SUF(z) (reference) across a well-watered day and a dried day:

![Day 1 (well-watered): RelSink / Ûp / SUF collapse onto a single curve across all six cycles](figures/day1_relsink_up_suf.png)
*Figure 1a. Day 1 — compensation negligible.*

![Day 6 (dried profile): RelSink fans out with time; Ûp develops a secondary deep peak diverging from SUF](figures/day6_relsink_up_suf.png)
*Figure 1b. Day 6 — compensation active, Ûp diverges from SUF at mid-depth.*

### 4.1 Why SUF is depth-dependent, and why this matters

SUF(z) is fixed by root architecture, confirmed by the near-identical green curves across all six cycles within each day. It is largest near the surface (≈2–5 cm) and falls to near-zero below ≈20–25 cm. Since Up(z,t) ≡ SUF(z)·Tact(t), this sets where compensation can dominate the sink: where SUF is large (top), the passive term starts large and compensation must work harder, proportionally, to matter; where SUF is small (bottom), even modest Us can dominate the local sink entirely.

### 4.2 Three-zone summary

| Zone | dθ/dt (observed) | S = Up+Us (true sink) | RSWF (soil-only, inferred) | Interpretation |
|---|---|---|---|---|
| **Top (≈5 cm)** | Strict decrease (\|·\|) | Decreasing | ≈ 0 | Nothing left nearby to redistribute from; observed signal directly tracks true Up decline. |
| **Middle (≈10–20 cm)** | ≈ unchanged (small rise-then-fall "bump") | Large, increasing | Large, decreasing — nearly cancels S | Near-total, near-instantaneous compensation; true uptake signal almost fully hidden from raw data. |
| **Bottom (>25 cm)** | Slightly increasing | Slightly increasing | ≈ 0 | Gradually absorbs uptake demand as upper layers deplete; least active zone until late drydown. |

RSWF(z,t) is not measured; it is obtained by subtraction: `RSWF = dθ/dt − S`, using the known true sink S from the mechanistic model.

### 4.3 Cycle-to-cycle (within-day) slope/intercept behaviour

| Depth | Slope Ûp (C1→C6) | Intercept Sr (C1→C6) |
|---|---|---|
| Top | Decreasing | Increasing |
| Middle | Increasing | Decreasing |
| Bottom | ≈ stable | ≈ stable |

The rising middle-depth Ûp is not a change in SUF (fixed) — it is compensation intensifying and migrating toward mid-depth as the top dries, captured by the regression because it is correlated with Tact in a consistent, cycle-over-cycle manner.

---

## 5. Time-Resolved Analysis (Within a Single Cycle)

![Day 4 regression at 10/20/30/40 cm — looped trajectory at high Tact, most visible at 20 cm](figures/day4_regression_10-40cm.png)
*Figure 2a. Day 4.*

![Day 6 regression — hook shape more pronounced, R² collapsing at deeper layers](figures/day6_regression_10-40cm.png)
*Figure 2b. Day 6.*

| Phase | Sink relation | vs. SUF | R² of fit | RSWF |
|---|---|---|---|---|
| **T-start** | S ≈ Up | S ≈ SUF | High | ≈ 0 |
| **T-middle** | S ≠ Up (transitional) | S deviates from SUF | Drops | Non-zero, uncertain sign |
| **T-end** | S ≈ Us (compensation-dominated) | S settles near Us, not SUF | Recovers | ≈ 0 again |

The fit is reliable at the two ends of a cycle. The middle of each cycle is where the lag between rising Up (demand) and rising Us (compensation catching up) produces a non-linear, non-single-valued (θ,Tact) trajectory.

$$(\hat{U}p - SUF) \approx [\text{component of } Us + RSWF \text{ that co-varies with } Tact(t) \text{ within the window}]$$

---

## 6. Cross-Cutting Synthesis

- **Depth axis:** reliable at the top (large SUF, RSWF≈0), worst at mid-depth (large, near-cancelling Us and RSWF).
- **Time axis:** reliable at cycle start/end (quasi-steady), worst mid-cycle (Up and Us dynamically chasing each other).

Worst case — mid-depth, mid-cycle, dried-down day — is precisely where the largest Ûp/SUF divergence and lowest R² occur.

| Comparison | Isolates |
|---|---|
| RelSink(z,t) − SUF(z) | True root compensation (Us) only — pure physics, no estimation |
| Ûp(z) − SUF(z) | Total regression error (Us leakage + soil-capillary/RSWF leakage combined) |
| Ûp(z) − RelSink(z,t) | RSWF leakage specifically (Us component cancels between the two terms) |

---

## 7. Final Answer — Validity of the SWaP Assumption

Controlled by soil hydraulic diffusivity `D(θ) = K(θ)/C(θ)` and characteristic redistribution timescale `τ_soil ≈ L²/D(θ)`; assumption holds when fitting window `Δt ≪ τ_soil`.

| Condition | Ûp(z) ≈ SUF(z)? | Physical reason |
|---|---|---|
| Well-watered profile (Day 1–2), any depth | Yes — near-exact | Us ≈ 0 everywhere; nothing to contaminate the fit. |
| Top layer, any day | Yes, throughout drydown | RSWF ≈ 0; observed signal directly tracks true Up. |
| Bottom layer, early–mid drydown | Yes | SUF(z) small, so absolute error stays small even under relative compensation. |
| Middle layer, once vertical Δh develops (Day 3+) | No — systematic deviation | Us and RSWF large, nearly cancel; regression signal mostly consumed by compensation. |
| Any depth, middle of a transpiration cycle | No — transient failure | Up/Us not yet quasi-steady; Sr not linear in t; loop/hook artefact. |
| Any depth, start/end of a transpiration cycle | Yes, even under drydown | Sink reaches quasi-steady regime; RSWF ≈ 0 at both bounds. |
| Very late, severe drydown (Day 7) | No, different reason | S dominated by Us, not Up — naive dθ/dt≈S≈Up fails structurally, independent of redistribution. |

**Summary:** the SWaP assumption is not simply "valid when wet, invalid when dry." It is valid wherever the true sink is either (a) genuinely dominated by the passive term (high SUF, low compensation — e.g. top layer throughout), or (b) has reached a quasi-steady state (cycle start/end). It fails wherever compensation is both large and actively evolving in step with Tact — intermediate depths, mid-cycle, intensifying with soil heterogeneity over time.

---

## 8. Conclusion and Recommendations

- Ûp(z) is a reliable proxy for SUF(z) only under a specific, identifiable combination of depth and cycle-phase conditions, not uniformly across the profile or day.
- Middle-depth, mid-cycle is the primary failure mode — near-total, dynamically-lagged compensation (Us), not noise or an arbitrary drift term.
- Restrict Ûp extraction to start/end phases of each transpiration cycle; interpret mid-cycle windows with caution or exclude them.
- Report the three-way diagnostic (RelSink−SUF, Ûp−SUF, Ûp−RelSink) alongside any future Ûp(z) estimate.

---

## References

1. Couvreur, V., Vanderborght, J., Javaux, M. (2012). A simple three-dimensional macroscopic root water uptake model based on the hydraulic architecture approach. *HESS*, 16, 2957–2971.
2. Couvreur, V., Vanderborght, J., Draye, X., Javaux, M. (2014a). Horizontal soil water potential heterogeneity: simplifying approaches for crop water dynamics models. *HESS*, 18, 1723–1743.
3. Couvreur, V., Vanderborght, J., Draye, X., Javaux, M. (2014b). Dynamic aspects of soil water availability for isohydric plants: focus on root hydraulic resistances. *Water Resources Research*, 50(11), 8891–8906.
4. van Dusschoten, D., Kochs, J., Kuppe, C.W., Sydoruk, V.A., Couvreur, V., Pflugfelder, D., Postma, J.A. (2020). Spatially Resolved Root Water Uptake Determination Using a Precise Soil Water Sensor. *Plant Physiology*, 184(3), 1221–1235.
5. Javaux, M., Couvreur, V., Vanderborght, J., Vereecken, H. (2013). Root Water Uptake: From Three-Dimensional Biophysical Processes to Macroscopic Modeling Approaches. *Vadose Zone Journal*, 12(4).
