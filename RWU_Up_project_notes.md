# Root Water Uptake — Up / Us / RSWF Project Notes

**Purpose of this document:** running reference for our ongoing analysis of plant-driven vs.
soil-driven root water uptake, built from a mechanistic Couvreur-model simulation treated as a
"virtual SWaP experiment." Re-upload this file at the start of future sessions to pick up where
we left off.

---

## 1. Project setup

- **Ground truth generator:** Vanderborght's mechanistic 1D Richards + Couvreur hydraulic
  architecture solver (SUF, Kr, kr, B, Krs computed from real root geometry/conductances).
  This produces `theta(z,t)` that contains *real* compensatory root water uptake, not just a
  passive proportional sink.
- **Virtual sensor experiment:** the simulated `theta(z,t)` is treated as if it were SWaP sensor
  output (van Dusschoten et al., 2020) — i.e. replicating their real experiment in silico.
- **Empirical extraction:** on top of that virtual data, `Up(z)` is estimated the way van
  Dusschoten's method does it — as the slope of a regression of `-dtheta/dt` vs. `Tact`
  (originally single-predictor; later extended with a linear-in-time drift term, and with
  selective point sampling instead of full-window continuous points).
- **Why this is a stronger test than the original paper's own validation:** van Dusschoten et
  al.'s in-silico check imposed a passive, non-compensating `Ûp(z)` shape as the sink. Our
  ground truth already contains real, mechanistically generated compensation — so this tests
  whether the empirical regression can recover the true `SUF(z)` / true instantaneous `RelSink`
  even when genuine compensation is active in the data-generating process.

---

## 2. Key definitions

| Symbol | Meaning | Source |
|---|---|---|
| `SUF(z)` / `SSF(z)` | Standard (Sink) Fraction — normalized uptake distribution under **uniform** soil water potential. Fixed property of root hydraulic architecture; independent of time, Tact, and soil water potential. In our script, approximated by normalized RLD (a valid simplification only when radial conductance is uniform per unit root length and axial resistance is negligible). | Couvreur, Vanderborght, Javaux (2012, HESS) |
| `Kcomp` | Compensatory RWU conductance — links local compensatory uptake to local deviation of soil water potential from the SUF-weighted mean (`Hs,eq`). | Couvreur et al. (2012) |
| `Up(z,t)` | **Plant-driven** RWU — the part of the sink that scales proportionally with total transpiration; by definition `Up(z,t) ≡ Tact(t)·SUF(z)`, i.e. compensation-free. | Couvreur et al. (2014a); van Dusschoten et al. (2020) |
| `Us(z,t)` | **Soil-driven** RWU redistribution — compensatory root flux; integrates to zero over depth by definition; independent of Tact instantaneously. | Couvreur et al. (2014a) |
| `RSWF(z,t)` | Pure soil capillary flux divergence (`∂/∂z[K(θ)(∂h/∂z − 1)]`) — **not** root-mediated. Kept separate from `Us` in our finer-grained decomposition (van Dusschoten's `Sr` lumps `Us` + `RSWF` together). | Extended Richards eq. (Couvreur 2014a); our own split |
| `Ûp(z)` | Normalized shape function: `Ûp(z) = Up(z,t)/Tact(t)`, integrates to 1 over depth. The regression target. | van Dusschoten et al. (2020) |
| `RelSink(z,t)` | `S(z,t)/Tact(t)` — **true, instantaneous** uptake fraction from the mechanistic sink (ground truth, no estimation). `RelSink(z,t) = SUF(z) + Us(z,t)/Tact(t)`. | Derived from Couvreur (2012) sink equation |
| `S(z,t)` | True mechanistic sink term = `Up(z,t) + Us(z,t)` (root-only, excludes RSWF). | Couvreur (2012) |

---

## 3. Governing equations

**Full Richards equation with mechanistic sink (Couvreur 2012):**
```
dtheta/dt = div[K * grad(Hs)] - S
S_k * V_k = Tact * SSF_k + Kcomp * (Hs_k - sum_j(Hs_j * SSF_j)) * SSF_k
```
Both capillary redistribution and compensatory sink solved exactly, simultaneously. **Nothing neglected.**

**Extended 1D form (Couvreur, Vanderborght, Beff, Javaux 2014a) — van Dusschoten's Eq. 1A:**
```
dtheta(z,t)/dt = Up(z,t) + Us(z,t) + d/dz[K(theta)*(dh/dz - 1)] + Es(t) at z=0
```
Still exact — `Up`/`Us` split made explicit, redistribution still fully included.

**Lumped, data-driven form (van Dusschoten et al. 2020) — Eq. 1B/2, actual fitting equation:**
```
dtheta(z,t)/dt = Up_hat(z) * Utot(t)/V + Sr(z,t) + Es(t) at z=0
Sr(z,t) ≈ p1(z) + p2(z)*t          <- linearized in time, valid only within short (<=6h) windows
```
`Sr` = `Us` + `RSWF` lumped together (can't be separated from theta(z,t) data alone).

**Our finer-grained board decomposition:**
```
dtheta/dt = Up + Us + RSWF
```

---

## 4. Redistribution — who neglects it

| Source | Redistribution treatment | Neglected? |
|---|---|---|
| Couvreur et al. (2012) | Full: capillary flux + compensatory sink both solved exactly | No — fully included |
| Couvreur et al. (2014a) | Same, `Up`/`Us` split made explicit inside Richards eq. | No — fully included |
| van Dusschoten et al. (2020) | `Us`+`RSWF` lumped into `Sr`, approximated linear-in-time within window | Partially — linearized, valid only for short windows (<=6h) |
| Our original single-predictor regression (`X=[1,Tact]`) | No drift term at all | Fully neglected (implicitly Sr=0) |

**Validity criterion (diffusion timescale argument):**
```
tau_soil ~ L^2 / D(theta),   D(theta) = K(theta)/C(theta)
Assumption valid when:  Delta_t << tau_soil
```
- Drier soil / lower K(theta) -> longer tau_soil -> assumption holds better
- Wetter soil / higher K(theta) -> shorter tau_soil -> assumption breaks down faster
- Short windows (4-6h) keep Delta_t << tau_soil even at moderate D(theta)
- Large established vertical gradients raise Sr magnitude/curvature -> more likely to violate linearity within a window

**Our advantage:** because theta(z,t) is simulated (not real sensor data), the TRUE Us(z,t) is
directly accessible — so the linear-in-time assumption can be tested directly per window/depth,
rather than only argued from the general diffusion-timescale reasoning above.

---

## 5. Three-way diagnostic (RelSink vs SUF vs Up_estimated)

```
RelSink(z,t) - SUF(z)         -> true root compensation (Us) only, pure physics
Up_estimated(z) - SUF(z)      -> total regression error (Us leakage + RSWF leakage)
Up_estimated(z) - RelSink(z,t)-> isolates RSWF leakage specifically (Us cancels out)
```
If (3) small even when (1),(2) large -> regression error is mostly real compensation leaking through (expected/fundamental limitation).
If (3) stays non-negligible even on low-compensation days -> soil capillary dynamics alone are already stressing the linear-in-t assumption (a distinct, non-root-related finding).

---

## 6. Whiteboard finding: dtheta/dt = Up + Us + RSWF, by depth

RSWF is not measured directly — it is inferred by subtraction:
```
RSWF(z,t) = dtheta(z,t)/dt - S(z,t)      where S = Up + Us (true mechanistic sink)
```

| Depth | dtheta/dt (observed) | S = Up+Us (root sink) | RSWF (inferred) |
|---|---|---|---|
| **TOP (5 cm)** | Strict decrease (abs. value) — mainly driven by Up itself declining | Decreasing | ~0 (nothing left nearby to redistribute from) |
| **MIDDLE** | ~0 — rises slightly then falls slightly (net ~unchanged); small "bump"/hook shape | Large, increasing | Large, decreasing — nearly exactly cancels S (near-total, near-instantaneous compensation) |
| **BOTTOM** | Slightly increasing | Slightly increasing | ~0 |

**Interpretation:**
- **TOP**: RSWF~0 -> observed signal directly, cleanly reflects true Up/S decline. No cancellation.
  This is where the observed dtheta/dt is the most honest reflection of true uptake.
- **MIDDLE**: S stays large and grows, but RSWF grows in the opposite direction, nearly perfectly
  cancelling it -> net signal looks deceptively calm while large two-way flux occurs underneath.
  The small "bump" (rise then fall) reflects imperfect/lagged synchronization between S and RSWF,
  not the plant's demand oscillating.
- **BOTTOM**: both drift up together gently, RSWF~0 — gradually picking up more uptake without
  much capillary support needed yet.

**Consequence:** the regression-based Up is most reliable at the TOP (where RSWF~0, signal
undisguised) and most compromised at the MIDDLE (where the true Up/S signal has been almost
entirely cancelled out of dtheta/dt before the regression ever sees it — a structurally
different failure mode than noise or slow drift).

**Open question (flagged, not yet resolved):** whether the bottom layer's RSWF~0 behavior
persists as the profile dries further on later days, or whether the bottom eventually develops
its own compensation signature once the middle layer depletes further.

---

## 7. Key literature (full citations)

1. Couvreur, V., Vanderborght, J., Javaux, M. (2012). A simple three-dimensional macroscopic
   root water uptake model based on the hydraulic architecture approach. *HESS* 16, 2957-2971.
2. Couvreur, V., Vanderborght, J., Draye, X., Javaux, M. (2014a — HESS). Horizontal soil water
   potential heterogeneity: simplifying approaches for crop water dynamics models. *HESS* 18,
   1723-1743.
3. Couvreur, V., Vanderborght, J., Draye, X., Javaux, M. (2014b — WRR). Dynamic aspects of soil
   water availability for isohydric plants: focus on root hydraulic resistances. *Water Resour.
   Res.* 50(11), 8891-8906.
4. van Dusschoten, D., Kochs, J., Kuppe, C.W., Sydoruk, V.A., Couvreur, V., Pflugfelder, D.,
   Postma, J.A. (2020). Spatially Resolved Root Water Uptake Determination Using a Precise Soil
   Water Sensor. *Plant Physiology* 184(3), 1221-1235.
5. Javaux, M., Couvreur, V., Vanderborght, J., Vereecken, H. (2013). Root Water Uptake: From
   Three-Dimensional Biophysical Processes to Macroscopic Modeling Approaches. *Vadose Zone
   Journal* 12(4).
6. Vanderborght, J., Couvreur, V., Javaux, M., et al. (2024). Mechanistically derived macroscopic
   root water uptake functions: The alpha and omega of root water uptake functions. *Vadose Zone
   Journal*.

---

## 8. Status / next steps (update this section each session)

- [x] Established Up definitions across Couvreur 2012 / 2014a / van Dusschoten 2020
- [x] Identified redistribution treatment differences (who neglects what)
- [x] Built RelSink vs SUF vs Up three-way diagnostic framework
- [x] Depth-wise Up/Us/RSWF sign analysis from simulation (top/middle/bottom)
- [ ] Resolve bottom-layer RSWF behavior on later days
- [ ] Quantify the "bump" phase lag between S and RSWF at the middle depth
- [ ] (Optional, parked) MATLAB implementation: selective-point regression + linear-time drift
      term, matching van Dusschoten's Eq. 2 exactly
