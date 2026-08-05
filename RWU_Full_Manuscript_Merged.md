# Root Water Uptake Beyond Root Length Density: Insights from SWaP Experiments and Process-Based Modelling

> **Status note (remove before submission):** This draft merges the original project structure/literature review with the depth- and time-resolved Ûp/Us/RSWF analysis developed on the artificial (mechanistic) scenario. Sections describing the real *Vicia faba* SWaP/MRI experiment are carried over from the original draft largely as written (marked *[from original draft]*) and still need their own dedicated results; sections describing the artificial scenario are fully populated with the analysis completed so far. Passages marked **[TODO]** were left incomplete in the source draft and still need content.

---

## Abstract

A major challenge in plant and soil science is understanding how soil, root system, and plant hydraulics interact to control root water uptake (RWU). This study combines high-resolution Soil Water Profiler (SWaP) measurements, MRI (Magnetic Resonance Imaging) of root systems, and a one-dimensional (1D) soil–plant hydraulic model to quantify the spatial and temporal patterns of RWU in *Vicia faba* grown in two substrates: sandy soil and sandy loam. SWaP provided high-resolution vertical soil water content (θ) profiles under alternating light intensities, allowing derivation of water depletion rates and a data-driven, plant-driven uptake fraction (Ûp). MRI supplied root length density (RLD) profiles used to parameterize root system hydraulics. A macroscopic 1D Richards-based model incorporating soil, rhizosphere, and root hydraulic compartments (Couvreur et al., 2012) was used to simulate RWU under a controlled artificial drydown scenario, in which the true, mechanistic sink term is known exactly, and against which the SWaP-style Ûp regression could be validated directly rather than assumed to be accurate.

We show that the reliability of Ûp as a proxy for the true root-density-based uptake distribution (SUF) is neither uniform across depth nor constant over time. Under uniform soil water potential, Ûp coincides with SUF, indicating uptake is governed purely by root architecture. As the profile dries, Ûp progressively diverges from SUF and instead tracks the model's true relative sink (relSink), reflecting an increasing soil-driven (compensatory) contribution, Us. This divergence is depth-dependent — greatest at intermediate depths, where a substantial passive uptake baseline coincides with a substantial, actively evolving compensatory flux — and time-dependent within each transpiration cycle, being smallest at the start and end of a light period and largest during the transient mid-cycle phase. These findings indicate that the validity of the SWaP assumption is not simply a function of overall soil water status, but of a joint depth × cycle-phase condition, with direct implications for interpreting Ûp as a proxy for root distribution under field or growth-chamber conditions.

---

## 1. Introduction

Soil texture plays a key role in the dynamics of plant–soil interactions, shaping how water moves through the soil and its availability to plants. Different textures affect soil hydraulic properties, including hydraulic conductivity and soil-water diffusivity, which are essential for understanding water retention and movement in the soil profile (Gardner, 1960). Studies on kaolinitic clay soils have shown that even small differences within a particular soil type can have a significant impact on hydraulic properties; variations in microstructure, including organic matter and clay composition, critically influence hydraulic conductivity and water-holding capacity, introducing substantial uncertainty into soil hydraulic model parameterization and, consequently, water balance predictions (De Jong van Lier and Pinheiro, 2021). Cai et al. (2021) showed that soil texture significantly affects RWU by influencing the water potential at the soil–root interface: in sandy soils, hydraulic conductivity decreases more rapidly with decreasing water content than in loamy soils, producing a steeper decline in soil–root interface water potential, earlier stomatal closure, and a more rapid decrease in transpiration. Wankmüller et al. (2024) further showed that once soil water potential falls below a critical threshold (ψcrit), stomatal closure is triggered by declining leaf water potential and soil–plant hydraulic conductance, downregulating transpiration — a threshold reached at a rate governed by the interplay of atmospheric conditions, plant traits, and soil properties.

RWU represents the key process linking soil water availability to plant transpiration and therefore plays a central role in soil–plant–atmosphere interactions. Quantifying the spatial distribution of RWU within the soil profile is essential for understanding soil drying dynamics, plant water use, and drought responses. However, direct measurement of RWU remains challenging, particularly when resolving uptake at different soil depths; consequently, many studies estimate RWU indirectly from soil water depletion measurements, or approximate its spatial distribution from the vertical distribution of roots.

Methods to quantify and model RWU have evolved over the past six decades from simple empirical descriptions to more mechanistic representations. In macroscopic models, RWU is represented as a sink term within the Richards equation, with its distribution across soil layers based on factors such as root length density (RLD), without explicitly solving flow toward individual roots. Early formulations used empirical or semi-empirical functions such as the Feddes model, in which uptake at each depth is proportional to RLD and modulated by an α-stress reduction function (Feddes et al., 1978). Foundational work by Gardner (1960), Molz and Remson (1970), and Hillel et al. (1976) laid the conceptual groundwork later formalized by Feddes et al. (1978). Owing to its simplicity, this framework became a cornerstone of hydrological and crop models (e.g., SWAP, HYDRUS). Subsequent refinements addressed some limitations — notably the absence of compensatory uptake — including the Feddes–Jarvis model (Jarvis, 1989), a continuous S-shaped reduction function (van Genuchten, 1987), optimized root distribution functions (Vrugt et al., 2001), and formulations incorporating additional soil hydraulic properties (de Jong van Lier et al., 2008; Li et al., 2006). All such formulations, however, assume uptake in each layer is strictly local and cannot exceed the potential uptake allocated to that layer (Jarvis, 1989; Šimůnek and Hopmans, 2009). The simplicity and low computational cost of this class of models is particularly suitable for large-scale applications.

A second, functional-structural class of models instead defines an explicit root system architectural domain, incorporating root hydraulic features to simulate water flow toward individual roots (Doussan et al., 1998; Javaux et al., 2008). Their complexity makes them well suited to questions of root growth–soil interactions (Pagès et al., 2004; Somma et al., 1998), foraging for soil resources (Lobet et al., 2014; Lynch, 2013; Pagès, 2011), and plant responses in heterogeneous environments (Couvreur et al., 2014a; Huber et al., 2014). Couvreur et al. (2012) bridged these two classes with a macroscopic model that retains the low computational cost of the Feddes-type approach while incorporating compensatory uptake mechanistically, through an explicit hydraulic-architecture-derived sink term. This model allows water uptake to be partitioned into a passive, root-density-proportional component and a compensatory component driven by soil water potential heterogeneity — the theoretical basis adopted throughout the present study.

Building on this framework, van Dusschoten et al. (2020) introduced a purely data-driven method — the Soil Water Profiler (SWaP) — to estimate the plant-driven fraction of uptake, Ûp(z), directly from high-resolution soil water content time series, without requiring prior knowledge of root hydraulic architecture. This approach assumes that soil water redistribution varies approximately linearly with time within a short (≤6 h) fitting window, such that a regression of the local depletion rate against total transpiration recovers the true root-density-based uptake distribution. It remains an open question, however, whether this assumption — and the resulting Ûp estimate — continues to hold once genuine, root-architecture-mediated compensation (rather than only passive soil redistribution) is active in the system, and whether its validity depends on soil depth or on the phase of the diurnal transpiration cycle.

**Objective.** This study assesses the validity of the SWaP assumption directly, using an artificial drydown scenario simulated with the mechanistic Couvreur (2012) model, in which the true sink term — and hence the true compensatory contribution — is known exactly. We ask: (i) how does the data-driven Ûp(z), estimated by light-modulation regression, differ from the true root-architecture distribution (SUF) and from the true soil-driven compensation (Us) as the profile transitions from wet to dry; (ii) is this difference uniform across the profile, or does it depend on soil depth; and (iii) is the underlying SWaP assumption — that redistribution is negligible or linear within a short window — soil-condition dependent, or does it additionally depend on the phase of the transpiration cycle at which the estimate is made.

---

## 2. Methodology

### 2.1 Plant growth conditions *[from original draft]*

We grew faba bean plants (*Vicia faba*) in sandy and sandy loam soil. Seeds were kept in the dark and germinated for 2–3 days, then placed into PVC pipes 50 cm high with an inner diameter of 8.1 cm, filled with soil. **[TODO — growth conditions to be completed: growth chamber environment, watering regime prior to the experiment, number of replicates per substrate.]**

### 2.2 SWaP measurement principle and light modulation *[from original draft]*

SWaP measures Ûp using sensors composed of two opposing copper plates (7 × 5 cm²) connected to a coil, forming a resonance circuit. The sensors move along the z-axis and partially enclose the pots containing the soil columns. The resonance frequency of the circuit depends on the electrical permittivity of the material between the plates, primarily influenced by soil water content (θ); measuring the resonance frequency therefore enables determination of θ. The sensors move vertically in 1 cm increments, and a full scan of two pots with a 45 cm soil column was completed in approximately 10 minutes.

Before measurements began, the sensors were calibrated: the soil substrate was uniformly mixed with a known volume of water and placed in 12 cm high pots, which were then scanned to determine the resonance frequency. This was repeated across soil samples ranging from 2% to 41% water content, and the resulting data were fitted to a second-order polynomial curve.

SWaP measurements were coupled with controlled lighting using a water-cooled LED panel. Light intensity alternated between a high level (1000 μmol m⁻² s⁻¹) and a low level (500 μmol m⁻² s⁻¹) every two hours, over a total of 14 hours of illumination per day. Measurements were taken on 4-week-old plants; to simulate drought and monitor changes in plant hydraulic status, plants underwent a period of water limitation. The underlying hypothesis is that a sufficiently high-frequency alternation in light intensity generates a fast transpiration response, allowing short-term RWU response to be resolved: light modulation alters the xylem water potential (Ψx) but not the soil water potential (Ψsoil), and therefore does not itself affect water redistribution through root and soil. Raw SWaP resonance frequencies were deconvoluted and calibrated via quadratic polynomial fitting, separately for each soil texture, to obtain θ.

**Soil water depletion and total root water uptake.** The soil water depletion rate (dθ/dt) was calculated from the gradient of θ at each depth. Total root water uptake rate (Utot) was obtained by summing depletion rates across all depths. **[TODO — results of this real experiment (Tact/Utot response to light alternation, depth-resolved dθ/dt) to be added as a dedicated Results subsection; not yet available in this draft.]**

### 2.3 Soil hydraulic and root water uptake model

Water flow was simulated using the one-dimensional Richards equation, solved numerically with an implicit finite-difference scheme (van Dam and Feddes, 2000), coupled with the nonlinear root water uptake formulation of Couvreur et al. (2012), which partitions the sink term into a passive, root-architecture-driven component and a compensatory component arising from soil–root hydraulic resistances. Soil hydraulic properties (water retention and unsaturated hydraulic conductivity) followed the van Genuchten (1980) parameterization.

### 2.4 Artificial drydown scenario

To isolate the hydraulic response of a fixed root–soil system from any confounding effect of root growth, an artificial 7-day scenario was simulated representing a loam soil, parameterized as: residual water content θᵣ = 0.038 cm³ cm⁻³, saturated water content θₛ = 0.364 cm³ cm⁻³, shape parameters α = 0.0159 cm⁻¹ and n = 2.024, and saturated hydraulic conductivity Ksat = 67.4 cm d⁻¹. The soil profile was discretized into 45 layers of 1 cm thickness (nodes at 0.5–44.5 cm depth). The initial condition was a linear pressure head profile from h = −200 cm at the surface to h = −155 cm at the base, representing a moderately dry initial state with a slight wetting gradient with depth.

At the lower boundary, a seepage-face condition allowed the profile to switch between free drainage and a fixed-head condition depending on the pressure head at the bottom node relative to a threshold (hbot), representing a pot with limited drainage. At the upper boundary, no evaporation or precipitation was imposed, so that water could leave the profile only through root water uptake.

The root system was prescribed using a fixed, non-evolving RLD profile measured on a single reference faba bean plant ("plant 8") and held constant throughout the simulation — no root growth or change in the Standard Uptake Fraction (SUF) was permitted, isolating hydraulic response from architectural change. The RLD profile was top-heavy, with most root length concentrated in the upper 10 cm and a progressive decline with depth, becoming negligible below approximately 40 cm.

The potential transpiration rate (Tpot) was prescribed as an artificial, repeating diurnal forcing. During nighttime (20:00–06:00), Tpot = 0. During daytime (06:00–20:00), Tpot alternated hourly between a high-demand plateau (≈1.1 × 10⁻³ cm min⁻¹) and a low-demand plateau (≈3.2 × 10⁻⁴ cm min⁻¹), repeated seven times per day, producing a characteristic "comb" pattern; each transition passed through a brief ramp rather than an instantaneous step. This diurnal pattern was identical across all seven simulated days, so that any change in plant water status over the simulation arose solely from progressive soil drying, not from variation in atmospheric demand. Actual transpiration (Tact) was obtained by integrating the simulated sink term over the profile at each time step.

### 2.5 Root distribution and root/soil–root hydraulic conductances

For each layer, SUF was computed as the depth-weighted, normalized RLD profile:

$$SUF(z) = \frac{RLD(z)\cdot \Delta z}{\sum RLD(z)\cdot \Delta z}, \qquad \sum_z SUF(z) = 1$$

SUF was set to zero at the uppermost and lowermost nodes. Local root radius (rroot) was fixed at 0.05 cm; the mean inter-root distance within a layer was estimated from local root density under a Gardner-type upscaling, giving the dimensionless geometric parameter:

$$B(z) = \frac{2(\rho^2 - 1)}{2\rho^2(\ln\rho - \ln(1/0.53)) + 1 - 0.53\rho^2}$$

where ρ is the ratio of the soil-cylinder radius attributed to a single root (rbulk) to the root radius (rroot). Whole-root-system hydraulic conductance (Krs) was computed as:

$$Krs = \sum_z RLD(z)\cdot k_{rs}\cdot \Delta z$$

with root radial conductivity krs treated as a free calibration parameter. Local root-to-soil conductance components, Kr(z) and kr(z), were derived from Krs and SUF(z) following Couvreur et al. (2012).

### 2.6 Sink term formulation

At each time step, the local sink term S(z) was computed iteratively following the nonlinear stress function of Couvreur et al. (2012). The effective root collar pressure head was constrained between the prescribed Tpot divided by Krs and a fixed wilting-point pressure head (h_W = −15000 cm):

$$H_{collar} = \max\big(h_W,\; \psi_{eff} - Tpot/Krs\big)$$

where ψ_eff is the SUF-weighted average soil–root interface pressure head across all layers. The soil–root interface pressure head, and the corresponding nonlinear soil-to-root conductance, were obtained via a matric-flux-potential formulation (Schröder et al., 2008) at each layer and time step, solved iteratively (20 fixed-point iterations per step) for convergence of S(z), the interface pressure head, and Hcollar. Total actual transpiration was obtained by integrating the sink term:

$$Tact = \sum_z S(z)\cdot \Delta z$$

### 2.7 Estimation of the plant-driven uptake fraction (Ûp)

To assess whether the modelled uptake distribution was consistent with the prescribed root system and the model's own hydraulic dynamics, an independent, purely data-driven estimate of the relative contribution of each layer to total uptake was derived from the simulated θ(z,t) and Tact(t) time series alone, without reference to the prescribed RLD or the sink term computed by the model. For each layer and time step:

$$\frac{d\theta(z)}{dt} = \frac{\theta(z,t) - \theta(z,t-\Delta t)}{\Delta t}$$

Within each non-overlapping analysis window (2-hour windows, corresponding to one high- and one low-demand half-cycle of Tpot), a linear regression across all time steps in the window related the negative local depletion rate to Tact:

$$-\frac{d\theta(z)}{dt} = \beta_0(z) + \hat{U}p(z)\cdot Tact + \varepsilon$$

The slope, Ûp(z), was interpreted as the plant-driven uptake fraction attributable to layer z over that window. This regression was repeated independently for each of the 45 depth layers and each 2-hour window over the 7-day simulation, excluding nighttime periods (Tpot = 0), when no signal is available to constrain the fit.

### 2.8 Comparison variables and diagnostic framework

For each window, Ûp(z) was compared against two alternative descriptors: (i) the normalized RLD profile, nRLD(z) = SUF(z), the purely root-architecture-driven expectation under uniform soil–root conductance; and (ii) the model's own relative sink, relSink(z,t) = S(z,t)/Σ S(z,t), which already incorporates local soil water availability and hydraulic limitation. The coefficient of determination (R²) between Ûp(z) and each of nRLD(z) and relSink(z,t) was computed across all 45 layers per window, yielding two time series, R²(Ûp, nRLD) and R²(Ûp, relSink), tracking how well each descriptor predicted the true depth distribution of uptake as the profile progressively dried.

The soil-driven redistribution component was calculated as the residual between the model's total sink and the passive term:

$$Us(z,t) = S(z,t) - Up(z,t), \qquad Up(z,t) \equiv SUF(z)\cdot Tact(t)$$

and the pure soil capillary redistribution term was obtained by subtraction from the observed depletion rate:

$$RSWF(z,t) = \frac{d\theta(z,t)}{dt} - S(z,t)$$

Three diagnostic comparisons were used to separate genuine root compensation from soil-capillary estimation artefacts: RelSink(z,t) − SUF(z) (true compensation only); Ûp(z) − SUF(z) (total regression error); and Ûp(z) − RelSink(z,t) (soil-capillary leakage specifically, since the compensation component cancels between the two terms).

---

## 3. Results

*(Results below refer to the artificial drydown scenario, §2.4–2.8, for which the full analysis has been completed. Results from the real* Vicia faba *SWaP/MRI experiment (§2.1–2.2) are not yet included in this draft — see TODO markers above.)*

### 3.1 General behaviour: Tact, Tpot, and the well-watered–to–dry transition

Tact closely tracked the imposed Tpot forcing through the early part of the simulation, with visible deviation emerging as the profile dried in later days (quantified further in §4). This deviation coincides with the depth- and time-resolved changes in Ûp described below.

### 3.2 Depth-resolved behaviour of Ûp

Comparing RelSink(z,t), Ûp(z), and SUF(z) across the profile shows that under well-watered conditions (Day 1–2) the three quantities coincide across all six diurnal cycles, indicating negligible compensatory activity. As the profile dries (from approximately Day 3 onward), RelSink(z,t) fans out with time within each cycle, and Ûp(z) develops a secondary peak at intermediate depth absent from the fixed SUF(z) reference (Table 1).

**Table 1.** Depth-resolved decomposition of the observed depletion rate. RSWF is obtained by subtraction, using the known true sink from the mechanistic simulation.

| Zone | dθ/dt (observed) | S = Up+Us (true) | RSWF (residual) | Interpretation |
|---|---|---|---|---|
| Top (≈5 cm) | Monotonic decrease | Decreasing | ≈ 0 | Observed signal tracks true uptake directly. |
| Middle (≈10–20 cm) | Near-unchanged | Large, increasing | Large, decreasing — cancels S | Near-total compensation; true signal largely hidden. |
| Bottom (>25 cm) | Slight increase | Slight increase | ≈ 0 | Gradual, low-magnitude uptake transfer. |

### 3.3 Structural basis for the depth-dependence

Since Up(z,t) = SUF(z)·Tact(t), the absolute magnitude of the passive term at any depth is set entirely by local root density. Near the surface, where SUF(z) is largest, Up starts large, so early compensation — small in absolute terms while vertical potential gradients remain modest — stays negligible by comparison, and the near-surface sink remains Up-dominated. At depth, where SUF(z) approaches zero, the passive baseline is itself small, so even a modest Us could dominate; but because the vertical gradient has not yet propagated that deep during early-to-intermediate drying, Us also stays small there, and the sink remains quiescent. Intermediate depths occupy a structurally distinct regime, where SUF(z) is large enough to sustain a substantial passive baseline while simultaneously falling within reach of the developing gradient — producing the near-complete cancellation of Up and Us observed at these depths.

This also explains the cycle-to-cycle evolution of the fitted slope and intercept within a single day (Table 2): Ûp decreases at the top and increases at intermediate depth across successive cycles, while the intercept (Sr) shows the opposite trend — a genuine migration of compensatory activity from the surface toward mid-depth as upper layers progressively deplete.

**Table 2.** Within-day evolution of the fitted regression slope and intercept.

| Depth | Slope Ûp (Cycle 1→6) | Intercept Sr (Cycle 1→6) |
|---|---|---|
| Top | Decreasing | Increasing |
| Middle | Increasing | Decreasing |
| Bottom | Approximately stable | Approximately stable |

### 3.4 Linearity of Ûp against nRLD and against relSink

Early in the drydown, R²(Ûp, nRLD) was high across most depths, consistent with uptake being governed by root architecture alone. As drying progressed, R²(Ûp, nRLD) declined — most sharply at intermediate and deeper layers — while R²(Ûp, relSink) remained comparatively high, indicating that Ûp increasingly tracked the model's own compensation-inclusive relative sink rather than the fixed root-density profile.

### 3.5 Time-resolved behaviour within a transpiration cycle

The regression scatter of −dθ/dt against Tact developed a visible loop at later, drier days, rather than a straight line, most pronounced at intermediate depth and degrading R² sharply at 30–40 cm. This corresponds to a three-phase within-cycle evolution of the true sink (Table 3).

**Table 3.** Within-cycle phase decomposition of sink behaviour.

| Phase | Sink relation | vs. SUF | Fit quality (R²) | RSWF |
|---|---|---|---|---|
| T-start | S ≈ Up | S ≈ SUF | High | ≈ 0 |
| T-middle | S ≠ Up (transitional) | Deviates from SUF | Degraded | Non-zero |
| T-end | S ≈ Us | Settles near Us, not SUF | Recovers | ≈ 0 |

Only the fraction of Sr(z,t) that co-varies with Tact(t) within the window biases the fitted slope:

$$(\hat{U}p - SUF) \approx [\text{component of } Us + RSWF \text{ co-varying with } Tact(t) \text{ within the window}]$$

### 3.6 Us and RSWF

Consistent with §3.2–3.3, Us is negligible at the top and bottom of the profile in early drydown and becomes the dominant term at intermediate depth as drying progresses, while RSWF (pure soil capillary redistribution, isolated by the Ûp − relSink diagnostic) remains comparatively small throughout, indicating that the dominant source of Ûp/SUF divergence at intermediate depth is genuine root-mediated compensation rather than soil-physics estimation artefacts.

---

## 4. Discussion

When there is no stress (uniform soil water potential) and Tact/Tpot is close to 1, uptake is purely rooting-related: where roots are present, the plant simply takes up water, governed by soil water availability and root conductance; atmospheric demand is sustained by the plant with Hcollar decreasing, stomata open, and water pulled from the soil.

**When stress starts.** Tact follows Tpot during the first four days, indicating the plant can meet atmospheric demand without soil limitation. From Day 5 onward, as the soil dries, a deviation between Tact and Tpot indicates the onset of plant stress, becoming most visible on Days 6–7, when the average water content of the soil column decreases to approximately 0.09 cm³ cm⁻³. This decline in water content reduces soil hydraulic conductivity (Ks) and increases resistance to water flow through the soil–plant system, requiring a stronger gradient between plant and soil to extract water; Hcollar consequently becomes progressively more negative, increasing the water potential gradient between soil and plant.

Before stress, under wet soil conditions, Ûp follows SUF, indicating water uptake is mainly controlled by root architecture, with Ks sufficiently high to supply the roots without limitation. Following Vanderborght et al. (2023), the sink term can be expressed as:

$$S(z) = Up(z) + Us(z) \qquad \text{(Eq. 1)}$$

where Us is the soil-dependent component, negligible when soil water potential is uniform. As stress develops and differences in soil water potential emerge along the profile:

$$S(z) = Up(z) + K_{rs}\, SUF(z)\,\big(H_s(z) - \overline{H_s}\big) \qquad \text{(Eq. 2)}$$

the second term is no longer negligible. Consequently, the relative sink no longer follows SUF/RLD alone but reflects the combined influence of root distribution and soil hydraulic condition; as drying progresses, Ûp tracks relSink rather than RLD/SUF, indicating that root architecture is no longer the sole determinant of the observed uptake pattern.

**What is the real meaning of Ûp?** The present analysis extends this wet-to-dry framing along a second, previously unresolved axis: depth and cycle-phase. The findings above show that the answer to "is Ûp equal to SUF, or to relSink?" is not simply a function of *how dry* the profile is overall, but of *where in the profile* and *at what point in the transpiration cycle* the estimate is taken. Ûp remains close to SUF at the surface throughout the drydown (where the passive baseline dominates) and at the start/end of each transpiration cycle (where the sink has reached a quasi-steady state), regardless of the mean soil water status. It diverges toward relSink specifically at intermediate depths and during the transient mid-cycle phase — precisely where compensatory uptake is both large in magnitude and actively evolving on a timescale comparable to the imposed light modulation itself, violating the linear-in-time redistribution assumption underlying the SWaP regression.

This refines the original hypothesis (Ûp = RLD when wet, Ûp = relSink when dry) into a joint condition: **Ûp ≈ SUF wherever the passive uptake term dominates the local sink or the sink has reached a quasi-steady state; Ûp ≈ relSink wherever compensatory uptake is both substantial and actively evolving within the fitting window.** Soil dryness is the main driver of *when* this second condition is met, but depth and cycle-phase determine *where within the profile and diurnal cycle* it is met at any given level of overall dryness.

---

## 5. Conclusion

Using a mechanistic root water uptake simulation as ground truth, this study shows that the SWaP-estimated Ûp(z) recovers the true root-density distribution (SUF) reliably only under a specific, identifiable combination of depth and transpiration-cycle phase, rather than uniformly across the profile or the diurnal cycle, and only for as long as the soil water potential profile remains close to uniform. The intermediate-depth, mid-cycle regime constitutes the principal failure mode, driven by near-total, dynamically lagged compensatory uptake (Us) rather than by measurement noise or an arbitrary drift term in the regression. We recommend restricting Ûp(z) extraction to the start and end phases of each transpiration cycle, and reporting the three-term diagnostic (RelSink−SUF, Ûp−SUF, Ûp−RelSink) proposed here alongside any future SWaP-based estimate, to distinguish genuine root compensation from soil-capillary estimation artefacts. **[TODO: extend this conclusion once the real *Vicia faba* SWaP/MRI results (§2.1–2.2, sandy vs. sandy loam) are analyzed, to state whether the same depth × cycle-phase structure is observed experimentally, and whether it differs between the two substrates.]**

---

## References

Cai, G., et al. (2021). *(full citation to be completed)* — soil texture effects on root water uptake via soil-root interface water potential.

Couvreur, V., Vanderborght, J., Javaux, M. (2012). A simple three-dimensional macroscopic root water uptake model based on the hydraulic architecture approach. *Hydrology and Earth System Sciences*, 16, 2957–2971.

Couvreur, V., Vanderborght, J., Draye, X., Javaux, M. (2014a). Horizontal soil water potential heterogeneity: simplifying approaches for crop water dynamics models. *Hydrology and Earth System Sciences*, 18, 1723–1743.

De Jong van Lier, Q., Pinheiro, E.A.R. (2021). *(full citation to be completed)* — uncertainty in soil hydraulic parameterization and water balance prediction.

de Jong van Lier, Q., et al. (2008). *(full citation to be completed)* — root water uptake formulation incorporating soil hydraulic properties.

Doussan, C., Pagès, L., Vercambre, G. (1998). *(full citation to be completed)* — functional-structural root water uptake model.

Feddes, R.A., Kowalik, P.J., Zaradny, H. (1978). *Simulation of field water use and crop yield.* Wiley.

Gardner, W.R. (1960). Dynamic aspects of water availability to plants. *Soil Science*, 89(2), 63–73.

Hillel, D., Talpaz, H., van Keulen, H. (1976). *(full citation to be completed)* — root uptake and soil water dynamics.

Huber, K., et al. (2014). *(full citation to be completed)* — root water uptake in heterogeneous environments.

Jarvis, N.J. (1989). A simple empirical model of root water uptake. *Journal of Hydrology*, 107(1–4), 57–72.

Javaux, M., Couvreur, V., Vanderborght, J., Vereecken, H. (2013). Root Water Uptake: From Three-Dimensional Biophysical Processes to Macroscopic Modeling Approaches. *Vadose Zone Journal*, 12(4).

Javaux, M., Schröder, T., Vanderborght, J., Vereecken, H. (2008). *(full citation to be completed)* — three-dimensional root water uptake model (R-SWMS).

Li, K.Y., et al. (2006). *(full citation to be completed)* — root water uptake formulation.

Lobet, G., Pagès, L., Draye, X. (2014). *(full citation to be completed)* — root system foraging.

Lynch, J.P. (2013). *(full citation to be completed)* — root traits for resource acquisition.

Molz, F.J., Remson, I. (1970). *(full citation to be completed)* — root water uptake and soil moisture flow.

Pagès, L. (2011). *(full citation to be completed)* — root system architecture models.

Pagès, L., et al. (2004). *(full citation to be completed)* — root growth and soil interaction.

Schröder, T., et al. (2008). *(full citation to be completed)* — matric-flux-potential formulation for soil-root interface.

Šimůnek, J., Hopmans, J.W. (2009). *(full citation to be completed)* — root water uptake models in HYDRUS.

Somma, F., Hopmans, J.W., Clausnitzer, V. (1998). *(full citation to be completed)* — root growth and soil interactions.

van Dam, J.C., Feddes, R.A. (2000). *(full citation to be completed)* — implicit finite-difference Richards equation scheme.

van Dusschoten, D., Kochs, J., Kuppe, C.W., Sydoruk, V.A., Couvreur, V., Pflugfelder, D., Postma, J.A. (2020). Spatially Resolved Root Water Uptake Determination Using a Precise Soil Water Sensor. *Plant Physiology*, 184(3), 1221–1235.

van Genuchten, M.Th. (1980). A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. *Soil Science Society of America Journal*, 44(5), 892–898.

van Genuchten, M.Th. (1987). *(full citation to be completed)* — S-shaped stress reduction function.

Vanderborght, J., et al. (2023). *(full citation to be completed — cited in original draft for Eq. 1–2 sink partitioning)*.

Vrugt, J.A., Hopmans, J.W., Šimůnek, J. (2001). *(full citation to be completed)* — root water uptake distribution functions.

Wankmüller, F.J.P., et al. (2024). *(full citation to be completed)* — critical soil water potential threshold and stomatal closure.

---

*Note on references marked "full citation to be completed": these were cited by author/year only in the original draft; I have not independently verified their full bibliographic details and have not fabricated them. Please confirm or supply the complete citations before submission.*
