# Session note: verifying the Us mechanism against couvreurnl.m

**Context.** From the supervisor meeting: Jan proposed that the deviation between
Ûp and relSink at cycle 3 is caused by a fast change in the soil water
gradient that drives Us, triggered by Tact. Task: verify this against the
simulation data rather than argue it from theory. This note is the working
log of that verification.

---

## 1. The exact formula behind Us, read from `couvreurnl.m`

The solver computes, per layer:

```
S(z) = SUF(z)·Krs·(ψeff − Hcollar) + Krs·SUF(z)·(ψint(z) − ψeff)
        \_____________________/     \_______________________/
              = SUF(z)·Tact                  = Us(z)
```

where `Hcollar = max(h_W, ψeff − Tp/Krs)`, so the first term is exactly
`SUF·Tact` (= `Up_true`), and the second term is `Us`.

**Confirmed numerically on the simulation output:**
```matlab
Us_true = S - SUF .* sum(S,1);
Us_form = Krs * SUF .* (HTINT - HTEFFINT(:)');
max(abs(Us_true(:)-Us_form(:))) / max(abs(Us_true(:)))   % = 1.68e-15
```
This is machine precision, not an approximation — `Us(z,t) =
Krs·SUF(z)·(HTINT(z,t) − HTEFFINT(t))` is exactly what the solver computes.

**Also confirmed:** `HTEFFINT = SUF' * HTINT` (SUF sums to 1, all four
arrays come from the same call) — mismatch ~7e-12 cm.

There is **no separate Kcomp parameter**. `Krs` (whole-root-system
conductance) plays that role directly, under the simplification the model
uses.

---

## 2. No rhizosphere storage - psi_int is solved algebraically, every minute

Read from `krsoilrootfunction2.m`: for each layer, `psi_int` is found by
solving

```
psi_int = (g_root*psi_x + g_soil*hT) / (g_root + g_soil)
```

with `lsqnonlin`, using only the **current** `Hsoil`, `phsoil`, `Hx` at
that same minute. Nothing from the previous minute enters this equation -
no volume, no water content, no time derivative for the rhizosphere
anywhere in this function.

**This ruled out an earlier, wrong guess** (that the fast response comes
from the rhizosphere being a small reservoir that fills/empties quickly).
That would require psi_int to depend on its own history. It doesn't. Us at
minute t is fully determined by the bulk soil state and the flux at
minute t, with no lag from storage.

---

## 3. The actual mechanism: g_soil collapses as the soil dries

Rearranging the equation above:

```
d = hT - psi_int = f * (hT - psi_x),      f = g_root / (g_root + g_soil)
```

- `g_root = rroot*kr` - fixed, a root property, does not change with soil
  wetness.
- `g_soil = B*K(h)` - depends on the soil's unsaturated hydraulic
  conductivity, which is highly nonlinear: near saturation it barely
  changes with h, but past a threshold it collapses by orders of magnitude
  for a small further drop in h.

**Wet soil:** `g_soil >> g_root` -> `f ~ 0` -> almost no drop across the
rhizosphere -> `psi_int ~ hT`.

**Dry soil:** `g_soil` collapses, becomes comparable to or smaller than
`g_root` -> `f -> 1` -> nearly the entire soil-to-xylem drop now happens
right at the root surface -> `psi_int` moves well below `hT`.

Because `f` is recomputed fresh every minute with no memory, once a layer's
local soil crosses into that steep, low-conductivity region, the same
Tact-driven flux produces a disproportionately large, **instantly
responsive** drop at the interface - this is the direct mechanism for why
the gradient that drives Us can start tracking Tact within a single
140-min cycle, rather than staying slow as Dagmar's assumption requires.

**Simple analogy:** same amount of water needs to move through a straw;
early on the straw is wide (g_soil large) and barely any suction is
needed; once the straw narrows (g_soil collapses), the same sip requires
much more suction (bigger potential drop) - not because more water is
being pulled, but because the same flow now needs much more "push"
through a worse conductor.

---

## 4. Three potentials - do not confuse these

| Symbol | Code variable | Meaning |
|---|---|---|
| `hT(z,t)` | `H(z,t) + Z(z)` | Bulk soil total head at depth z - the soil's own potential, untouched by the root |
| `psi_int(z,t)` | `HTINT(z,t)` | Potential right at the root surface at depth z, after the drop across the rhizosphere |
| `psi_eq(t)` (= psi_eff) | `HTEFFINT(t)` | SUF-weighted **average** of psi_int across all depths - one number, no depth - "what the plant as a whole feels" |

**Us compares psi_int to psi_eq, not to hT:**
```
Us(z) = Krs * SUF(z) * (psi_int(z) - psi_eq)
```
A layer can have Us > 0 (compensating) even if its own bulk soil isn't the
wettest in the profile, as long as its psi_int sits above the plant-wide
average psi_eq.

---

## 5. Important pipeline finding: H is stored one timestep ahead of S/HTINT/THETA

While checking `HTEFF` against `SUF'*(H+Z)` (off by up to 0.31 cm - not
negligible), a lag test found:

```matlab
% lag +1 gives max|diff| = 5.7e-13 (machine precision); all other lags don't
```

**Cause, confirmed from `richards1Dcouvreurnl.m`:** each timestep computes
`S`, `HTINT`, `HTEFF`, `HTEFFINT`, `HCOLLAR`, `THETA` from the head at the
**start** of the step, then solves for the new head and only afterward
sets `H = hp` (the new, end-of-step head). So at column index k, `H(:,k)`
is one minute **later** than the other six arrays at that same index.

**Not affected (same call, same index already):**
- `Us_true = S - SUF*Tact`
- The regression of `-dtheta/dt` on `Tact` (confirmed: `-dtheta_rate_change(:,k)`
  correctly pairs with `S(:,k)`, verified via the column mass-balance
  residual: shift 0 gives residual ~9e-16, shift 1 gives ~1e-3)
- Ûp, RelSink, RSWF, the whole existing report - all built on this correct
  pairing, nothing needs to be redone.

**Affected - needs a one-step correction before use:**
- Any calculation mixing `H` directly with `S`, `THETA`, `HTINT`, etc. at
  the "same" column index, e.g. the bulk-vs-rhizosphere gradient split
  (`hT - HTEFF`, `d - dbar`). Correct form: use `H(:,k-1)+Z` (or
  equivalently `hT_start = [H0+Z, H(:,1:end-1)+Z]`) as the head that was
  actually seen by step k.

**Still to verify:** whether `THETA` has the same alignment as the other
five (expected yes, since it comes from the same pre-step call as `S`);
and the value of the column-summed residual `sum(-dtheta/dt) - sum(S)`,
which should be constant in time if the domain conserves mass with no
top/bottom flux - if that constant is 0, it re-opens the earlier
RSWF-at-depth explanation (see point 6).

---

## 6. Open correction to an earlier RSWF explanation

The initial condition `H0` (linear from -200 to -155 cm) plus `Z` gives a
**uniform total head** (hydrostatic equilibrium) within about 1 cm across
the whole profile. This means the earlier explanation for the large,
day-1 RSWF signal at depth ("the initial condition is not at equilibrium,
so the profile is relaxing toward equilibrium from minute one") was
**wrong** - the initial condition already is in equilibrium.

The corrected, not-yet-confirmed hypothesis: the early RSWF at depth is
water moving **toward** the drying root zone once transpiration starts,
not an artifact of a non-equilibrium start. This fits the RSWF sign
pattern already documented (positive/large at depth from day 1, negative
higher up where roots are actively depleting). Needs the column-residual
check above to confirm the domain is closed (no boundary flux) before this
can be stated with confidence.

---

## 7. Status / next steps

- [x] Confirmed exact Us formula from the solver code (1e-15 match)
- [x] Ruled out rhizosphere storage as the fast-response mechanism
- [x] Identified g_soil collapse (via K(h) nonlinearity) as the actual
      mechanism, through `f = g_root/(g_root+g_soil)`
- [x] Found and diagnosed the H-vs-rest one-step time offset; confirmed it
      does **not** affect Ûp/RelSink/Us_true/the existing report
- [ ] Run the corrected bulk-vs-rhizosphere gradient split (`hT_start`
      aligned version) to determine whether the fast D3C3 signal comes
      from the rhizosphere drop specifically, or from bulk-soil
      heterogeneity too
- [ ] Confirm THETA alignment and the column-residual value (point 5/6)
- [ ] Once confirmed: re-run the mechanism plot (bulk part vs rhizosphere
      part vs Tact, days 1/3/5) with corrected alignment
- [ ] Supervisor's separate next step: simulate a different soil (steeper
      or flatter K(h) / retention curve) and check whether the D3C3 onset
      day shifts, as the g_soil-collapse mechanism predicts
