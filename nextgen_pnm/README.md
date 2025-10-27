# Next-Generation CaCO3→CaO Pore Network Model

This folder contains a self-contained rewrite of the spherical pore network
model tailored for CaCO₃ calcination with explicit handling of newly formed
micropores.  The original MATLAB scripts are left intact at the repository root;
all updated logic now lives under `nextgen_pnm/` so the new workflow can be
tested side-by-side with the legacy implementation.

## Highlights of the rewrite

- **Material-driven initialisation** – `buildSphericalPNM` now consumes a
  `materialSpec` structure (radius, BET areas, pore volume, etc.) and converts
  the experimental data into geometric targets via `materialToPNMParams`.
- **Multiscale pore classes** – micropores, mesopores, and macropores are
  tracked explicitly.  Class-specific coordination profiles and throat
  statistics are derived from SSA contributions, enabling denser micro-networks
  without losing the coarse backbone.
- **Mass-conserving calcination** – `simulateCalcination` follows the
  CaCO₃→CaO conversion, keeps a running solid mass budget, and spawns new
  micropores using literature-derived 10–100 nm radii when local CaCO₃ is
  depleted.
- **Network updates on-the-fly** – `newPoreInsertion` attaches newborn pores to
  their parent neighbourhood, re-computes throat geometry, and maintains
  connectivity while observing throat-radius bounds.
- **Extended QC** – `validatePNM` reports per-class statistics (coordination,
  porosity/SSA gaps) before and after reaction, making deviations from target
  measurements transparent.

Refer to the inline comments of each function for details on the new
assumptions and how they differ from the legacy scripts.

## Usage

```matlab
addpath('nextgen_pnm');
history = run_reaction_sim_V12();
```

The driver constructs the CaCO₃ network from the supplied material data,
simulates 300 s of calcination, and prints the final conversion, porosity,
and surface area while returning the full time history and the updated PNM.
