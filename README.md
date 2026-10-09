# Hybrid Cold Spray + Friction Stir Processing (CS + FSP) of an Al Coating on Cu: A FreeFEM++ 3D Eulerian (CEL-type) Thermo-Mechanical Study

<p align="center">
  <img src="https://img.shields.io/badge/FreeFEM++-Simulation-blue?style=for-the-badge&logo=gnu&logoColor=white"/>
  <img src="https://img.shields.io/badge/Al%20on%20Cu-Bimetallic-silver?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Cold%20Spray-Deposition-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/FSP-Stir%20Processing-red?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/3D-Eulerian%20%2F%20CEL-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Norton--Hoff-Viscoplastic-darkgreen?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/ParaView-VTK%20Export-green?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey?style=for-the-badge"/>
</p>

<p align="center">
  A <b>two-stage hybrid process</b> on a single fixed 3D mesh:
  a <b>cold spray nozzle</b> rasters over a <b>pure Cu substrate</b> and
  builds an <b>AA1100 Al coating</b> layer by layer, the part cools during
  a short <b>dwell</b>, and a <b>friction stir tool</b> then plunges through
  the as-sprayed coating into the Cu, traverses and lifts out. The script
  tracks coating growth, temperature, Al/Cu mixing, porosity closure and
  accumulated plastic strain, and exports everything to ParaView and CSV.
</p>

<img width="1008" height="772" alt="cs fsp" src="https://github.com/user-attachments/assets/9c635ef6-7b9c-4ec0-bfd9-9ea627d75d94" />


---

## Concept

Cold spray (CS) deposits metal powder in the solid state: particles
accelerated to several hundred m/s bond by plastic deformation on
impact, so the coating is built without melting. The as-sprayed deposit
is, however, **porous**, **work-hardened** and has **weak
inter-particle bonding**, and the coating/substrate interface is a sharp
mechanical bond.

Friction stir processing (FSP) applied **after** cold spray stirs the
deposit with a rotating tool: it closes pores, refines the grains,
breaks up particle boundaries and, if the pin reaches the substrate,
mixes coating and substrate across the interface. For Al on Cu this is
how Al–Cu composite or intermetallic-reinforced layers are produced.

This script is a numerical testbed for that sequence. Both stages run on
one Eulerian mesh that contains the substrate **plus an empty air layer
above it**; cold spray fills the air layer from below, and FSP then stirs
whatever material is there.

| Aspect | Description |
|---|---|
| Process | Cold spray raster deposition → dwell → FSP (plunge, traverse, lift) |
| Substrate | Pure Cu, 2 mm thick |
| Coating | AA1100 Al powder (or an Al + Cu blend via `fAlPow`) |
| Domain | Substrate + 3.5 mm void layer the coating grows into |
| Kinematics | 3D, fixed Eulerian mesh, tool and nozzle move through it |
| Flow model | Incompressible Stokes with Norton–Hoff viscosity (FSP stage only) |
| Tool model | Brinkman penalisation of a rotating, translating rigid pin + shoulder |
| Thermal model | θ-scheme advection–diffusion with mixed Al / Cu / void properties |
| Material tracking | Advected solid fraction, Al fraction, porosity, plastic strain |
| Elements | P1b–P1 (velocity–pressure), P1 for all scalar fields |
| Output | `.vtu` / `.pvd` fields for ParaView, CSV time history |

---

## Process Sequence

```
  Stage A  COLD SPRAY      nozzle rasters 4 layers x 11 passes (zig-zag),
                           coating grows from 0 to ~2.3 mm, no flow solved
  Stage B  DWELL           20 s cooling, thermal only
  Stage C  FSP  PLUNGE     pin tip from coating surface down to 0.5 mm into Cu
                TRAVERSE   tool moves in +x at 90 mm/min while rotating
                LIFT       tool retracts

  Pin length and plunge stroke are NOT fixed in advance: they are set
  after Stage A from the measured mean coating thickness on the FSP path.
```

---

## Geometry

```
      z
      |  ............................................  z = Hdom (top of void, air)
      |  .            void / air layer              .
      |  .     coating grows upward into this       .
      |  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  as-sprayed Al coating (~2.3 mm)
      |  ============================================  z = HzCu = 2 mm
      |  |              Cu substrate                |
      |  ============================================  z = 0  anvil (clamped, cooled)
    0 +-------------------------------------------------- x
        -Lx                                         +Lx

  Plan view (x-y):
        +Ly  +--------------------------------------------+
             |                                            |
             |   ==== coated band, width Wband = 30 mm ====|   raster passes along x,
       y=0   |   o--------------------------------------> |   hatch 3 mm, zig-zag
             |   ==========================================|
             |  xStart                              xEnd  |   FSP track along y = 0
        -Ly  +--------------------------------------------+

  Plate: 200 x 60 mm,  substrate 2 mm,  domain height 5.5 mm
```

The y-coordinate is mapped with a cubic stretch so elements are about
twice as fine at the plate centre (the FSP track) as at the edges.

---

## Physics

### Stage A: Cold spray deposition

Each sub-step the nozzle sweeps a segment `[xa, xb]` of length
`Lseg = vNoz * dtCS`. The Gaussian spray footprint is averaged over that
segment exactly with `erf`, so deposited mass is conserved for any time
step:

```
  gJet(x,y) = exp(-(y-yP)^2 / 2 s^2) * sqrt(2 pi) s / Lseg
            * 0.5 [ erf((x-xa)/(sqrt2 s)) - erf((x-xb)/(sqrt2 s)) ]

  dh/dt     = DE * mdot / (rho_coat (1 - poro0)) / (2 pi s^2) * gJet
  alpha     = clamp( (HzCu + h - z)/dz + 1 , 0, 1 )        solid fraction
```

Heat input during cold spray:

```
  Particle impact  q_KE  = etaKE * DE * mdot * vp^2 / 2 / (2 pi s^2) * gJet   [W/m^2]
  Gas jet (Robin)  q_gas = hJet0 * gJet * (Tgas - T)                          [W/m^2]
  New material     enters at Tpart (enthalpy mixed into T on deposition)
```

Both surface fluxes are applied to the top solid node layer of the
deposit (band of thickness one element).

### Stage C: FSP flow

```
  Stokes + Brinkman:
     -div(2 eta eps(u)) + grad p + rho/dt u + lam chiT (u - u_tool) = rho/dt u0 - rho (u0.grad)u0
      div u = 0

  Tool velocity:  u_tool = ( -omega y + vx_tool ,  omega (x - x_tool) ,  vz_tool )

  Viscosity (metal):  eta = K(T,phi) / ( 2 edot^(1-m) + eps )   clipped to [etaMin, etaMax]
  Viscosity (void):   eta = etaVoid
  Mixed:              eta_eff = alpha eta_metal + (1 - alpha) etaVoid

  Flow stress:        K(T,phi) = (phi K0Al + (1-phi) K0Cu) * exp( Q/R (1/T - 1/Tref) )
```

### FSP heat generation

```
  Shoulder (total)  P_sh  = delta mu sigma omega (2 pi / 3)(Rsh^3 - Rpin^3)
  Pin side          P_pin = delta mu sigma omega Rpin * 2 pi Rpin * L_immersed
  Viscous           q_v   = alpha (1 - chiT) 2 eta edot^2

  Both frictional terms are distributed on the mesh and re-normalised so the
  integrated power equals P_sh and P_pin; heat is applied only where alpha > 0,
  so the shoulder produces no heat until it touches the coating.
```

### Material state transport (FSP stage)

```
  dc/dt + u.grad(c) = Dart lap(c)       for c = alpha, phi, poro, epsAcc

  Porosity closure:   poro  <- poro * exp( -Cden * edot * dt )
  Plastic strain:     epsAcc <- epsAcc + edot * dt
```

### Thermal problem (both stages)

```
  rhoCp (dT/dt + u.grad T) = div(k grad T) + Q

  rhoCp = alpha (1-poro)[phi rhoCp_Al + (1-phi) rhoCp_Cu] + (1-alpha) rhoCp_void
  k     = alpha (1-1.5 poro)[phi k_Al + (1-phi) k_Cu]     + (1-alpha) k_void
```

### Boundary conditions

| Boundary | Label | Flow (FSP) | Thermal |
|---|---|---|---|
| Bottom, z = 0 | 1 | Clamped, u = 0 | Anvil, h = 300 W/m²K |
| Top of void, z = Hdom | 2 | Free | Air, h = 10 W/m²K |
| Front / back, y = ∓Ly | 3 / 4 | Clamped, u = 0 | Air, h = 10 W/m²K |
| Left / right, x = ∓Lx | 5 / 6 | Clamped, u = 0 | Air, h = 10 W/m²K |
| Tool volume | (penalised) | u = u_tool | Friction heat at contact |

Face labels are assigned explicitly with `change(Th, flabel=...)`, so
they do not depend on the default labelling of `cube()`.

---

## Material Properties

| Property | AA1100 Al | Pure Cu | Void (numerical air) | Unit |
|---|---|---|---|---|
| Density ρ | 2710 | 8960 | 10 | kg/m³ |
| Specific heat Cp | 900 | 385 | (ρCp = 1200) | J/kg·K |
| Conductivity k | 222 | 401 | 0.026 | W/m·K |
| Solidus / melting | 916 | 1358 | n/a | K |
| Norton–Hoff K0 | 1.5 × 10⁸ | 3.5 × 10⁸ | n/a | Pa·s^m |
| Strain-rate sensitivity m | 0.18 | 0.18 | n/a | n/a |
| Activation term Q | 1.2 × 10⁴ | 1.2 × 10⁴ | n/a | J/mol |

The temperature is capped at the Al solidus (916 K), since Al melting is
the first limit reached in this pair. The void density and heat capacity
are numerical values chosen for conditioning, not physical air data.

---

## Process Parameters

### Cold spray

| Parameter | Symbol | Default | Unit |
|---|---|---|---|
| Al fraction in powder | `fAlPow` | 1.0 | n/a |
| Powder feed rate | `mdotPow` | 3.3 × 10⁻⁴ (≈ 20 g/min) | kg/s |
| Deposition efficiency | `DE` | 0.70 | n/a |
| Particle impact velocity | `vPart` | 650 | m/s |
| Particle temperature at impact | `Tpart` | 400 | K |
| KE-to-heat fraction | `etaKE` | 0.9 | n/a |
| Gas temperature | `TgasJet` | 623 | K |
| Peak jet HTC | `hJet0` | 1500 | W/m²K |
| Spray spot Gaussian radius | `sigSpot` | 3 | mm |
| Nozzle traverse speed | `vNoz` | 50 | mm/s |
| Hatch spacing | `hatch` | 3 | mm |
| Coated band width | `Wband` | 30 | mm |
| Raster layers | `Nlayer` | 4 | n/a |
| Sub-steps per pass | `NsubCS` | 6 | n/a |
| As-sprayed porosity | `poro0` | 0.03 | n/a |
| Dwell time | `tDwell` | 20 | s |

### FSP

| Parameter | Symbol | Default | Unit |
|---|---|---|---|
| Shoulder radius | `Rsh` | 10 | mm |
| Pin radius | `Rpin` | 3 | mm |
| Pin penetration into Cu | `dPinSub` | 0.5 | mm |
| Shoulder plunge below coating | `dShPlunge` | 0.1 | mm |
| Rotation speed | `omegaTool` | 1200 | rpm |
| Plunge / traverse / lift speed | `vPlunge`, `vFeed`, `vLift` | 1.0 / 1.5 / 2.0 | mm/s |
| Friction coefficient | `muFric` | 0.35 | n/a |
| Slip fraction | `deltaFric` | 0.4 | n/a |
| Nominal contact pressure | `sigmaNnom` | 40 | MPa |
| FSP time steps | `NstepFSP` | 200 | n/a |
| Output interval | `Nsave` | 10 | steps |

### Numerical

| Parameter | Symbol | Default | Unit |
|---|---|---|---|
| Mesh divisions | `Nx`, `Ny`, `Nz` | 24, 12, 11 | n/a |
| Viscosity bounds | `etaMin`, `etaMax` | 10³, 10⁹ | Pa·s |
| Porosity closure rate | `Cden` | 1.0 | per unit strain |
| Artificial diffusion (scalars) | `Dart` | 10⁻⁷ | m²/s |
| Thermal θ-scheme | `thetaT` | 0.5 | n/a |
| Penalty coefficient | `lamPen` | 100 · etaMax / hmin² | Pa·s/m² |

---

## Derived Design Values (computed from the default inputs)

These follow directly from the parameters above and are printed to the
console at run time. They are **not** simulation results.

| Quantity | Value |
|---|---|
| Passes per layer | 11 |
| Time per pass | 3.8 s |
| CS sub-step / segment length | 0.633 s / 31.7 mm |
| Total cold spray time | 167.2 s |
| Expected mean coating thickness | ≈ 2.34 mm |
| Particle impact heat input | ≈ 44 W |
| Shoulder frictional power | ≈ 1.43 kW |
| Pin length used (if coating ≈ 2.34 mm) | ≈ 2.84 mm |
| FSP plunge / traverse / lift | ≈ 2.9 / 113.3 / 1.5 s |
| FSP time step | ≈ 0.59 s |
| Saved frames | 69 (1 initial + 44 CS + 4 dwell + 20 FSP) |

The pin length and FSP timings depend on the actual as-sprayed thickness
along y = 0, so the console values after Stage A are authoritative.

---

## Mesh

| Quantity | Value |
|---|---|
| Vertices | 3,900 |
| Tetrahedra | 19,008 |
| Stokes unknowns (P1b³ × P1) | ≈ 72,600 |
| Thermal / scalar unknowns (P1) | 3,900 |
| Element size, x | 8.3 mm |
| Element size, y (centre / edge) | 2.5 / 10 mm |
| Element size, z | 0.5 mm |

Treat the console line `Mesh: ... nodes ... tets` as authoritative; it
changes with any mesh parameter.

**Resolution caveat:** an 8.3 mm element in x barely resolves a 3 mm pin.
The penalised pin radius is inflated to about half an element so the
pin is always captured by at least one node column. The default mesh is
sized to run on a workstation; for quantitative mixing results use
`Nx ≥ 60` with local refinement along the FSP track.

---

## Results

Run the script and record your values here. The console and
`fsam_data.csv` contain all of them.

### Cold spray stage

| Quantity | Value |
|---|---|
| Mean coating thickness on FSP track | |
| Coating thickness at plate centre (0, 0) | |
| Peak temperature during spraying | |
| Mean as-sprayed porosity | |

### FSP stage

| Quantity | Value |
|---|---|
| Peak temperature | |
| Temperature under shoulder (traverse, steady) | |
| Mean coating porosity after FSP | |
| Approximate axial force Fz | |
| Approximate torque | |

### What to look for

| Field | Expected behaviour |
|---|---|
| `CoatThick` | Ripple of period ≈ hatch across y, flat plateau along x, ramps at the pass ends |
| `TempK` (CS) | Mild, moving hot spot under the nozzle (cold spray is a low-heat process) |
| `TempK` (FSP) | Hot zone around the tool, elongated behind it during traverse |
| `AlFrac` | 1 in coating, 0 in Cu; intermediate values in the stir zone across the interface |
| `Porosity` | `poro0` in the coating, dropping to near zero in the stirred track |
| `PlasticStrain` | Concentrated under the shoulder and around the pin, left behind as a track |

---

## Output Files: `D:\freefem++\cs_fsp_3d\`

| File | Description |
|---|---|
| `fsam.pvd` | ParaView collection for all three stages; open this for the animation |
| `fsam_0000.vtu` ... `fsam_0068.vtu` | Workpiece frames (substrate + coating + void) |
| `tool_0000.vtu` ... | FSP tool (shoulder + pin), parked above the plate during CS and dwell |
| `nozzle_0000.vtu` ... | CS nozzle body + spray plume, parked above the plate during FSP |
| `fsam_data.csv` | Time history of process variables |

### VTU Data Fields

| Field | Description |
|---|---|
| `Vx`, `Vy`, `Vz` | Velocity components (m/s), non-zero only during FSP |
| `TempK` | Temperature (K) |
| `AlFrac` | Al fraction of the solid: 1 Al, 0 Cu |
| `SolidFrac` | Solid fraction: 1 metal, 0 void |
| `Porosity` | Local porosity |
| `PlasticStrain` | Accumulated equivalent plastic strain |
| `CoatThick` | As-sprayed coating thickness (m) |
| `EtaEff` | Effective viscosity (Pa·s) |

### CSV columns

```
fsam_data.csv : time_s, phase, xNoz_mm, yNoz_mm, hCoat_centre_mm, xTool_mm, ztip_mm,
                Tmax_K, Tprobe_K, Fz_N, Torque_Nm, poroMean

phase         : INIT | COLD_SPRAY | DWELL | FSP_PLUNGE | FSP_TRAVERSE | FSP_LIFT
Tprobe_K      : under the nozzle (CS), at the plate centre (dwell),
                under the shoulder at r = Rsh/2 (FSP)
```

---

## Repository Structure

```
cs_fsp_hybrid_3d.edp                   # Main FreeFEM++ simulation script
README.md                              # This file

D:\freefem++\cs_fsp_3d\
├── fsam.pvd                           # Master animation (open in ParaView)
├── fsam_0000.vtu
├── ...
├── fsam_0068.vtu
├── tool_0000.vtu ... tool_0068.vtu
├── nozzle_0000.vtu ... nozzle_0068.vtu
└── fsam_data.csv
```

---

## How to Run

### Requirements
* FreeFEM++ 4.x with the `iovtk` plugin: https://freefem.org
* ParaView: https://www.paraview.org
* Several GB of free RAM (3D UMFPACK factorisation of the Stokes system)

### Step 1: Check the output path

`outdir` is hardcoded to `D:\freefem++\cs_fsp_3d\` and the folder is
created with `system("mkdir ...")`. On Linux or macOS, edit `outdir` and
the `system` line to a writable path.

### Step 2: Run the simulation

```bash
FreeFem++ -nw cs_fsp_hybrid_3d.edp
```

The script will:
1. Build the substrate + void mesh and relabel its faces
2. Raster the cold spray nozzle and grow the coating (thermal only)
3. Let the part cool for `tDwell`
4. Measure the coating thickness on the FSP path and set the pin length
5. Run FSP plunge, traverse and lift with coupled flow, heat and material transport
6. Write 69 frames plus the time history

For the cold spray stage only (fast, no Stokes solves), set
`NstepFSP` to a small value or stop the script after the dwell loop.

### Step 3: Open in ParaView

1. `File > Open`, select `fsam.pvd`, click Apply
2. `Filters > Threshold`, Scalars = `SolidFrac`, lower = 0.5, click Apply
3. Colour by `TempK`, `AlFrac` or `Porosity`
4. Click "Rescale to Data Range Over All Timesteps", then press Play

---

## Visualisation in ParaView

### Option A: Watch the coating grow

```
Open fsam.pvd
Threshold  SolidFrac > 0.5
Colour by CoatThick (or TempK)
Show part 2 (nozzle) as Surface
Play through the COLD_SPRAY frames: layers build up pass by pass
```

### Option B: Watch the stir zone form

```
Open fsam.pvd
Threshold  SolidFrac > 0.5
Filters > Clip, plane normal (0,1,0) through y = 0 to see the cross-section
Colour by AlFrac (mixing) or Porosity (densification)
Show part 1 (tool) as Surface
```

### Option C: Temperature field around the tool

```
Open fsam.pvd, go to an FSP_TRAVERSE frame
Threshold  SolidFrac > 0.5
Colour by TempK, add Contour at 500, 600, 700, 800 K
Add Glyph on (Vx,Vy,Vz) to see the material flow around the pin
```

---

## Key Numerical Insights

### The coating is built on a fixed mesh
No elements are added during cold spray. The void layer is part of the
mesh from the start, and `SolidFrac` marks which part of it is metal. The
deposit surface is therefore resolved to one element (0.5 mm) in z.

### Deposited mass is time-step independent
The `erf` segment average means a coarse `NsubCS` changes the timing of
heat input but not the amount of material deposited.

### Cold spray heating is small by design
With the defaults, particle impact contributes about 44 W against about
1.4 kW from the FSP shoulder. This is consistent with cold spray being a
solid-state, low-heat process; the dwell mainly matters if `hJet0` or
`TgasJet` are raised.

### FSP is quasi-steady, not rotation-resolved
One FSP time step (≈ 0.59 s) spans about 12 tool revolutions. The flow
field is meaningful as a quasi-steady stirring pattern, but `AlFrac`
mixing should be read qualitatively.

### Pin displacement is not modelled
Penalisation forces the material inside the tool volume to rotate with
the tool; it is not pushed out of the pin's volume. Mass conservation is
approximate during the plunge.

### Porosity closure is strain-driven
`Cden = 1` closes about 95 % of the porosity after an accumulated strain
of 3. Temperature dependence of densification is not included.

### Interface is perfectly bonded
The coating/substrate interface has no separate bonding, oxide or
contact model; it is represented only by the jump in `AlFrac`.

### Output buffering
FreeFEM++ buffers `ofstream` output, so CSV files may appear to stall
during a run and complete only when the script finishes. This is normal.

---

## Common Errors and Fixes

| Error / symptom | Cause | Fix |
|---|---|---|
| `mkdir` fails / output folder missing | Windows path on a non-Windows system | Edit `outdir` and the `system` line, or create the folder by hand |
| `load "iovtk"` fails | VTK plugin not installed or not on the load path | Install the full FreeFEM++ package |
| Error at `change(..., flabel=...)` | Older FreeFEM++ without 3D face relabelling | Print `labels(Th)` and map the `cube()` default labels instead |
| Macro expansion errors | A `//` comment inside a macro body ends the macro early | Keep comments outside `macro ... // EOM` blocks |
| Fewer passes than expected | `int(Wband/hatch)` rounds down on floating-point error | Use `int(Wband/hatch + 0.5) + 1` |
| Out of memory in the FSP stage | 3D direct factorisation of ~73k unknowns every step | Reduce `Nx`/`Ny`, or switch to an iterative / parallel solver |
| Coating not visible in ParaView | Void cells shown as well | Apply Threshold on `SolidFrac > 0.5` |
| Pin heat or flow flickers during traverse | Pin narrower than the x element size | Refine `Nx` or keep the inflated `Rpen` |
| Colours flicker during playback | Colour range rescaled per frame | "Rescale to Data Range Over All Timesteps" |
| ParaView cannot find frames | `.vtu` files moved away from `fsam.pvd` | Keep all files in the same folder |

---

## Extending the Model

| Extension | What to change |
|---|---|
| Al–Cu composite powder | Set `fAlPow` below 1 (e.g. 0.7 for 30 % Cu) |
| Thicker or thinner coating | Change `Nlayer`, `mdotPow` or `hatch` |
| Cold FSP start | Increase `tDwell` |
| Hotter cold spray | Raise `TgasJet`, `hJet0` or `Tpart` |
| FSP confined to the coating | Set `dPinSub` negative (pin tip stays above the interface) |
| Multi-pass FSP | Repeat Stage C with a lateral offset of the tool path |
| Temperature-dependent densification | Multiply `Cden` by an Arrhenius factor in `Temp` |
| Intermetallic growth at the interface | Add a diffusion–reaction field driven by `TempK` and `AlFrac` |
| Converged mixing results | Refine along y = 0 (`Nx ≥ 60`) and shorten the FSP time step |
| Separate tool material | Give the penalised region steel thermal properties |

---

## Citation

Mishra, A. (2026). *Hybrid Cold Spray + Friction Stir Processing (CS +
FSP) of an Al Coating on Cu: A FreeFEM++ 3D Eulerian (CEL-type)
Thermo-Mechanical Study* [Computer software]. Zenodo.
https://doi.org/10.5281/zenodo.23266913

```bibtex
@software{mishra2026csfsp,
  author       = {Mishra, Akshansh},
  title        = {Hybrid Cold Spray + Friction Stir Processing (CS + FSP)
                   of an Al Coating on Cu: A FreeFEM++ 3D Eulerian
                   (CEL-type) Thermo-Mechanical Study},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23266913},
  url          = {https://doi.org/10.5281/zenodo.23266913}
}
```

---

## License

<p>
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
<img alt="Creative Commons Licence" style="border-width:0"
  src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png"/>
</a>
<br/>
This work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
Creative Commons Attribution NonCommercial 4.0 International License</a>.
</p>

You are free to:
* **Share**: copy and redistribute in any medium or format
* **Adapt**: remix, transform, and build upon the material

Under the following terms:
* **Attribution**: give appropriate credit and link to this repository
* **NonCommercial**: not for commercial use without permission

Copyright 2026 Akshansh Mishra.
