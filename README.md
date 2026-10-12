# Tensor-Geometry-Activated Spectral Energy Transfer in 2-D Flows — Data & Code

Companion data and code for the paper:

> **Tensor-Geometry-Activated Spectral Energy Transfer at Reynolds Numbers of Order Unity in Two-Dimensional Flows**
> Ziyue Yu, Xinyu Si, Lei Fang
> *Journal of Fluid Mechanics*, JFM-2026-0511 (accepted 2026-08-26)

Electromagnetically driven thin-layer 2-D flow experiments. A 4 × 4 rod array
(spacing `L_i = 1.5 cm`) perturbs a large-scale background flow (cellular or
shear forcing). The stress / strain-rate tensor geometry activates a
turbulent-like cross-scale energy flux at Reynolds numbers of order unity.

This repository is a **compact demo**: representative velocity fields, the
matching raw camera frames, and the routines that compute the energy flux and
the stress–strain alignment angle. It is meant to let a reader run the method
end-to-end, not to reproduce every figure. Full-length records (1501 / 1200
frames per set, 40+ sets per experiment day) and complete raw image sequences
are kept offline by the authors and are not part of this release. Researchers
interested in collaborating on this data are welcome to contact the authors.

---

## Quick start

```matlab
% from the repository root
demo_energy_flux
```

`demo_energy_flux.m` loads the four sets of one experiment day (cellular by
default; set `cfg = 'shear'` for the other) and produces

1. the filter-space total energy flux `Π(L)` for `L = 1…12 cm`, one curve per
   configuration (pure perturbation, inverse, forward, neutral);
2. the PDF of the stress–strain alignment angle `θ^(L)` at `L = 7 cm`, one
   curve per configuration;
3. for inverse, forward and neutral, one frame of the underlying fields: the
   local flux `Π^(L)(x,y)`, the largest eigenvalue and extensional eigenvector
   of the filtered rate of strain, and the largest eigenvalue and eigenvector
   of the sub-grid stress (eigenvectors drawn as bars, every 6th grid point).

It runs in under a minute. The formulas behind each step, and the code lines
that implement them, are listed in `Formula.md`.

## Repository layout

```
demo_energy_flux.m       Runnable example (flux vs L, PDF of theta, flux / strain / stress fields)
Formula.md               Every formula of the analysis and where it is implemented
1_experiment_records/    Lab logs (xlsx) for the two experiment days
2_data/                  Velocity fields, 8 configurations, 21 frames each
3_photos/rawdata/        Raw camera frames, one per velocity frame (PNG)
3_photos/processed/      Same frames after PTV pre-processing (background & rods removed)
4_functions/             The three MATLAB routines used by the demo
5_figures/               Figures of the published paper (main text + SI)
6_videos/                Particle-trajectory movies of the flow configurations (mp4)
```

## 1. Experiment records

| File | Date | Background forcing |
|---|---|---|
| `250328_cellular.xlsx` | 2025-03-28 | Cellular (array of large eddies) |
| `250530_shear.xlsx` | 2025-05-30 | Shear |

Per-set logs: current, frame rate, frame count, purpose. Camera: 1600 × 1600 px,
field of view 20 cm × 20 cm (80 px/cm), 60 frames/s.

## 2. Data (`2_data/`)

Eulerian velocity fields obtained by particle tracking velocimetry (PTV) and
cubic interpolation of the particle velocities onto a uniform grid covering the
20 cm × 20 cm field of view. Each file keeps the native grid of that set and
**21 frames, every 10th camera frame from 101 to 301** (3.3 s of flow).

| Configuration | File | Original set | Grid |
|---|---|---|---|
| Cellular, pure perturbation (rod array only) | `cellular/Cellular_PurePerturbation_set21.mat` | 250328 set21 | 127 × 127 |
| Cellular, inverse (rod axis at θ ≈ 0) | `cellular/Cellular_Inverse_set27.mat` | 250328 set27 | 135 × 135 |
| Cellular, forward (θ ≈ π/2) | `cellular/Cellular_Forward_set32.mat` | 250328 set32 | 122 × 122 |
| Cellular, neutral (θ ≈ π/4) | `cellular/Cellular_Neutral_set39.mat` | 250328 set39 | 134 × 134 |
| Shear, pure perturbation (rod array only) | `shear/Shear_PurePerturbation_set24.mat` | 250530 set24 | 146 × 146 |
| Shear, inverse | `shear/Shear_Inverse_set30.mat` | 250530 set30 | 154 × 154 |
| Shear, forward | `shear/Shear_Forward_set39.mat` | 250530 set39 | 155 × 155 |
| Shear, neutral | `shear/Shear_Neutral_set48.mat` | 250530 set48 | 135 × 135 |

Variables in each `.mat` (MATLAB v7, loadable with `scipy.io.loadmat`):

| Variable | Size | Meaning |
|---|---|---|
| `px`, `py` | n × n | Grid coordinates (cm) |
| `pvx`, `pvy` | n × n × 21 | Velocity components (cm/s) |
| `frames` | 1 × 21 | Camera frame number of each velocity frame (`101:10:301`) |
| `fps` | scalar | 60 frames/s |
| `note` | char | Provenance note |

`pvx(:,:,k)` is the velocity field at camera frame `frames(k)` and corresponds
to `3_photos/rawdata/<config>/<set>/Image<frames(k)>.png` and
`3_photos/processed/<config>/<set>/Image<frames(k)>.png`.

**Reynolds number.** The flow is dominated by bottom friction rather than
viscosity, so the paper uses the friction-based Reynolds number

    Re = u_rms / (L_i α)

with `u_rms` the rms velocity over the analysis subdomain (dashed box in
figure 1 of the paper), `L_i = 1.5 cm` the rod spacing and `α` the linear
friction coefficient measured from kinetic-energy decay
(`α = 0.12 s⁻¹` cellular, `α = 0.055 s⁻¹` shear; see the Supplementary
Material). The conventional viscous Reynolds number is

    Re_ν = u_rms L_i / ν,   ν = 0.0125 cm²/s

so that `Re_ν ≈ 21.6 Re` (cellular) and `≈ 9.9 Re` (shear). Both can be
evaluated from the archived fields with these formulas.

Quick look:

```matlab
S = load('2_data/cellular/Cellular_Inverse_set27.mat');
quiver(S.px, S.py, S.pvx(:,:,1), S.pvy(:,:,1)); axis equal tight
```

## 3. Photos (`3_photos/`)

`rawdata/`: raw Basler camera frames (1600 × 1600, 8-bit), converted **losslessly** from
tif to PNG. For each of the 8 configurations: `Image0101.png … Image0301.png`
(every 10th frame, one-to-one with the velocity frames in `2_data/`) plus the
`background.png` reference image used for background subtraction.

`processed/`: the same frames after the PTV pre-processing step, i.e. the
images actually handed to the particle tracker. The background image is
subtracted and the rod array (with its shadows) is masked by replacing every
connected bright region larger than 100 px with the background. Each folder
also holds the `background.png` of the pre-processed sequence. Processing
parameters follow the original pipeline (binarisation threshold 0.01, minimum
object area 100 px).

## 4. Functions (`4_functions/`)

| File | Purpose |
|---|---|
| `fst_geo_efficiency_2D.m` | Filter-space total energy flux `Π(L) = −τ_ij S̄_ij` on a 2-D velocity field for cut-off length `R` (cm). Optional outputs return the Leonard / cross / sub-grid contributions separately. |
| `Cal_eigvector_eigvalue_strainandstress.m` | Eigenvalues and eigenvector angle of the filtered rate-of-strain tensor (`type = 1`) or of the sub-grid stress tensor (`type = 0`); the difference of the two angles is the alignment angle θ^(L) of the paper. |
| `sharpfilt2.m` | Sharp spectral low-pass filter with Gaussian conditioning of the cut-off (X. Si); the filter kernel used by both routines. |

Usage pattern (see `demo_energy_flux.m`):

```matlab
addpath('4_functions');
Q = fst_geo_efficiency_2D(S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), L, 0.5);   % L in cm
Pi_L = mean(Q, 'all');
```

## 5. Figures (`5_figures/`)

- `main/` — the 4 main-text figures (experimental setup; flux enhancement;
  PDF of θ; enhancement vs Re).
- `SI/` — the 13 Supplementary Material figures (setup details, perturbation
  field, α statistics, flux history, 3-D checks, scale locality,
  linear-superposition controls Π(L), filter / resolution / domain robustness).

## 6. Videos (`6_videos/`)

Particle-trajectory movies (particle streaks over a few frames, 20 fps
playback, 3.7–4.8 s) of the main flow configurations, made from the same
records as `2_data/`. Files ending in `_R` are a rotated and cropped close-up
of the rod-array region (1200 × 1200 px); the others show the full field of
view (1500 × 1500 px).

| File | Configuration | Set |
|---|---|---|
| `Set_8_cellular.mp4` | Cellular background flow, no rod array | 250328 set8 |
| `Set_21_cellular_perb.mp4` | Cellular, pure perturbation (rod array only) | 250328 set21 |
| `Set_27_cellular_back_perb.mp4` | Cellular, inverse (θ ≈ 0) | 250328 set27 |
| `Set_32_cellular_for_perb.mp4` | Cellular, forward (θ ≈ π/2) | 250328 set32 |
| `Set_13_shear.mp4`, `Set_13_shear_R.mp4` | Shear background flow, no rod array | 250530 set13 |
| `Set_24_shear_perb.mp4`, `Set_24_shear_perb_R.mp4` | Shear, pure perturbation (rod array only) | 250530 set24 |
| `Set_30_shear_back_perb.mp4`, `Set_30_shear_back_perb_R.mp4` | Shear, inverse | 250530 set30 |
| `Set_39_shear_for_perb.mp4`, `Set_39_shear_for_perb_R.mp4` | Shear, forward | 250530 set39 |

---

*Demo archive assembled 2026-10-12 from the full archive kept by the authors.*
