# Formulas and their implementation

This note lists every formula used in the paper's analysis chain and shows where
each one is implemented. Line numbers refer to the files in `4_functions/` and
to `demo_energy_flux.m`. Steps of the original pipeline that are *not* shipped
in this demo (particle tracking, gridding, Reynolds number, friction
coefficient) are documented here with the exact formulas and parameters so
that they can be reproduced.

Notation: $u_i$ is the measured velocity (cm/s), $x_i$ the position (cm),
superscript $(L)$ a low-pass filter with cut-off length $L$, repeated indices
are summed, $\langle\cdot\rangle$ is a spatial average, $\overline{\cdot}$ a
time average. MATLAB arrays: `px, py` grid (cm), `pvx, pvy` velocity (cm/s).

---

## 1. From images to velocity fields (original pipeline, not in demo)

**Image pre-processing** (`3_photos/processed/`). With the raw frame $I$ and
the background $B$ (both 8-bit, 1600 × 1600):

$$I' = I - B,\qquad
M = \mathrm{bwareaopen}\big(\mathrm{imbinarize}(\mathrm{imadjust}(I'),\,0.01),\,100\big),\qquad
I_{\rm out} = \big(I\,[M \leftarrow B]\big) - B .$$

Every connected bright region larger than 100 px (the rods and their shadows)
is replaced by the background before the subtraction. Parameters:
`threshold = 0.01`, `arealim = 100` px.

**PTV.** Particle positions and velocities (px, px/frame) come from the tracker
of Ouellette et al. (2006). Particles with a speed above
$10\times\mathrm{IQR}(\text{speed})$ in a frame are discarded.

**Gridding.** Particle velocities are interpolated onto a uniform grid with
MATLAB `griddata(...,'cubic')`; grid points with no data are set to 0.

**Units.** With `ppcm = 80` px/cm and `fps = 60` frames/s:

$$x\,[\mathrm{cm}] = \frac{x\,[\mathrm{px}]}{80},\qquad
u\,[\mathrm{cm/s}] = \frac{u\,[\mathrm{px/frame}]}{80}\times 60 .$$

The files in `2_data/` are already in cm and cm/s on the native grid of each
set (grid spacing ≈ 0.13–0.15 cm).

---

## 2. Low-pass filter

Paper (SI, Eq. *Filtering*):

$$u_i^{(L)}(\mathbf{x}) = \int G^{(L)}(\mathbf{r})\,u_i(\mathbf{x}-\mathbf{r})\,d\mathbf{r},$$

with $G^{(L)}$ a sharp spectral cut-off at wavelength $L$, conditioned by a
Gaussian window to suppress ringing.

Implementation: `sharpfilt2.m` (`y = sharpfilt2(x, Wn, sigma, lowpass, padding)`).

| Quantity | Formula | Code |
|---|---|---|
| grid spacing | $d_0 = (x_{\rm end}-x_1)/N$ | `fst_geo_efficiency_2D.m` line 4: `d0=(X(end)-X(1))/length(X)` |
| cut-off (fraction of Nyquist) | $W_n = 2d_0/L$, i.e. $k_c = W_n\,\pi/d_0 = 2\pi/L$ | `fst_geo_efficiency_2D.m` line 5: `Wp=2*d0/R` |
| ideal low-pass | $\hat H(\mathbf k)=1$ for $\lvert\mathbf k\rvert\le k_c$, 0 otherwise | `sharpfilt2.m` line 51: `Hd(freqamp>Wn)=0` |
| Gaussian conditioning | $h(\mathbf r)\leftarrow h(\mathbf r)\,e^{-\lvert\mathbf k\rvert^2/2\sigma^2}$, $\sigma=0.5$ | `sharpfilt2.m` lines 56 and 72: `win=exp(-(freqamp.^2)/2/sigma^2)`, `h=h.*win` |
| convolution, zero padding | $u^{(L)} = h * u$ | `sharpfilt2.m` lines 98 and 100: `padImage(A,hsize,padding)`, `ifftn(fftn(A).*fftn(h))` |

Call used everywhere: `sharpfilt2(field, Wn, 0.5, 1, 0)` (low-pass, zero
padding). In `Cal_eigvector_eigvalue_strainandstress.m` the same cut-off is
built from the arguments `L` and `d0`: line 5 `Wn=2*d0/L`.

The analysis scale of the paper is $L = 0.7\,W = 7$ cm with $W = 10$ cm (half
the measurement domain). Resolution criterion: $L/d_0 \ge 6$ is satisfied for
all $L \ge 1$ cm.

---

## 3. Sub-grid stress and rate of strain

Paper (main text §2, SI *Energy flux calculation and decomposition*):

$$\tau_{ij}^{(L)} = (u_iu_j)^{(L)} - u_i^{(L)}u_j^{(L)}
= \underbrace{(u_i^{(L)}u_j^{(L)})^{(L)} - u_i^{(L)}u_j^{(L)}}_{\tau_{L}}
+ \underbrace{(u_i^{(L)}u_j')^{(L)} + (u_j^{(L)}u_i')^{(L)}}_{\tau_{C}}
+ \underbrace{(u_i'u_j')^{(L)}}_{\tau_{S}},
\qquad u_i' = u_i - u_i^{(L)} .$$

$$s_{ij}^{(L)} = \tfrac12\Big(\frac{\partial u_i^{(L)}}{\partial x_j} + \frac{\partial u_j^{(L)}}{\partial x_i}\Big).$$

Implementation in `fst_geo_efficiency_2D.m`:

| Quantity | Code (line) |
|---|---|
| $u^{(L)}, v^{(L)}$ | 9–10 `u=sharpfilt2(velocitydatax,Wp,sigma,lowpass,0)`, same for `v` |
| $u', v'$ | 13–14 `du = velocitydatax - u`, `dv = velocitydatay - v` |
| $\tau_{S,11},\tau_{S,12},\tau_{S,22}$ | 17–20 `t11=sharpfilt2(du.^2,…)`, `t12=sharpfilt2(du.*dv,…)`, `t22=sharpfilt2(dv.^2,…)` |
| $\partial_x u^{(L)}, \partial_y u^{(L)}, \dots$ | 28–29 `[u_x,u_y] = gradient(u,dx,dy)`, `[v_x,v_y] = gradient(v,dx,dy)` |
| $s_{11}, s_{12}=s_{21}, s_{22}$ | 31–34 `s11=u_x`, `s12=1/2*(u_y+v_x)`, `s22=v_y` |
| $\tau_L$ (Leonard) | 52–55 `u2_bar - u.^2`, `uv_bar - u.*v`, `v2_bar - v.^2` |
| $\tau_C$ (cross) | 58–61 `sharpfilt2(u.*du + du.*u,…)`, … |

Note `gradient(F, dx, dy)` returns $(\partial_x F, \partial_y F)$ for `F(y,x)`
arrays, so `u_y` is $\partial u/\partial y$.

---

## 4. Energy flux across scale $L$

Paper (main text, SI Eq. *MotionEnergy*):

$$\Pi^{(L)} = -\tau_{ij}^{(L)}\,s_{ij}^{(L)},
\qquad \Pi^{(L)} = \Pi_L^{(L)} + \Pi_C^{(L)} + \Pi_S^{(L)},\quad
\Pi_X^{(L)} = -\tau_{X,ij}^{(L)}\,s_{ij}^{(L)} .$$

$\Pi^{(L)}<0$: inverse transfer (towards scales larger than $L$);
$\Pi^{(L)}>0$: forward transfer. Units cm$^2$ s$^{-3}$. The paper reports
$\Pi_S^{(L)}$ (the sub-grid-stress part) and calls it $\Pi^{(L)}$; this is also
what the demo calls the *total energy flux*.

Implementation, `fst_geo_efficiency_2D.m`:

| Output | Code (line) |
|---|---|
| `Q` $=\Pi_S^{(L)}$ | 37 `Q = -(s11.*t11 + s12.*t12 + s21.*t21 + s22.*t22)` |
| `Q_L` $=\Pi_L^{(L)}$ | 70 |
| `Q_C` $=\Pi_C^{(L)}$ | 71 |
| `Q_S` $=\Pi_S^{(L)}$ (identical to `Q`) | 72 |

Spatial and temporal average used in the figures:

$$\overline{\langle\Pi^{(L)}\rangle} = \frac{1}{N_t}\sum_t \frac{1}{\lvert A\rvert}\sum_{\mathbf x\in A}\Pi^{(L)}(\mathbf x,t),$$

with $A$ the analysis subdomain (dashed box in the figures) and $t$ running
over the frames of the record. `demo_energy_flux.m` lines 27–36 do the same
with $A$ = full field and the 21 archived frames.

The original pipeline evaluated `fst_geo_efficiency_2D` in pixel units and
converted afterwards: $\Pi[\mathrm{cm^2 s^{-3}}] = \Pi[\mathrm{px^2\,frame^{-3}}]\times 60^3/80^2$.
The demo passes cm and cm/s directly, which is equivalent.

In the figures $\Pi^{(L)}$ is normalised by the frictional dissipation
$\alpha\,\overline{\langle\lvert\mathbf u\rvert^2\rangle}$ (see §7) and lengths by $W=10$ cm.

---

## 5. Tensor geometry: eigenvalues, eigenvectors, alignment angle

Paper (main text Eq. 2.1, SI *Tensor geometry*). With $\sigma^{(L)}$ the
extensional eigenvalue of $s_{ij}^{(L)}$, $\gamma^{(L)}$ that of the deviatoric
stress $\tau_{ij,{\rm Dev}}^{(L)} = \tau_{ij}^{(L)} - \tfrac12\tau_{kk}^{(L)}\delta_{ij}$,
and $\theta^{(L)}$ the angle between the two extensional eigenvectors
$\hat\sigma^{(L)}, \hat\gamma^{(L)}$:

$$\Pi^{(L)} = -2\,\sigma^{(L)}\gamma^{(L)}\cos\big(2\theta^{(L)}\big),
\qquad \eta^{(L)} = \cos\big(2\theta^{(L)}\big)\ \ (\text{efficiency}).$$

$\theta^{(L)}<\pi/4$: inverse flux, $\theta^{(L)}>\pi/4$: forward flux,
$\theta^{(L)}=\pi/4$: no flux.

Implementation, `Cal_eigvector_eigvalue_strainandstress.m`
(`[x,y,ex,ey,lambda,angle] = Cal_eigvector_eigvalue_strainandstress(px,py,pvx,pvy,type,L,d0,filter)`):

| Quantity | Code (line) |
|---|---|
| filtered velocity | 7–8 `pvx=sharpfilt2(pvxs,Wn,0.5,1,0)`, same for `pvy` |
| $\partial_j u_i^{(L)}$ | 18–19 `[U2x,U2y] = gradient(pvx,…)`, `[V2x,V2y] = gradient(pvy,…)` |
| $\tau_{S,ij}$ | 27–30 `subt11=sharpfilt2((pvxs-pvx).^2,…)`, … |
| $s_{ij}$ matrix (`type=1`) | 36–37 `strain=[U2x, 0.5*(U2y+V2x); 0.5*(U2y+V2x), V2y]` |
| $\tau_{S,ij}$ matrix (`type=0`) | 47–48 `stress=[subt11,subt12; subt21,subt22]` |
| eigen-decomposition | 52 `[V,D] = eig(strain)` |
| extensional eigenvalue ($\sigma$, or $\lambda_{\max}(\tau_S)$) | 54 `eig_val(k,l) = max(diag(D))` |
| extensional eigenvector $(\hat e_x,\hat e_y)$, sign fixed so $\hat e_x\ge0$ | 64–67 |
| eigenvector angle (deg) | 71 `eigvec_angle(k,l) = atand(eigvec_y/eigvec_x)` |

Alignment angle and efficiency, `demo_energy_flux.m` lines 52–54 (same as the
original `d_Main_analysis_global.m`):

```matlab
th = abs(ang_stress - ang_strain);   % degrees, both angles in (-90, 90]
th(th > 90) = 180 - th(th > 90);     % eigen-directions are axial: fold to [0, 90]
theta = th/180;                      % theta^(L) in units of pi, in [0, 1/2]
eta   = cosd(2*th);                  % efficiency cos(2 theta^(L))
```

Remark. The routine diagonalises the full $\tau_S$, not its deviatoric part.
The eigenvectors (hence $\theta^{(L)}$ and $\eta^{(L)}$) are identical for both,
because subtracting $\tfrac12\tau_{kk}\delta_{ij}$ does not change eigenvectors;
the deviatoric eigenvalue is $\gamma^{(L)} = \lambda_{\max}(\tau_S) - \tfrac12\tau_{kk}$.

---

## 6. Reynolds numbers

Paper (main text §2, SI *Transformation between* $Re_{\rm traditional}$ *and* $Re_\alpha$):

$$\mathrm{Re} \equiv \mathrm{Re}_\alpha = \frac{u_{\rm rms}}{L_i\,\alpha},\qquad
\mathrm{Re}_\nu = \frac{u_{\rm rms}\,L_i}{\nu},\qquad
\frac{\mathrm{Re}_\nu}{\mathrm{Re}_\alpha} = \frac{L_i^2\alpha}{\nu}
= \begin{cases}21.6 & \text{cellular}\\ 9.9 & \text{shear}\end{cases}$$

with $L_i = 1.5$ cm (rod spacing), $\nu = 0.0125$ cm$^2$/s,
$\alpha = 0.12$ s$^{-1}$ (cellular) and $0.055$ s$^{-1}$ (shear), and

$$u_{\rm rms} = \sqrt{\overline{\langle u_x^2 + u_y^2\rangle_A}}$$

evaluated over the analysis subdomain $A$ and the frames of the record.

Not computed in the demo; from an archived file:

```matlab
S = load('2_data/cellular/Cellular_Inverse_set27.mat');
L_i = 1.5; alpha = 0.12; nu = 0.0125;               % cellular; shear: alpha = 0.055
u_rms    = sqrt(mean(S.pvx.^2 + S.pvy.^2, 'all'));   % full field, archived frames
Re_alpha = u_rms/(L_i*alpha);
Re_nu    = u_rms*L_i/nu;
```

In the original pipeline the subdomain is a polygon mask (`inpolygon` on the
four corners stored per set) and frames whose $\langle\lvert\mathbf u\rvert^2\rangle$
exceeds 3× the record mean (tracking glitches) are replaced by the mean before
averaging.

---

## 7. Friction coefficient and dissipation

Linear bottom friction enters the momentum equation as $-\alpha\,u_i$. After
the forcing is switched off, $u \propto e^{-\alpha t}$, hence

$$\langle\lvert\mathbf u\rvert^2\rangle(t) = \langle\lvert\mathbf u\rvert^2\rangle_0\,e^{-2\alpha t}.$$

Original pipeline (`d_Main_analysis_global.m`, decay sets): the time series of
$\langle\lvert\mathbf u\rvert^2\rangle$ over the analysis subdomain is fitted with
MATLAB `fit(t, E, 'exp1')` ($E = a\,e^{bt}$) and $\alpha = -b/2$. The values
adopted in the paper are $\alpha = 0.12$ s$^{-1}$ (cellular, record 250328)
and $0.055$ s$^{-1}$ (shear, record 250530).

Frictional dissipation rate used to normalise $\Pi^{(L)}$ in the figures:

$$\varepsilon_\alpha = \alpha\,\langle u_x^2+u_y^2\rangle$$

(original `zFunction_cal_dissipation(..., "friction", alpha)`:
`Dissipation = friction*(pvx.^2+pvy.^2)`).

---

## 8. Summary of parameters

| Symbol | Value | Where |
|---|---|---|
| field of view | 20 cm × 20 cm, 1600 × 1600 px, 80 px/cm | README §1 |
| frame rate | 60 frames/s | `fps` in each `.mat` |
| $L_i$ | 1.5 cm | §6 |
| $\alpha$ | 0.12 s$^{-1}$ (cellular), 0.055 s$^{-1}$ (shear) | §6, §7 |
| $\nu$ | 0.0125 cm$^2$/s | §6 |
| filter conditioning $\sigma$ | 0.5 | `fst_geo_efficiency_2D(...,R,0.5)`, `Cal_eigvector...` line 7 |
| analysis scale | $L = 7$ cm ($0.7\,W$, $W = 10$ cm) | `demo_energy_flux.m` line 45 |
| flux scan | $L = 1,\dots,12$ cm | `demo_energy_flux.m` line 27 |
| image pre-processing | threshold 0.01, min. area 100 px | §1 |
| PTV outlier cut | speed > 10 × IQR | §1 |
