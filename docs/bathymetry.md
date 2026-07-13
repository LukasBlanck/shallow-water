# Numerical scheme for `Reconstruction::MUSCL` + `Riemann::HLL` + well-balanced hydrostatic reconstruction + dry threshold

---

## PDE

We consider the 2D shallow water equations with bottom topography:

$$
U=
\begin{bmatrix}
h\\
hu\\
hv
\end{bmatrix},
\qquad
b=b(x,y),
\qquad
\eta = h+b.
$$

The PDE is

$$
\partial_t U + \partial_x F(U) + \partial_y G(U) = S(U,b),
$$

with fluxes

$$
F(U)=
\begin{bmatrix}
hu\\
hu^2+\frac12 g h^2\\
huv
\end{bmatrix},
\qquad
G(U)=
\begin{bmatrix}
hv\\
huv\\
hv^2+\frac12 g h^2
\end{bmatrix},
$$

and source term

$$
S(U,b)=
\begin{bmatrix}
0\\
-gh\,\partial_x b\\
-gh\,\partial_y b
\end{bmatrix}.
$$

Here:
- $h(x,y,t)$ .... water height
- $u(x,y,t)$ .... x-velocity
- $v(x,y,t)$ .... y-velocity
- $b(x,y)$ .... bottom topography
- $\eta=h+b$ .... free surface

---

## Why hydrostatic reconstruction?

For non-flat bottom, an important steady state is the **lake at rest** state:

$$
u=v=0,
\qquad
h+b=\text{const}.
$$

A naive discretization of the flux gradient plus source term usually does **not** preserve this exactly. Then one gets artificial velocities although the exact solution should stay at rest.

Hydrostatic reconstruction is designed to:
- preserve the lake-at-rest equilibrium much better,
- avoid negative reconstructed water heights,
- behave robustly near wet/dry fronts.

---

## Grid

Domain

$$
\Omega = [0,L_x]\times[0,L_y].
$$

$$
\Delta x=\frac{L_x}{N_x}, \qquad \Delta y=\frac{L_y}{N_y}.
$$

##### Cell centers

$$
x_i=\left(i+\frac12\right)\Delta x,
\qquad
i=0,\dots,N_x-1,
$$

$$
y_j=\left(j+\frac12\right)\Delta y,
\qquad
j=0,\dots,N_y-1.
$$

##### Control volumes

$$
C_{ij}=
[x_{i-\frac12},x_{i+\frac12}]
\times
[y_{j-\frac12},y_{j+\frac12}],
$$

with area

$$
|C_{ij}|=\Delta x\,\Delta y.
$$

---

## Cell averages

FV stores cell averages:

$$
\bar U_{ij}(t)=
\frac{1}{\Delta x\,\Delta y}
\int_{C_{ij}} U(x,y,t)\,dx\,dy=
\begin{bmatrix}
\bar h_{ij}\\
\overline{hu}_{ij}\\
\overline{hv}_{ij}
\end{bmatrix}.
$$

Also store the cell average of the bottom:

$$
\bar b_{ij}=
\frac{1}{\Delta x\,\Delta y}
\int_{C_{ij}} b(x,y)\,dx\,dy.
$$

Define also the cell-average free surface

$$
\bar\eta_{ij} = \bar h_{ij} + \bar b_{ij}.
$$

---

## Dry threshold

Introduce a small dry tolerance

$$
h_{\mathrm{dry}} > 0.
$$

Typical use:
- if $h < h_{\mathrm{dry}}$, the cell is treated as dry,
- in dry cells set

$$
hu=0,\qquad hv=0.
$$

Also when computing velocities:

$$
u=\frac{hu}{h},\qquad v=\frac{hv}{h}
$$

should only be used if $h \ge h_{\mathrm{dry}}$.

Otherwise use

$$
u=0,\qquad v=0.
$$

This avoids division by tiny water heights.

---

## MUSCL reconstruction

We reconstruct piecewise linear states.

For hydrostatic reconstruction, the important point is that we reconstruct the **free surface**

$$
\eta = h+b
$$

directly, instead of reconstructing $h$ first and only afterwards forming $\eta$.

We therefore reconstruct the variables

$$
W =
\begin{bmatrix}
\eta\\
hu\\
hv
\end{bmatrix}.
$$

---

### Slopes in x-direction

A simple choice is the **minmod** limiter applied componentwise:

$$
\sigma^x_{ij}=
\operatorname{minmod}\!\left(
\bar W_{ij}-\bar W_{i-1,j},
\bar W_{i+1,j}-\bar W_{ij}
\right).
$$

### Slopes in y-direction

$$
\sigma^y_{ij}=
\operatorname{minmod}\!\left(
\bar W_{ij}-\bar W_{i,j-1},
\bar W_{i,j+1}-\bar W_{ij}
\right).
$$

The minmod function is

$$
\operatorname{minmod}(a,b)=
\begin{cases}
\operatorname{sign}(a)\min(|a|,|b|), & ab>0,\\[1mm]
0, & ab\le 0.
\end{cases}
$$

---

### Reconstructed interface values

At an x-interface $\left(i+\frac12,j\right)$:

$$
W_{i+\frac12,j}^{-}=
\bar W_{ij}+\frac12\sigma^x_{ij},
\qquad
W_{i+\frac12,j}^{+}=
\bar W_{i+1,j}-\frac12\sigma^x_{i+1,j}.
$$

At a y-interface $\left(i,j+\frac12\right)$:

$$
W_{i,j+\frac12}^{-}=
\bar W_{ij}+\frac12\sigma^y_{ij},
\qquad
W_{i,j+\frac12}^{+}=
\bar W_{i,j+1}-\frac12\sigma^y_{i,j+1}.
$$

Write these components as

$$
W^- =
\begin{bmatrix}
\eta^-\\
(hu)^-\\
(hv)^-
\end{bmatrix},
\qquad
W^+ =
\begin{bmatrix}
\eta^+\\
(hu)^+\\
(hv)^+
\end{bmatrix}.
$$

---

### Bottom reconstruction

For the bottom topography, I use a **piecewise constant reconstruction**. No bottom slopes are computed.

At an x-interface $\left(i+\frac12,j\right)$, the bottom values on the two sides are

$$
b_{i+\frac12,j}^{-}=\bar b_{ij},
\qquad
b_{i+\frac12,j}^{+}=\bar b_{i+1,j}.
$$

At a y-interface $\left(i,j+\frac12\right)$, use

$$
b_{i,j+\frac12}^{-}=\bar b_{ij},
\qquad
b_{i,j+\frac12}^{+}=\bar b_{i,j+1}.
$$

Thus, each cell contributes its cell-average bottom value directly to its adjacent interfaces.

For the hydrostatic reconstruction, define the common interface bottom as the maximum of the two one-sided values:

$$
b_{i+\frac12,j}^{*}
=
\max\left(
\bar b_{ij},
\bar b_{i+1,j}
\right),
$$

and

$$
b_{i,j+\frac12}^{*}
=
\max\left(
\bar b_{ij},
\bar b_{i,j+1}
\right).
$$

Equivalently,

$$
b_{i+\frac12,j}^{*}
=
\max\left(
b_{i+\frac12,j}^{-},
b_{i+\frac12,j}^{+}
\right),
$$

$$
b_{i,j+\frac12}^{*}
=
\max\left(
b_{i,j+\frac12}^{-},
b_{i,j+\frac12}^{+}
\right).
$$

This implementation does not reconstruct a piecewise linear bottom and therefore does not require a slope limiter for $b$. It provides the standard piecewise constant bottom treatment used by the first-order well-balanced hydrostatic reconstruction.

Although the flow variables may be reconstructed with MUSCL, the bottom remains piecewise constant in this implementation.

---

## Reconstruct free surface and velocities

From the reconstructed variables, define

$$
\eta^- = (W^-)_{1}, \qquad \eta^+ = (W^+)_{1},
$$

and the reconstructed bottom values $b^-, b^+$ from the bottom reconstruction.

Define the provisional water heights

$$
h^- = \eta^- - b^-,
\qquad
h^+ = \eta^+ - b^+.
$$

Velocities are computed safely:

$$
u^-=
\begin{cases}
\dfrac{(hu)^-}{h^-}, & h^- \ge h_{\mathrm{dry}},\\[2mm]
0, & h^- < h_{\mathrm{dry}},
\end{cases}
\qquad
v^-=
\begin{cases}
\dfrac{(hv)^-}{h^-}, & h^- \ge h_{\mathrm{dry}},\\[2mm]
0, & h^- < h_{\mathrm{dry}},
\end{cases}
$$

and similarly for the right/top states.

---

## Hydrostatic reconstruction at x-interfaces

At interface $\left(i+\frac12,j\right)$, define the interface bottom by

$$
b_{i+\frac12,j}^{*}=
\max\left(b_{i+\frac12,j}^{-},\,b_{i+\frac12,j}^{+}\right).
$$

Then define corrected water heights

$$
h_{i+\frac12,j}^{*, -}=
\max\left(0,\eta_{i+\frac12,j}^{-}-b_{i+\frac12,j}^{*}\right),
$$

$$
h_{i+\frac12,j}^{*, +}=
\max\left(0,\eta_{i+\frac12,j}^{+}-b_{i+\frac12,j}^{*}\right).
$$

The corrected left and right states are then

$$
U_{L}^{*}=
\begin{bmatrix}
h_{i+\frac12,j}^{*,-}\\
h_{i+\frac12,j}^{*,-}\,u^-\\
h_{i+\frac12,j}^{*,-}\,v^-
\end{bmatrix},
\qquad
U_{R}^{*}=
\begin{bmatrix}
h_{i+\frac12,j}^{*,+}\\
h_{i+\frac12,j}^{*,+}\,u^+\\
h_{i+\frac12,j}^{*,+}\,v^+
\end{bmatrix}.
$$

This guarantees nonnegative interface water heights.

---

## Hydrostatic reconstruction at y-interfaces

At interface $\left(i,j+\frac12\right)$, define

$$
b_{i,j+\frac12}^{*}=
\max\left(b_{i,j+\frac12}^{-},\,b_{i,j+\frac12}^{+}\right).
$$

Then corrected heights

$$
h_{i,j+\frac12}^{*,-}=
\max\left(0,\eta_{i,j+\frac12}^{-}-b_{i,j+\frac12}^{*}\right),
$$

$$
h_{i,j+\frac12}^{*,+}=
\max\left(0,\eta_{i,j+\frac12}^{+}-b_{i,j+\frac12}^{*}\right).
$$

Corrected bottom/top states:

$$
U_{B}^{*}=
\begin{bmatrix}
h_{i,j+\frac12}^{*,-}\\
h_{i,j+\frac12}^{*,-}\,u^-\\
h_{i,j+\frac12}^{*,-}\,v^-
\end{bmatrix},
\qquad
U_{T}^{*}=
\begin{bmatrix}
h_{i,j+\frac12}^{*,+}\\
h_{i,j+\frac12}^{*,+}\,u^+\\
h_{i,j+\frac12}^{*,+}\,v^+
\end{bmatrix}.
$$

---

## HLL flux

The HLL flux is applied to the **hydrostatically corrected states**.

---

### x-interface flux

At interface $\left(i+\frac12,j\right)$, use

$$
U_L = U_{L}^{*},
\qquad
U_R = U_{R}^{*}.
$$

Define

$$
c_L=\sqrt{g h_L}, \qquad c_R=\sqrt{g h_R}.
$$

Estimate wave speeds

$$
s_L=\min(u_L-c_L,\;u_R-c_R),
\qquad
s_R=\max(u_L+c_L,\;u_R+c_R).
$$

Then

$$
\hat F_{i+\frac12,j}^{\mathrm{HLL}}=
\begin{cases}
F(U_L), & s_L\ge 0,\\[2mm]
\dfrac{s_R F(U_L)-s_L F(U_R)+s_L s_R (U_R-U_L)}{s_R-s_L},
& s_L<0<s_R,\\[4mm]
F(U_R), & s_R\le 0.
\end{cases}
$$

---

### y-interface flux

At interface $\left(i,j+\frac12\right)$, use

$$
U_B = U_{B}^{*},
\qquad
U_T = U_{T}^{*}.
$$

Define

$$
c_B=\sqrt{g h_B}, \qquad c_T=\sqrt{g h_T}.
$$

Estimate

$$
s_B=\min(v_B-c_B,\;v_T-c_T),
\qquad
s_T=\max(v_B+c_B,\;v_T+c_T).
$$

Then

$$
\hat G_{i,j+\frac12}^{\mathrm{HLL}}=
\begin{cases}
G(U_B), & s_B\ge 0,\\[2mm]
\dfrac{s_T G(U_B)-s_B G(U_T)+s_B s_T (U_T-U_B)}{s_T-s_B},
& s_B<0<s_T,\\[4mm]
G(U_T), & s_T\le 0.
\end{cases}
$$

---

## Hydrostatic pressure correction

Hydrostatic reconstruction modifies the interface states. To retain the correct balance with bottom topography, one usually adds a pressure correction term.

For x-direction, define

$$
\Phi_x(U,h^*)=
\begin{bmatrix}
0\\
\frac12 g \left(h^2-(h^*)^2\right)\\
0
\end{bmatrix}.
$$

Then the corrected interface fluxes are

$$
\hat F_{i+\frac12,j}^{+}=
\hat F_{i+\frac12,j}^{\mathrm{HLL}}+
\Phi_x(U_{i+\frac12,j}^{+},\,h_{i+\frac12,j}^{*,+}),
$$

$$
\hat F_{i+\frac12,j}^{-}=
\hat F_{i+\frac12,j}^{\mathrm{HLL}}+
\Phi_x(U_{i+\frac12,j}^{-},\,h_{i+\frac12,j}^{*,-}),
$$

Similarly in y-direction define

$$
\Phi_y(U,h^*)=
\begin{bmatrix}
0\\
0\\
\frac12 g \left(h^2-(h^*)^2\right)
\end{bmatrix}.
$$

Then

$$
\hat G_{i,j+\frac12}^{+}=
\hat G_{i,j+\frac12}^{\mathrm{HLL}}+
\Phi_y(U_{i,j+\frac12}^{+},\,h_{i,j+\frac12}^{*,+}),
$$

$$
\hat G_{i,j+\frac12}^{-}=
\hat G_{i,j+\frac12}^{\mathrm{HLL}}+
\Phi_y(U_{i,j+\frac12}^{-},\,h_{i,j+\frac12}^{*,-}).
$$

---

## Finite volume operator

The semidiscrete FV operator becomes

$$
\frac{d}{dt}\bar U_{ij}=
-\frac{1}{\Delta x}
\left(
\hat F_{i+\frac12,j}^{-}-
\hat F_{i-\frac12,j}^{+}
\right)
-\frac{1}{\Delta y}
\left(
\hat G_{i,j+\frac12}^{-}-
\hat G_{i,j-\frac12}^{+}
\right).
$$

So define

$$
L(\bar U,\bar b)_{ij}=
-\frac{1}{\Delta x}
\left(
\hat F_{i+\frac12,j}^{-}-
\hat F_{i-\frac12,j}^{+}
\right)
-\frac{1}{\Delta y}
\left(
\hat G_{i,j+\frac12}^{-}-
\hat G_{i,j-\frac12}^{+}
\right).
$$

This is the well-balanced hydrostatic-reconstruction form.

---

## Dry-state cleanup after each stage

After each RK stage, apply a cleanup:

### Height positivity

$$
h_{ij} \leftarrow \max(h_{ij},0).
$$

### Dry-cell reset

If

$$
h_{ij}<h_{\mathrm{dry}},
$$

then set

$$
h_{ij}=0,\qquad (hu)_{ij}=0,\qquad (hv)_{ij}=0.
$$

This is a very common practical rule.

---

## Reflecting walls

At boundaries use ghost cells.

At a vertical wall:

$$
h_{ghost}=h_{inside},
\qquad
(hu)_{ghost}=-(hu)_{inside},
\qquad
(hv)_{ghost}=(hv)_{inside},
\qquad
b_{ghost}=b_{inside}.
$$

At a horizontal wall:

$$
h_{ghost}=h_{inside},
\qquad
(hu)_{ghost}=(hu)_{inside},
\qquad
(hv)_{ghost}=-(hv)_{inside},
\qquad
b_{ghost}=b_{inside}.
$$

If the inside cell is dry, the ghost cell should also be dry.

---

## SSP-RK3

Use the same time stepping as before.

**Stage 1**

$$
U^{(1)} = U^n + \Delta t\,L(U^n,b)
$$

Then apply dry-state cleanup.

**Stage 2**

$$
U^{(2)} = \frac34 U^n + \frac14\left(U^{(1)}+\Delta t\,L(U^{(1)},b)\right)
$$

Then apply dry-state cleanup.

**Stage 3**

$$
U^{n+1} = \frac13 U^n + \frac23\left(U^{(2)}+\Delta t\,L(U^{(2)},b)\right)
$$

Then apply dry-state cleanup again.

---

## CFL

A standard CFL restriction is

$$
\Delta t=
\mathrm{CFL}\min_{i,j}
\left(
\frac{1}{
\dfrac{|u_{ij}|+\sqrt{g h_{ij}}}{\Delta x}+
\dfrac{|v_{ij}|+\sqrt{g h_{ij}}}{\Delta y}
}
\right),
$$

where dry cells are ignored or treated with zero velocity and zero wave speed.

Near wet/dry fronts one often uses a slightly smaller CFL number for robustness.

---

## Practical remarks

- **HLL + hydrostatic reconstruction** is a very robust default for shallow water with topography.
- **Hydrostatic reconstruction** is much better than a naive centered bottom source.
- **Dry threshold** is essential in practice.
- In dry or almost dry cells, always avoid computing $u=hu/h$ and $v=hv/h$ directly.
- After every RK stage, enforce

  $$
  h\ge 0
  $$

  and reset momentum in dry cells.
- For a first code, piecewise constant bottom reconstruction is usually enough.
- First order plus hydrostatic reconstruction is easier to debug than full MUSCL near wet/dry fronts.
- After that works, add MUSCL slopes.

---

## Implementation outline

For each RK stage:

1. fill ghost cells,
2. mark dry cells using $h_{\mathrm{dry}}$,
3. reconstruct MUSCL states,
4. assign piecewise constant bottom values on both sides of each interface,
5. compute free surface $\eta$ on both sides from the reconstructed variables,
6. apply hydrostatic reconstruction to get corrected nonnegative interface states,
7. compute HLL fluxes from corrected states,
8. add hydrostatic pressure correction,
9. build the FV operator,
10. update the RK stage,
11. clip negative heights,
12. reset momentum in dry cells.

---

## In one sentence

This scheme combines:
- **MUSCL** for second-order spatial reconstruction,
- **HLL** for robust fluxes,
- **hydrostatic reconstruction** for well-balanced bottom treatment,
- **dry threshold + cleanup** for stable wet/dry handling.