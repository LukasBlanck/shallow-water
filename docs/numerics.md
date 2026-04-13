## PDE
$$
U = \begin{bmatrix}
h \\
hu \\
hv
\end{bmatrix}
$$

### Homogeneous, flat bottom
$
{\partial d_t} U  + {\partial d_x} F(U)   + {\partial d_y} G(U) = 0    $


$$
F(U) =
\begin{bmatrix}
hu \\
hu^2 + \frac{1}{2}gh^2 \\
huv
\end{bmatrix}
$$

$$
G(U) =
\begin{bmatrix}
hv \\
huv  \\
hv^2 + \frac{1}{2}gh^2
\end{bmatrix}
$$



h(x,y,t) .... water height
u(x,y,t) .... x-velocity
v(x,y,t) .... y-velocity

## Reflecting Walls

normal velocity components at walls are set to zero, tangential is mirrored?

## Grid
Domain $$\Omega = [0, L_x][0, L_y]$$

$$
\Delta x = \frac{L_x}{N_x} \qquad\Delta y = \frac{L_y}{N_y}
$$
$N_{x,y}$ ... number of cells

##### Cell centers

$$
x_i = (i + \frac{1}{2})\Delta x, \qquad i = 0, ..., N_x - 1
$$
$$
y_i = (i + \frac{1}{2})\Delta y, \qquad i = 0, ..., N_y - 1
$$

##### Control volumes
$$
C_{ij} = [x_{i-\frac{1}{2}}, x_{i+\frac{1}{2}}][y_{i-\frac{1}{2}}, y_{i+\frac{1}{2}}]
$$
area: $$|C_{ij}| = \Delta x \Delta y$$

#### Cell averages

FV does not store points values but cell averages:
$$
\bar{U}_{ij}(t)=\frac{1}{\Delta x\,\Delta y}\int_{C_{ij}} U(x,y,t)\,dx\,dy
$$

$$
\bar{U}_{ij}(t)=
\begin{bmatrix}
\bar{h}_{ij}(t) \\
\overline{hu}_{ij}(t) \\
\overline{hv}_{ij}(t)
\end{bmatrix}
$$

## Finite Volume
$$
\frac{d}{dt} \bar{U}_{ij}(t) = - \frac{1}{\Delta x} (\hat{F}_{i+\frac{1}{2},j} - \hat{F}_{i-\frac{1}{2},j}) -  \frac{1}{\Delta y} (\hat{G}_{i,j+\frac{1}{2}} - \hat{G}_{i,j-\frac{1}{2}}) = L(\bar{U})_(ij)
$$
$\hat{F}$ ... flux in x-direction (average over y)

## Riemann Integrators

#### Rusanov: 
###### simple robust and easy implementation
###### x-interface flux

At interface $\left(i+\tfrac{1}{2},\,j\right)$, left and right states $U_L$ and $U_R$. Then

$$
\hat{F}_{i+\frac{1}{2},j} = \frac{1}{2}\Bigl(F(U_L)+ F(U_R)\Bigr) - \frac{1}{2}a_{i+\frac{1}{2},j}(U_R-U_L)
$$

with dissipation speed

$$
a_{i+\frac{1}{2},j} = \max\left(|u_L|+\sqrt{g h_L},\,|u_R|+\sqrt{g h_R}\right).
$$

###### y-interface flux

Similarly,

$$
\hat{G}_{i,j+\frac{1}{2}}
=\frac{1}{2}\Bigl(G(U_B)+G(U_T)\Bigr)-\frac{1}{2}b_{i,j+\frac{1}{2}}(U_T-U_B)
$$

with

$$
b_{i,j+\frac{1}{2}} =\max\left(|v_B|+\sqrt{g h_B},\,|v_T|+\sqrt{g h_T}\right).
$$

Here:
- $U_B$ = state below the horizontal interface
- $U_T$ = state above it

for now I implement picewise constant reconstruction (first order):

$$U_L = \bar{U}_{ij} \qquad U_R = \bar{U}_{i+1,j}$$
$$U_B = \bar{U}_{ij} \qquad U_T = \bar{U}_{i,j+1}$$

So we have a spatial operator $L$:
$$
L(\bar{U})_{ij} = - \frac{1}{\Delta x} (\hat{F}_{i+\frac{1}{2},j} - \hat{F}_{i-\frac{1}{2},j}) -  \frac{1}{\Delta y} (\hat{G}_{i,j+\frac{1}{2}} - \hat{G}_{i,j-\frac{1}{2}})
$$

## Reflecting Walls

At boundary ghost cells are introduced, that mirror the momentum.  $U_{ghost}$


### SSP-RK3

**Stage 1**
$$
U^{(1)} = U^n + \Delta t\,L(U^n)
$$

**Stage 2**
$$
U^{(2)} = \frac{3}{4}U^n + \frac{1}{4}\left(U^{(1)} + \Delta t\,L(U^{(1)})\right)
$$

**Stage 3**
$$
U^{n+1} = \frac{1}{3}U^n + \frac{2}{3}\left(U^{(2)} + \Delta t\,L(U^{(2)})\right)
$$

## CFL
$$
\Delta t
=\mathrm{CFL}\min_{i,j}\left(\frac{1}{\frac{|u_{ij}|+\sqrt{g h_{ij}}}{\Delta x}+\frac{|v_{ij}|+\sqrt{g h_{ij}}}{\Delta y}}\right)
$$

# Numerical scheme for `Reconstruction::MUSCL` and `Riemann::HLL` / `Riemann::ROE`


---

## Reconstruction

### Piecewise constant (first order)
I had

$$
U_L = \bar U_{ij}, \qquad U_R = \bar U_{i+1,j},
$$

$$
U_B = \bar U_{ij}, \qquad U_T = \bar U_{i,j+1}.
$$

---

### MUSCL (piecewise linear, second order in space)

Instead of a constant state in each cell, I reconstruct a linear profile from cell averages.

For each cell $C_{ij}$, define limited slopes in $x$- and $y$-direction.

#### x-slopes
A simple choice is the **minmod** limiter applied componentwise:

$$
\sigma^x_{ij}=
\operatorname{minmod}\!\left(
\bar U_{ij}-\bar U_{i-1,j},
\bar U_{i+1,j}-\bar U_{ij}
\right).
$$

#### y-slopes
Similarly,

$$
\sigma^y_{ij}=
\operatorname{minmod}\!\left(
\bar U_{ij}-\bar U_{i,j-1},
\bar U_{i,j+1}-\bar U_{ij}
\right).
$$

The minmod function is applied **componentwise**:
$$
\operatorname{minmod}(a,b)=
\begin{cases}
\operatorname{sign}(a)\min(|a|,|b|), & ab>0,\\[1mm]
0, & ab\le 0.
\end{cases}
$$

---


#### At an x-interface $\left(i+\tfrac12,j\right)$

Left state from cell $C_{ij}$:
$$
U_{i+\frac12,j}^{-}=
\bar U_{ij}
+\frac12\,\sigma^x_{ij}.
$$

Right state from cell $C_{i+1,j}$:
$$
U_{i+\frac12,j}^{+}=
\bar U_{i+1,j}
-\frac12\,\sigma^x_{i+1,j}.
$$

So for the Riemann solver at this interface:
$$
U_L = U_{i+\frac12,j}^{-},\qquad U_R = U_{i+\frac12,j}^{+}.
$$

#### At a y-interface $\left(i,j+\tfrac12\right)$

Bottom state from cell $C_{ij}$:
$$
U_{i,j+\frac12}^{-}=
\bar U_{ij}
+\frac12\,\sigma^y_{ij}.
$$

Top state from cell $C_{i,j+1}$:
$$
U_{i,j+\frac12}^{+}=
\bar U_{i,j+1}
-\frac12\,\sigma^y_{i,j+1}.
$$

So:
$$
U_B = U_{i,j+\frac12}^{-}, \qquad U_T = U_{i,j+\frac12}^{+}.
$$

---

## Riemann Integrators

## HLL

#### x-interface flux

At interface $\left(i+\tfrac12,j\right)$, let
$$
U_L = U_{i+\frac12,j}^{-},
\qquad
U_R = U_{i+\frac12,j}^{+}.
$$

Define
$$
c_L=\sqrt{g h_L}, \qquad c_R=\sqrt{g h_R}.
$$

Estimate the left and right wave speeds by
$$
s_L = \min(u_L-c_L,\;u_R-c_R),
\qquad
s_R = \max(u_L+c_L,\;u_R+c_R).
$$

Then the HLL flux is

$$
\hat F_{i+\frac12,j}=
\begin{cases}
F(U_L), & s_L \ge 0,\\[2mm]
\dfrac{s_R F(U_L)-s_L F(U_R)+s_L s_R (U_R-U_L)}{s_R-s_L},
& s_L < 0 < s_R,\\[4mm]
F(U_R), & s_R \le 0.
\end{cases}
$$

---

#### y-interface flux

At interface $\left(i,j+\tfrac12\right)$, let
$$
U_B = U_{i,j+\frac12}^{-},
\qquad
U_T = U_{i,j+\frac12}^{+}.
$$

Define
$$
c_B=\sqrt{g h_B}, \qquad c_T=\sqrt{g h_T}.
$$

Estimate
$$
s_B = \min(v_B-c_B,\;v_T-c_T),
\qquad
s_T = \max(v_B+c_B,\;v_T+c_T).
$$

Then

$$
\hat G_{i,j+\frac12}=
\begin{cases}
G(U_B), & s_B \ge 0,\\[2mm]
\dfrac{s_T G(U_B)-s_B G(U_T)+s_B s_T (U_T-U_B)}{s_T-s_B},
& s_B < 0 < s_T,\\[4mm]
G(U_T), & s_T \le 0.
\end{cases}
$$

---

## ROE

The Roe flux has the form

$$
\hat F_{i+\frac12,j}=
\frac12\Bigl(F(U_L)+F(U_R)\Bigr)-
\frac12 \sum_{k=1}^3 |\lambda_k|\,\alpha_k\,r_k
$$

in x-direction, and analogously in y-direction.

---

### Roe flux in x-direction

At interface $\left(i+\tfrac12,j\right)$, with states $U_L,U_R$,  primitive variables are

$$
u_L = \frac{(hu)_L}{h_L}, \qquad v_L = \frac{(hv)_L}{h_L},
$$

$$
u_R = \frac{(hu)_R}{h_R}, \qquad v_R = \frac{(hv)_R}{h_R}.
$$

#### Roe averages

$$
\tilde u=
\frac{\sqrt{h_L}\,u_L+\sqrt{h_R}\,u_R}{\sqrt{h_L}+\sqrt{h_R}},
\qquad
\tilde v=
\frac{\sqrt{h_L}\,v_L+\sqrt{h_R}\,v_R}{\sqrt{h_L}+\sqrt{h_R}}.
$$

$$
\tilde h = \frac{h_L+h_R}{2},
\qquad
\tilde c = \sqrt{g\tilde h}.
$$

#### Eigenvalues

$$
\lambda_1 = \tilde u-\tilde c,
\qquad
\lambda_2 = \tilde u,
\qquad
\lambda_3 = \tilde u+\tilde c.
$$

#### Right eigenvectors

$$
r_1 =
\begin{bmatrix}
1\\
\tilde u-\tilde c\\
\tilde v
\end{bmatrix},
\qquad
r_2 =
\begin{bmatrix}
0\\
0\\
1
\end{bmatrix},
\qquad
r_3 =
\begin{bmatrix}
1\\
\tilde u+\tilde c\\
\tilde v
\end{bmatrix}.
$$

#### Wave strengths

Let
$$
\Delta U = U_R-U_L=
\begin{bmatrix}
\Delta h\\
\Delta(hu)\\
\Delta(hv)
\end{bmatrix}.
$$

Then
$$
\alpha_1=
\frac12\left(
\Delta h-
\frac{\Delta(hu)-\tilde u\,\Delta h}{\tilde c}
\right),
$$

$$
\alpha_3=\frac12\left(
\Delta h+
\frac{\Delta(hu)-\tilde u\,\Delta h}{\tilde c}\right),
$$

$$
\alpha_2=\Delta(hv)-\tilde v\,\Delta h.
$$

Hence

$$
\hat F_{i+\frac12,j}=
\frac12\Bigl(F(U_L)+F(U_R)\Bigr)-\frac12\left(
|\lambda_1|\alpha_1 r_1+
|\lambda_2|\alpha_2 r_2+
|\lambda_3|\alpha_3 r_3
\right).
$$

---

### Roe flux in y-direction

At interface $\left(i,j+\tfrac12\right)$, with states $U_B,U_T$, define Roe averages

$$
\tilde u
= \frac{\sqrt{h_B}\,u_B+\sqrt{h_T}\,u_T}{\sqrt{h_B} +\sqrt{h_T}},\qquad \tilde v=\frac{\sqrt{h_B}\,v_B+\sqrt{h_T}\,v_T}{\sqrt{h_B}+\sqrt{h_T}}.
$$

$$
\tilde h = \frac{h_B+h_T}{2},
\qquad
\tilde c = \sqrt{g\tilde h}.
$$

#### Eigenvalues

$$
\mu_1 = \tilde v-\tilde c,
\qquad
\mu_2 = \tilde v,
\qquad
\mu_3 = \tilde v+\tilde c.
$$

#### Right eigenvectors

$$
s_1 =
\begin{bmatrix}
1\\
\tilde u\\
\tilde v-\tilde c
\end{bmatrix},
\qquad
s_2 =
\begin{bmatrix}
0\\
1\\
0\end{bmatrix},\qquad s_3 =\begin{bmatrix}
1\\
\tilde u\\
\tilde v+\tilde c
\end{bmatrix}.
$$

#### Wave strengths

Let
$$
\Delta U = U_T-U_B=\begin{bmatrix}
\Delta h\\
\Delta(hu)\\
\Delta(hv)
\end{bmatrix}.
$$

Then
$$
\beta_1=\frac12\left(\Delta h-\frac{\Delta(hv)-\tilde v\,\Delta h}{\tilde c}\right),
$$

$$
\beta_3=\frac12\left(\Delta h+\frac{\Delta(hv)-\tilde v\,\Delta h}{\tilde c}\right),
$$

$$
\beta_2=\Delta(hu)-\tilde u\,\Delta h.
$$

Therefore

$$
\hat G_{i,j+\frac12}=\frac12\Bigl(G(U_B)+G(U_T)\Bigr)-\frac12\left(|\mu_1|\beta_1 s_1+|\mu_2|\beta_2 s_2+|\mu_3|\beta_3 s_3
\right).
$$

---

## Finite volume operator with MUSCL + HLL / Roe

The semi-discrete scheme keeps the same form:

$$
\frac{d}{dt} \bar{U}_{ij}(t)=- \frac{1}{\Delta x}
\left(
\hat{F}_{i+\frac{1}{2},j}-
\hat{F}_{i-\frac{1}{2},j}
\right)-\frac{1}{\Delta y}
\left(\hat{G}_{i,j+\frac{1}{2}}-
\hat{G}_{i,j-\frac{1}{2}}
\right)=L(\bar U)_{ij}.
$$

The only difference is:

- **first order**: interface states are piecewise constant,
- **MUSCL**: interface states are reconstructed linearly,
- **Rusanov / HLL / Roe**: different numerical flux functions.

Then SSP-RK3 is applied exactly as before.

---


## Practical remarks

- **MUSCL + HLL** is usually a very good robust default.
- **Roe** is sharper than HLL, but usually needs an **entropy fix** near transonic rarefactions.
- For shallow water with very small $h$, **HLL is safer** than Roe.
- With MUSCL, it is common to limit slopes in primitive variables $(h,u,v)$ instead of conserved variables $U$, but componentwise limiting of $U$ is the simplest starting point. $(h,u,v)$ avoids a mixing of $(h,u,v)$ in $(h,hu,hv)$, so that if h changes strongly and u doesnt, hu still changes strongly.
