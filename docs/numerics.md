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

or derive on my own...