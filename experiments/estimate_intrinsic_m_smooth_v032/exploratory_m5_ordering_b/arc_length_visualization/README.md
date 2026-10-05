# True B trajectories before and after arc-length reparameterization

The figures show two true feature functions in the original latent coordinate t
and the same functions reparameterized by arc length v. The feature values and
geometric path stay unchanged; the spacing along the parameter axis changes.

V13 and V14 are the first two lexicographically sorted non-anchor B features in
`main_M5_S4_r001`. Selection is independent of method performance. Their true
functions are reconstructed from the saved sine/cosine coefficients at frequencies
2:4, angular frequencies pi*k, attenuation 1/k^2, and saved centers/scales. They
match the saved 300-sample signal exactly. No noisy observations or fitted method
trajectories enter this visualization; the functions are evaluated over t in [0,1].

For f(t)=(f_V13(t),f_V14(t)), analytic derivatives give speed a(t)=||f'(t)||.
Cumulative trapezoidal integration on 20001 regular t nodes gives
v(t)=integral_0^t a(u)du. Linear inversion followed by evaluation of the exact
functions gives g(v)=f(t(v)). The minimum speed is 2.48566, so the coordinate map
is strictly increasing; maximum speed is 14.27311. Total length is L=8.83929.
The arc-length functions satisfy ||g'(v)||=1. With a [0,1] coordinate w=v/L,
||d g(Lw)/dw||=L instead. Feature amplitudes are not normalized again.

Arc length is defined here in the selected two-feature plane, using the original
standardized feature units. Arc length of the full 12-feature B trajectory would
be different: it would include squared derivatives of all 12 features. This
illustration demonstrates the reparameterization of the two-dimensional curve,
not the latent coordinates from an actual principal-curve fit.

- `transformation_overview.png` / `.pdf`: original and reparameterized feature
  functions, then the same geometric curve with equal-t versus equal-v markers.
- `true_functions_t_vs_arc_length.png` / `.pdf`: adds the coordinate map and speed.
- `normalized_arc_length_functions.png`: plots g(Lw) against w in [0,1].
- `true_t_and_arc_length.csv`, `reparameterized_functions.csv`: numerical curves.
- `provenance.rds`: selected generating coefficients, input hash, grid and metric.
- `verification.csv`: signal reconstruction, quadrature, inverse-map and unit-speed checks.

Equal arc-length markers are equally spaced along the path, not necessarily in
straight-line distance between markers. Reparameterization need not make each
individual feature function smoother: for example the V14 trough becomes visibly
sharper on the v scale. Unit speed constrains the joint vector derivative norm,
not each feature derivative separately, and does not enforce constant curvature.

Verification: reconstruction error is zero; total-length difference from a
10001-node grid is 8.8e-8; function round-trip error is below 1.7e-7; analytic
unit-speed error is below 2.3e-16. The figures were visually inspected. Reproduce
from the InferOrder root with the study `run_r.sh` wrapper and this directory's
`render.R`. No website or package sources are modified.
