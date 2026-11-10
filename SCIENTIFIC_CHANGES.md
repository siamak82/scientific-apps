# Corrections following the independent review — 2026-11-09

The user supplied Claude's review and its welded-interface check script. This snapshot includes that script in `checks/born_independent_check.py`, preserving its computations. It is run against this app; it is not described as a test independently authored by this app.

## Verified before/after

- Baseline independent elastic reflection/transmission test: PASS, maximum error 1.1e-6.
- Baseline viscoacoustic normal incidence: FAIL, maximum error 2.5e-2.
- Corrected elastic test: PASS, maximum error 1.1e-6.
- Corrected viscoacoustic normal incidence: PASS, maximum error 5.0e-7.
- Exact stdout: `checks/review_baseline.txt` and `checks/review_corrected.txt`.

## Implemented corrections

1. Complex homogeneous-wave scalar slowness:
   k_i k_o / omega² = rho / (sqrt(M_i*) sqrt(M_o*)), replacing 1/(v_i v_o).
   Square-root branches have positive real part, with passive negative-imaginary moduli under exp(-i omega t). Propagation directions and homogeneous isotropic wave polarizations remain real. Inhomogeneous interface-generated wave propagation remains outside this point-scatterer app.

2. Passivity rejection, retaining existing Qp/Qs coordinates:
   M/Qp >= 4 mu/(3 Qs), or Qp <= (3/4)(Vp/Vs)² Qs.
   The app rejects invalid backgrounds and clears prior results. It does NOT switch to bulk/shear attenuation coordinates. Modulus derivatives still hold Qp/Qs fixed, and can therefore change the derived bulk Q. At a passivity boundary, admissible perturbations would require constrained inversion; this app does not perform inversion.

3. Single-mode matrices by default, following the plotted mode. Stacking is opt-in.
   Choices: equalized kernels, or unit incident displacement in each selected experiment.
   Displacement mode weights are |Mref / M_out*|²; in the lossless limit this gives (Vp/v_out)^4. Mref=rho Vp² is the real reference modulus. Phase of the common outgoing prefactor cancels from each mode's real-parameter Gram products. Source radiation, receiver-component projections, path losses and survey noise weights remain omitted. This is an ideal weighting convention, not a survey-calibrated displacement Hessian.

4. Selectable angular measure: uniform dtheta on the plotted angular section, or sin(theta) dtheta for the corresponding ideal 3-D solid-angle integral. Samples remain at one-degree increments, including endpoints; no broadband integration is added.

5. CSV exports now identify the rows included in the Gram matrix and each row's weight. Plot-component selection remains separate from the Gram metric.

6. Tiny real/imaginary components below 1e-12 of a column's full complex magnitude are suppressed only in the plot to prevent roundoff normalization artifacts. Kernel values, Gram calculations and CSV values are not clipped by this display rule.

## Reproduced mode-weighting example

Lossless Vp=3, Vs=1.5, rho=2, PP+PS, full 0–180 degrees, uniform angle measure:

| Weights | C(Vp,Vs) | C(Vp,rho) | C(Vs,rho) |
| --- | ---: | ---: | ---: |
| Equalized kernels | -0.533 | 0.573 | -0.051 |
| Unit incident displacement | -0.172 | 0.226 | 0.382 |

## Q interpretation that needs qualification

The app retains the linearized modulus M*=M[1+q L], with L=(2/pi)ln(f/fref)-i. With corrected kinematics the logarithmic Q column is q L/(1+q L) in the acoustic PP velocity-density-Q basis. At f=fref it has an O(q²) real part and an O(q) imaginary part. It is mainly, but not exactly, imaginary. This matches the supplied exact-impedance test using the same constitutive model.

An exact Kjartansson model is a different constitutive law and needs its own consistent derivatives. It was not silently substituted into this app. The acoustic Vp and density columns in the velocity–density basis are real, including in an attenuating background for this homogeneous-wave angular factor after the complex stiffness/slowness cancellation.

Magnitude-only overlap is an upper bound on the absolute phase-aware real-parameter correlation under identical nonnegative weights; it can equal it. It does not always strictly increase it.

## Checks

Run:

    node checks/physics.cjs
    node checks/parameterizations.cjs
    node checks/review_regressions.cjs
    python3 checks/born_independent_check.py dist/physics.js

The supplied Python reference calculations are unchanged. The script now exits with status 1 if either test exceeds its tolerance, and 0 only when both pass. The existing 3,180 chain-rule finite-difference checks now use independently constructed complex reference slownesses, frozen when differentiating the perturbation operator.

The independent script verifies lossless elastic reflection and selected forward transmission, plus viscoacoustic normal incidence. It does not validate general oblique viscoelastic scattering at inhomogeneous interfaces. That remains an explicit limit of this verification.


## Follow-up review changes

The user interface now explicitly identifies impedance–Q at fixed density as equivalent to velocity–Q. A dedicated basis note explains that fixing Qp and Qs in the ln μ derivative changes bulk QK. A full-angle note identifies the symmetric lossless velocity–density versus impedance–density correlation equivalence. Magnitude-only overlap is labeled as an upper bound on absolute phase-aware correlation for identical weights. A visible low-loss approximation warning appears when max over P and S of |q Re L| exceeds 0.05; this is a diagnostic threshold, not a physical instability boundary, and results remain available for inspection. The existing hard passivity rejection remains separate.
