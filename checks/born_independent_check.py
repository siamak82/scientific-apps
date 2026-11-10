"""Independent check supplied by the user from Claude's review.
Reference: own welded-interface plane-wave boundary-condition solve.
Usage: python3 checks/born_independent_check.py dist/physics.js
Needs node and numpy. Original user-supplied calculations are unchanged.
"""
import json, subprocess, sys
import numpy as np
PHYS = sys.argv[1]
DUMP = r"""
const {kernel,parameters}=require(require("path").resolve(process.argv[1]));
const s=JSON.parse(process.argv[2]),m=process.argv[3],a=+process.argv[4];
console.log(JSON.stringify(parameters(s).map(p=>kernel(s,a,m,p))));
"""
def app_kernel(state, mode, theta_deg):
    out = subprocess.run(["node", "-e", DUMP, PHYS, json.dumps(state), mode, str(theta_deg)],capture_output=True, text=True, check=True).stdout
    return [complex(*c) for c in json.loads(out)]
def wave(med, kind, p, sgn):
    r, a, b = med
    v = a if kind == "P" else b
    q = sgn * np.sqrt(1 / v**2 - p**2 + 0j)
    s = np.array([p, 0, q]); n = s * v
    d = n if kind == "P" else np.array([n[2], 0, -n[0]])
    lam, mu = r * (a**2 - 2 * b**2), r * b**2
    return np.array([d[0], d[2], mu * (s[0] * d[2] + s[2] * d[0]), lam * (s @ d) + 2 * mu * s[2] * d[2]])
def psv(m1, m2, p, inc):
    A = np.column_stack([-wave(m1, "P", p, -1), -wave(m1, "S", p, -1),wave(m2, "P", p, +1), wave(m2, "S", p, +1)])
    return np.linalg.solve(A, wave(m1, inc, p, +1))
def sh(m1, m2, p):
    (r1, _, b1), (r2, _, b2) = m1, m2
    q1, q2 = np.sqrt(1 / b1**2 - p**2 + 0j), np.sqrt(1 / b2**2 - p**2 + 0j)
    return (r1 * b1**2 * q1 - r2 * b2**2 * q2) / (r1 * b1**2 * q1 + r2 * b2**2 * q2)
vp, vs, rho, eps = 3.2, 1.4, 2.2, 1e-6
med = (rho, vp, vs)
state = dict(elastic=True, visco=False, basis="velocity", qcoord="quality",vp=vp, vs=vs, rho=rho, qp=80, qs=55, f=10, fref=10)
def perturbed(k):
    d = [0.0, 0.0, 0.0]; d[k] = eps
    return (rho * np.exp(d[2]), vp * np.exp(d[0]), vs * np.exp(d[1]))
worst = 0.0
for mode in ["PP", "PS", "SP", "SS", "HH"]:
    inc = "P" if mode[0] == "P" else "S"
    vi = vp if mode[0] == "P" else vs
    vo = vp if mode[1] == "P" else vs
    for ideg in [5, 15, 25, 35]:
        i = np.radians(ideg); p = np.sin(i) / vi
        if p * vo >= 1: continue
        j = np.arcsin(p * vo)
        B = app_kernel(state, mode, 180 - ideg - np.degrees(j))
        born = [-b.real / (2 * vo * np.cos(j) * (np.cos(i) / vi + np.cos(j) / vo)) for b in B]
        if mode == "HH": ex = [(sh(med, perturbed(k), p) - sh(med, med, p)).real / eps for k in range(3)]
        else:
            idx = 0 if mode[1] == "P" else 1
            ex = [(psv(med, perturbed(k), p, inc)[idx] - psv(med, med, p, inc)[idx]).real / eps for k in range(3)]
        worst = max(worst, max(abs(x - y) for x, y in zip(born, ex)))
for ideg in [3, 8, 15, 20, 25]:
    i = np.radians(ideg); p = np.sin(i) / vs; j = np.arcsin(p * vp)
    B = app_kernel(state, "SP", np.degrees(j - i))
    born = [-b.real / (2 * vp * np.cos(j) * (np.cos(i) / vs - np.cos(j) / vp)) for b in B]
    ex = [(psv(med, perturbed(k), p, "S")[2] - psv(med, med, p, "S")[2]).real / eps for k in range(3)]
    worst = max(worst, max(abs(x - y) for x, y in zip(born, ex)))
print(f"Lossless elastic, 5 modes, reflection + forward transmission: max error {worst:.1e}" f"  -> {'PASS' if worst < 1e-5 else 'FAIL'}")
Q, f, fr = 20.0, 10.0, 10.0
M = rho * vp**2
L = 2 / np.pi * np.log(f / fr) - 1j
Ms = M * (1 + L / Q)
def R(r2, M2):
    z1, z2 = np.sqrt(rho * Ms), np.sqrt(r2 * M2)
    return (z2 - z1) / (z1 + z2)
vstate = dict(elastic=False, visco=True, basis="velocity", qcoord="quality",vp=vp, vs=0, rho=rho, qp=Q, qs=50, f=f, fref=fr)
B = app_kernel(vstate, "PP", 180)
cols = [(0, 2 * Ms), (rho, Ms), (0, -M / Q * L)]
bad = 0.0
for name, (dr, dM), b in zip(["ln Vp", "ln rho", "ln Qp"], cols, B):
    exact = (R(rho + eps * dr, Ms + eps * dM) - R(rho, Ms)) / eps
    bad = max(bad, abs(-b / 4 - exact))
    print(f"  visco {name:6s}: exact {exact:.5f}   app {-b / 4:.5f}")
print(f"Viscoacoustic normal incidence (Q=20): max error {bad:.1e}  -> {'PASS' if bad < 1e-5 else 'FAIL'}")

# Fail CI on a printed failure without changing the reference calculations.
sys.exit(0 if worst < 1e-5 and bad < 1e-5 else 1)
