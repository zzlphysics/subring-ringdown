"""
Photon-Ring Subring Ringdown — Core Module (v2.0)
==================================================

CONVENTIONS (definitive):

  Two observables (q = 1: half-orbit adjacent, q = 2: same-family):

    R^{pix}_{n,q} = ΔP_{n+q}(φ, t_o) / ΔP_n(φ, t_o)
      → exp[-q γ] exp[i q (ω τ - m δ)]        [same screen angle]

    R^{align}_{n,q} = ΔP_{n+q}(φ+qδ, t_o) / ΔP_n(φ, t_o)
      → exp[-q γ] exp[i q ω τ]                [pattern-aligned screen position]

  Complex QNM frequency: ω = ω_R - i ω_I

  Therefore:
    |R_raw|  = exp[-q γ + q ω_I τ]       ← full raw magnitude (includes ringdown envelope)
    |R_geom| = exp[-q γ]                 ← geometric demagnification only
    
    arg(R^pix)   = q (ω_R τ - m δ)  (mod 2π)
    arg(R^align) = q ω_R τ           (mod 2π)

  The phase is measured modulo 2π. The integer branch k is fixed by a Kerr/QNM prior.

  BASIS CONVENTION: R^{align} compares different screen positions (φ→φ+qδ) but
  Q+iU is ALWAYS projected onto a fixed Cartesian camera tetrad at the observer.
  The Stokes basis does NOT rotate with the image position.
  Otherwise a spin-2 basis phase would contaminate arg(R^align).

REFERENCES:
  [GL2020]  Gralla & Lupsasca, PRD 101, 044031 (1910.12873)
  [H2026]   Huang, Hou, Zhong, Guo, Chen (2605.11499)
  [ZC2025]  Zhong, Cardoso, Chen, PRL 134, 211402 (2408.10303)
"""
import numpy as np
from scipy.special import ellipk, ellipe
from numpy import sqrt, pi, cos

M_UNIT = 1.0

# ═══════════════════════════════════════════════════════════════
# 1. Photon Ring Critical Parameters (Gralla+Lupsasca 2020)
# ═══════════════════════════════════════════════════════════════

def photon_shell_bounds(a):
    a_abs = abs(a)/M_UNIT
    r_pro = 2*M_UNIT*(1 + cos(2/3*np.arccos(-a_abs)))
    r_ret = 2*M_UNIT*(1 + cos(2/3*np.arccos(a_abs)))
    return r_pro, r_ret

def r0_polar(a):
    M = M_UNIT
    if abs(a) < 1e-10: return 3*M
    p = a**2 - 3*M**2
    if p >= 0: return 2*M
    q = 2*a**2*M - 2*M**3
    arg = max(-1.0, min(1.0, -q/2*(-3/p)**1.5))
    phi = np.arccos(arg)
    roots = [M + 2*sqrt(-p/3)*cos(phi/3 + k*2*pi/3) for k in [0,1,2] if M + 2*sqrt(-p/3)*cos(phi/3 + k*2*pi/3) > 0]
    rmin, rmax = photon_shell_bounds(a)
    valid = [r for r in roots if rmin <= r <= rmax]
    return valid[0] if valid else 2*M

def b_tilde_polar(r0, a):
    M = M_UNIT
    if abs(a) < 1e-8: return 3*sqrt(3)*M
    Delta = r0**2 - 2*M*r0 + a**2
    return sqrt(max(r0**3/a**2*(4*M*Delta/(r0-M)**2 - r0) + a**2, 1e-10))

def critical_exponents_polar(a):
    """Returns (gamma, delta[rad], tau[M]) per half-orbit for face-on observer."""
    M = M_UNIT
    r0 = r0_polar(a)
    bt = b_tilde_polar(r0, a)
    Delta = r0**2 - 2*M*r0 + a**2
    k_sq = a**2/(a**2 - bt**2) if abs(bt-a) > 1e-10 else -1e10
    Kk, Ek = ellipk(k_sq), ellipe(k_sq)
    denom = sqrt(max(bt**2 - a**2, 1e-10))
    gamma = (4*r0/denom) * sqrt(max(1 - M*Delta/(r0*(r0-M)**2), 0)) * Kk
    delta = pi + (2*a/denom)*(r0+M)/(r0-M)*Kk
    tau   = (2/denom)*(r0**2*(r0+3*M)/(r0-M)*Kk - 2*a**2*Ek)
    return {'a':a, 'r0':r0, 'b_tilde':bt, 'gamma':gamma, 'delta':delta,
            'delta_deg':np.degrees(delta)%360, 'tau':tau, 'e_minus_gamma':np.exp(-gamma)}

# ═══════════════════════════════════════════════════════════════
# 2. QNM Frequencies
# ═══════════════════════════════════════════════════════════════

_QNM_l2m2 = {
    0.00:(0.373672,0.088962), 0.10:(0.380041,0.088578), 0.20:(0.388944,0.088352),
    0.30:(0.400927,0.088240), 0.40:(0.416864,0.088236), 0.50:(0.437973,0.088303),
    0.60:(0.466188,0.088361), 0.70:(0.504577,0.088318), 0.75:(0.529340,0.088125),
    0.80:(0.558058,0.088030), 0.85:(0.591859,0.087638), 0.90:(0.635237,0.087132),
    0.94:(0.686507,0.085863), 0.95:(0.704543,0.084744), 0.97:(0.749830,0.082989),
    0.98:(0.771320,0.081409), 0.99:(0.790639,0.078695),
}
def qnm_l2m2_fund(a):
    a_abs = abs(a); keys = sorted(_QNM_l2m2)
    wR = np.interp(a_abs, keys, [_QNM_l2m2[k][0] for k in keys])
    wI = np.interp(a_abs, keys, [_QNM_l2m2[k][1] for k in keys])
    return wR, wI

# ═══════════════════════════════════════════════════════════════
# 3. Subring Ratio — CORRECTED
# ═══════════════════════════════════════════════════════════════

def subring_ratio(a, m=2, q=1, aligned=False):
    """Compute the predicted subring ratio in the near-critical limit.

    Returns dict with fields:
      R_raw        : complex, exp[-q γ + i q ω τ] or exp[-q γ + i q(ω τ - m δ)]
      R_geom       : complex, same phase but |R_geom| = exp[-q γ] only (no ω_I)
      |R_raw|      : exp[-q γ + q ω_I τ]    ← full raw magnitude
      |R_geom|     : exp[-q γ]              ← geometric demagnification only
      phase_mod_deg: arg(R) modulo 360°
      phase_unwrapped_deg: continuous (unwrapped) total phase
      omega_R_tau_total_deg: q ω_R τ total (unwrapped), useful for branch k
      branch_k_min : minimum k such that ω_R > 0
    """
    crit = critical_exponents_polar(a)
    g, d, t = crit['gamma'], crit['delta'], crit['tau']
    wR, wI = qnm_l2m2_fund(a)
    omega = wR - 1j*wI
    
    gq, tq, dq = q*g, q*t, q*d
    wR_tq = wR * tq          # total phase (unwrapped), can be > 2π
    wI_tq = wI * tq
    mdq = m * dq
    
    if aligned:
        phase = wR_tq
    else:
        phase = wR_tq - mdq
    
    R_raw = np.exp(-gq + wI_tq + 1j*phase)
    R_geom = np.exp(-gq + 1j*phase)
    
    # Branch info
    # For R^align:  phase_mod = q ω_R τ (mod 2π)
    #   → ω_R = (phase_mod + 2πk) / (qτ)
    # For R^pix:    phase_mod = q(ω_R τ - mδ) (mod 2π)
    #   → ω_R = (phase_mod + 2πk + q m δ) / (qτ)
    phase_mod = phase % (2*pi)
    k_min = int(np.floor(wR_tq / (2*pi)))
    
    result = {
        'label': 'align' if aligned else 'pix',
        'q': q, 'aligned': aligned,
        'gamma_q': gq, 'tau_q': tq, 'delta_q': dq,
        'omega_R': wR, 'omega_I': wI,
        'omega_R_tau_q_rad': wR_tq,
        'omega_R_tau_q_deg': np.degrees(wR_tq),
        'omega_I_tau_q': wI_tq,
        'm_delta_q_rad': mdq,
        'm_delta_q_deg': np.degrees(mdq),
        # Raw ratio (full complex QNM)
        'R_raw': R_raw,
        'amp_raw': abs(R_raw),
        # Geometric only (envelope-removed)
        'R_geom': R_geom,
        'amp_geom': abs(R_geom),
        # Phase
        'phase_mod_deg': np.degrees(phase_mod),
        'phase_unwrapped_deg': np.degrees(phase),
        'phase_total_deg': np.degrees(phase),
        'branch_k_min': k_min,
        # True ω_R recovery candidates
        'omega_R_candidates': _omega_candidates(phase_mod, tq, dq, m, q, aligned),
    }
    return result


def _omega_candidates(phase_mod, tau_q, delta_q, m, q, aligned,
                      k_range=range(-2, 6)):
    """Recover possible ω_R values from a modulo-2π phase measurement.
    
    For R^align: ω_R = (φ_mod + 2πk) / (qτ)
    For R^pix:   ω_R = (φ_mod + 2πk + q m δ) / (qτ)
    
    Returns list of (k, ω_R) tuples with ω_R > 0.
    """
    candidates = []
    for k in k_range:
        if aligned:
            wR = (phase_mod + 2*np.pi*k) / tau_q
        else:
            wR = (phase_mod + 2*np.pi*k + q*m*delta_q) / tau_q
        if wR > 0:
            candidates.append((k, wR))
    return candidates


# ═══════════════════════════════════════════════════════════════
# 4. Convenience
# ═══════════════════════════════════════════════════════════════

def print_ratios(a=0.7, m=2):
    """Print all ratio quantities for given spin."""
    print(f"\n{'='*60}")
    print(f"SUBHING RATIOS: a={a}, m={m}")
    print(f"{'='*60}")
    crit = critical_exponents_polar(a)
    wR, wI = qnm_l2m2_fund(a)
    print(f"  γ={crit['gamma']:.4f}  τ={crit['tau']:.4f}M  δ={crit['delta_deg']:.1f}°")
    print(f"  ω_R={wR:.4f}/M  ω_I={wI:.4f}/M")
    print(f"  ω_R τ = {np.degrees(wR*crit['tau']):.1f}° total (q=1)")
    print(f"  ω_I τ = {wI*crit['tau']:.4f}")
    print()
    for aligned in [False, True]:
        for q in [1, 2]:
            r = subring_ratio(a, m=m, q=q, aligned=aligned)
            lbl = r['label']
            cand = r['omega_R_candidates']
            print(f"  R^{lbl} (q={q}):")
            print(f"    |R_raw|  = {r['amp_raw']:.4f}   (exp[-q*gamma + q*omega_I*tau])")
            print(f"    |R_geom| = {r['amp_geom']:.4f}   (exp[-q*gamma] only)")
            print(f"    arg mod 360 = {r['phase_mod_deg']:.1f} deg")
            print(f"    arg unwrap  = {r['phase_unwrapped_deg']:.1f} deg")
            print(f"    omega_R candidates (k, value):")
            for k, w in cand:
                marker = " <-- true" if abs(w - wR) < 0.001 else ""
                print(f"      k={k:+d}: {w:.4f}/M{marker}")
            print()


if __name__ == '__main__':
    # Verify Schwarzschild
    c = critical_exponents_polar(0.0)
    assert abs(c['gamma'] - pi) < 1e-5
    assert abs(c['tau'] - 3*sqrt(3)*pi) < 1e-3
    assert abs(c['delta'] - pi) < 1e-5
    print("Schwarzschild check: PASS")
    
    # Demo
    print_ratios(0.7, m=2)
    print_ratios(0.0, m=2)
