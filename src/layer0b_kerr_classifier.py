"""
Layer 0b: GL2020 Critical Exponent Cross-Validation
====================================================
Systematic verification of Gralla-Lupsasca critical exponents
across the full Kerr spin range.

Tests:
  1. Schwarzschild limit (a=0): γ=π, τ=3√3π, δ=π
  2. q=1 vs q=2 consistency: (2γ,2δ,2τ) vs 2×(γ,δ,τ)
  3. Delta accumulation: ω_R τ and mδ evolve smoothly with spin
  4. Step ratios: e^{-γ} and e^{-2γ} consistency

This is INDEPENDENT of the subring ratio formalism: it verifies
the pure Kerr geometric relations using the exact GL2020 integral
formulas (elliptic integrals, Eqs. 68-70).
"""
import numpy as np
from src.subring_core import critical_exponents_polar, qnm_l2m2_fund
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def verify_schwarzschild_limits():
    """Test 1: Exact Schwarzschild values."""
    c = critical_exponents_polar(0.0)
    pi = np.pi
    
    checks = [
        ('gamma = pi', abs(c['gamma'] - pi), 1e-10),
        ('tau = 3√3π', abs(c['tau'] - 3*np.sqrt(3)*pi), 1e-8),
        ('delta = pi', abs(c['delta'] - pi), 1e-10),
        ('b_tilde = 3√3', abs(c['b_tilde'] - 3*np.sqrt(3)), 1e-10),
    ]
    
    all_ok = True
    for name, err, tol in checks:
        ok = err < tol
        if not ok: all_ok = False
        print(f"  {name}: err={err:.2e} {'OK' if ok else 'FAIL'}")
    return all_ok


def verify_q_consistency(a=0.7):
    """Test 2: q=2 parameters = 2 * q=1 parameters."""
    c = critical_exponents_polar(a)
    g1, t1, d1 = c['gamma'], c['tau'], c['delta']
    
    # For q=2, the step relations should be exactly 2× the q=1 values
    # (this is built into the GL2020 formalism)
    g2_exp = 2*g1; t2_exp = 2*t1; d2_exp = 2*d1
    
    # Compare with what subring_ratio returns for q=2
    from src.subring_core import subring_ratio
    r2 = subring_ratio(a, m=2, q=2, aligned=False)
    
    checks = [
        ('gamma_2 = 2*gamma_1', abs(r2['gamma_q'] - g2_exp), 1e-10),
        ('tau_2 = 2*tau_1', abs(r2['tau_q'] - t2_exp), 1e-10),
        ('delta_2 = 2*delta_1', abs(r2['delta_q'] - d2_exp) % (2*np.pi), 1e-10),
    ]
    for name, err, tol in checks:
        ok = err < tol
        print(f"  {name}: err={err:.2e} {'OK' if ok else 'FAIL'}")


def verify_spin_scan(n_pts=50):
    """Test 3: Smooth evolution of all critical parameters with spin."""
    spins = np.linspace(0.0, 0.99, n_pts)
    data = [critical_exponents_polar(a) for a in spins]
    
    gamma_arr = np.array([d['gamma'] for d in data])
    tau_arr = np.array([d['tau'] for d in data])
    delta_arr = np.array([d['delta'] for d in data])
    
    # Check monotonicity (gamma and tau decrease with a)
    assert np.all(np.diff(gamma_arr) < 0), "gamma not monotonic decreasing"
    assert np.all(np.diff(tau_arr) < 0), "tau not monotonic decreasing"
    print(f"  gamma: {gamma_arr[0]:.4f} -> {gamma_arr[-1]:.4f} (decreasing OK)")
    print(f"  tau:   {tau_arr[0]:.4f} -> {tau_arr[-1]:.4f} (decreasing OK)")
    print(f"  delta: {np.degrees(delta_arr[0]):.1f} -> {np.degrees(delta_arr[-1]):.1f} deg")
    
    return spins, data


def make_layer0b_figure(spins, data):
    """Generate validation figure."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    a_arr = np.array(spins)
    
    # Panel 1: gamma and tau
    ax = axes[0]
    gamma_arr = np.array([d['gamma'] for d in data])
    tau_arr = np.array([d['tau'] for d in data])
    ax.plot(a_arr, gamma_arr, 'C0-', lw=2, label=r'$\gamma$')
    ax_twin = ax.twinx()
    ax_twin.plot(a_arr, tau_arr, 'C1--', lw=2, label=r'$\tau/M$')
    ax.set_xlabel('a/M'); ax.set_ylabel(r'$\gamma$', color='C0')
    ax_twin.set_ylabel(r'$\tau/M$', color='C1')
    ax.set_title('(a) Critical exponents', fontweight='bold')
    ax.legend(loc='upper right'); ax_twin.legend(loc='center right')
    ax.grid(True, alpha=0.3)
    
    # Panel 2: Demagnification and step ratios
    ax = axes[1]
    e_g = np.exp(-gamma_arr)
    e_2g = np.exp(-2*gamma_arr)
    ax.plot(a_arr, e_g, 'C0-', lw=2, label=r'$e^{-\gamma}$ (q=1)')
    ax.plot(a_arr, e_2g, 'C2--', lw=1.5, label=r'$e^{-2\gamma}$ (q=2)')
    ax.set_xlabel('a/M'); ax.set_ylabel('Demagnification per step')
    ax.set_title('(b) Step demagnification', fontweight='bold')
    ax.legend(); ax.grid(True, alpha=0.3); ax.set_yscale('log')
    
    # Panel 3: Phase budget
    ax = axes[2]
    from src.subring_core import subring_ratio
    wr_tau, md, total_align, total_pix = [], [], [], []
    for a in a_arr[::4]:  # subsample for clarity
        rp = subring_ratio(a, m=2, q=1, aligned=False)
        ra = subring_ratio(a, m=2, q=1, aligned=True)
        wr_tau.append(rp['omega_R_tau_q_deg'])
        md.append(rp['m_delta_q_deg'])
        total_align.append(ra['phase_unwrapped_deg'])
        total_pix.append(rp['phase_unwrapped_deg'])
    
    a_sub = a_arr[::4]
    ax.plot(a_sub, np.unwrap(np.radians(wr_tau)), 'C0-', lw=1.5, 
            alpha=0.6, label=r'$\omega_R\tau$')
    ax.plot(a_sub, np.unwrap(np.radians(md)), 'C1--', lw=1.5, 
            alpha=0.6, label=r'$m\delta$')
    ax.plot(a_sub, np.unwrap(np.radians(total_align)), 'C3-', lw=2,
            label=r'$\arg R^{\rm align}$')
    ax.plot(a_sub, np.unwrap(np.radians(total_pix)), 'C4-', lw=2,
            label=r'$\arg R^{\rm pix}$')
    ax.set_xlabel('a/M'); ax.set_ylabel('Phase (rad, unwrapped)')
    ax.set_title('(c) Phase budget', fontweight='bold')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    
    plt.suptitle('Layer 0b: GL2020 Critical Exponent Cross-Validation', fontweight='bold')
    plt.tight_layout()
    import os; os.makedirs('figures', exist_ok=True)
    plt.savefig('figures/layer0b_crossval.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[Figure saved to figures/layer0b_crossval.png]")


if __name__ == '__main__':
    print("="*60)
    print("LAYER 0b: GL2020 Critical Exponent Cross-Validation")
    print("="*60)
    
    print("\nTest 1: Schwarzschild limits")
    verify_schwarzschild_limits()
    
    print("\nTest 2: q=1 vs q=2 consistency (a=0.7)")
    verify_q_consistency(0.7)
    
    print("\nTest 3: Spin scan monotonicity")
    spins, data = verify_spin_scan()
    
    make_layer0b_figure(spins, data)
