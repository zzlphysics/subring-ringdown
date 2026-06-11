"""
Corrected demo figure for subring-ringdown formalism.
Shows both R^{pix} and R^{align} observables clearly separated.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os, sys
sys.path.insert(0, '/home/zhangzelin/research/projects/subring-ringdown')
from src.subring_core import (
    critical_exponents_polar, qnm_l2m2_fund, subring_ratio
)

rcParams.update({'font.size': 11, 'axes.labelsize': 12,
                 'legend.fontsize': 9, 'figure.figsize': (13, 9), 'figure.dpi': 150,
                 'axes.grid': True, 'grid.alpha': 0.3})

C_PIX   = '#2166AC'
C_ALIGN = '#B2182B'
C_GEO   = '#FF7F00'
C_GW    = '#4DAF4A'

def make_figure():
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
    
    # --- Panel (a): Critical parameters vs spin ---
    ax = axes[0, 0]
    spins = np.linspace(0.0, 0.99, 200)
    gamma_arr, tau_arr = [], []
    for a in spins:
        c = critical_exponents_polar(a)
        gamma_arr.append(c['gamma']); tau_arr.append(c['tau'])
    ax.plot(spins, gamma_arr, color=C_GEO, linewidth=2, label=r'$\gamma$ (Lyapunov)')
    ax_twin = ax.twinx()
    ax_twin.plot(spins, tau_arr, '--', color=C_GEO, linewidth=2, alpha=0.6,
                 label=r'$\tau\,/\,M$ (time delay)')
    ax.set_xlabel('Black hole spin $a/M$')
    ax.set_ylabel(r'Lyapunov exponent $\gamma$', color=C_GEO)
    ax_twin.set_ylabel(r'Time delay $\tau\,/\,M$', color=C_GEO)
    ax.set_title('(a)  Photon-ring critical parameters (per half-orbit)', fontweight='bold')
    ax.legend(loc='upper right'); ax_twin.legend(loc='center right')
    
    # --- Panel (b): Phase decomposition ---
    ax = axes[0, 1]
    a_scan = np.linspace(0.0, 0.98, 100)
    pix_ph, align_ph, wr_ph, md_ph = [], [], [], []
    for a in a_scan:
        rp = subring_ratio(a, m=2, q=1, aligned=False)
        ra = subring_ratio(a, m=2, q=1, aligned=True)
        wr_ph.append(rp['omega_R_tau_q_deg'])
        md_ph.append(rp['m_delta_q_deg'])
        pix_ph.append(rp['phase_mod_deg'])
        align_ph.append(ra['phase_mod_deg'])
    pix_uw   = np.unwrap(np.radians(pix_ph))
    align_uw = np.unwrap(np.radians(align_ph))
    wr_uw    = np.unwrap(np.radians(wr_ph))
    md_uw    = np.unwrap(np.radians(md_ph))
    
    ax.plot(a_scan, np.degrees(wr_uw), '--', color=C_GW, lw=1.5, alpha=0.7,
            label=r'$\omega_R\tau$ (GW phase)')
    ax.plot(a_scan, np.degrees(md_uw), ':', color=C_GEO, lw=1.5, alpha=0.7,
            label=r'$m\delta$ (geometric)')
    ax.plot(a_scan, np.degrees(align_uw), '-', color=C_ALIGN, lw=2.5,
            label=r'$\arg R^{\rm align} = \omega_R\tau$')
    ax.plot(a_scan, np.degrees(pix_uw), '-', color=C_PIX, lw=2.5,
            label=r'$\arg R^{\rm pix} = \omega_R\tau - m\delta$')
    ax.set_xlabel('Black hole spin $a/M$')
    ax.set_ylabel('Phase advance per subring step (degrees)')
    ax.set_title('(b)  ' + r'$R^{\rm pix}$ vs $R^{\rm align}$', fontweight='bold')
    ax.legend(loc='upper left', ncol=2)
    
    # --- Panel (c): Magnitude |R| ---
    ax = axes[1, 0]
    mag_raw, mag_geom = [], []
    for a in a_scan:
        rp = subring_ratio(a, m=2, q=1, aligned=False)
        mag_raw.append(rp['amp_raw'])
        mag_geom.append(rp['amp_geom'])
    ax.plot(a_scan, mag_raw, '-', color=C_PIX, lw=2,
            label=r'$|R_{\rm raw}| = e^{-\gamma + \omega_I\tau}$')
    ax.plot(a_scan, mag_geom, '--', color='gray', lw=1.5, alpha=0.6,
            label=r'$|R_{\rm geom}| = e^{-\gamma}$')
    ax.axhline(y=1.0, color='gray', linestyle=':', alpha=0.3)
    ax.set_xlabel('Black hole spin $a/M$')
    ax.set_ylabel(r'$|R|$')
    ax.set_title('(c)  Subring ratio magnitude', fontweight='bold')
    ax.legend(); ax.set_yscale('log')
    
    # --- Panel (d): Convergence at a=0.7 ---
    ax = axes[1, 1]
    a_demo = 0.7
    rp_p = subring_ratio(a_demo, m=2, q=1, aligned=False)
    ra_p = subring_ratio(a_demo, m=2, q=1, aligned=True)
    crit  = critical_exponents_polar(a_demo)
    
    # Generate synthetic field with monodromy algebra (Layer 1a demo)
    ns_test = np.arange(1, 8)
    DeltaP = []
    np.random.seed(42)
    I_emit = 1.0 + 0.05*(np.random.random()-0.5 + 1j*(np.random.random()-0.5))
    I_obs  = 1.0 + 0.05*(np.random.random()-0.5 + 1j*(np.random.random()-0.5))
    M_fac  = np.exp(-crit['gamma'] + 1j*(rp_p['omega_R']*crit['tau'] - 2*crit['delta']))
    for n in ns_test:
        dp = I_emit * (M_fac**n) * I_obs
        dp *= 1.0 + 0.02*(np.random.random()-0.5 + 1j*(np.random.random()-0.5))
        DeltaP.append(dp)
    DeltaP = np.array(DeltaP)
    
    # Ratios
    R_pix = DeltaP[1:] / DeltaP[:-1]
    R_align = R_pix * np.exp(1j * 2 * crit['delta'])  # correct geometric rotation
    
    ax.plot(ns_test[:-1], np.degrees(np.angle(R_pix)) % 360, 'o', color=C_PIX, ms=8,
            label=r'$R^{\rm pix}$ numerical')
    ax.axhline(y=rp_p['phase_mod_deg'], color=C_PIX, lw=2, alpha=0.7,
               label=f'pred: {rp_p["phase_mod_deg"]:.1f}$^\\circ$')
    ax.plot(ns_test[:-1], np.degrees(np.angle(R_align)) % 360, 's', color=C_ALIGN, ms=8,
            label=r'$R^{\rm align}$ numerical')
    ax.axhline(y=ra_p['phase_mod_deg'], color=C_ALIGN, lw=2, alpha=0.7,
               label=f'pred: {ra_p["phase_mod_deg"]:.1f}$^\\circ$')
    ax.set_xlabel('Subring index $n$')
    ax.set_ylabel(r'$\arg R_{n,1}$ (degrees)')
    ax.set_title(f'(d)  Convergence ($a={a_demo}$)', fontweight='bold')
    ax.legend(fontsize=8)
    
    # Text box
    g_s= f'{crit["gamma"]:.2f}'; t_s= f'{crit["tau"]:.1f}'
    d_s= f'{crit["delta_deg"]:.1f}'; w_s= f'{rp_p["omega_R"]:.4f}'
    wt_s=f'{rp_p["omega_R_tau_q_deg"]:.1f}'; md_s=f'{rp_p["m_delta_q_deg"]:.1f}'
    textstr = (f'$a={a_demo}$  $m=2$  $q=1$\n'
               f'$\\gamma={g_s}$  $\\tau={t_s}M$\n'
               f'$\\delta={d_s}^\\circ$\n$M\\omega_R={w_s}$\n'
               f'$\\omega_R\\tau={wt_s}^\\circ$\n$m\\delta={md_s}^\\circ$')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.03, 0.97, textstr, transform=ax.transAxes, fontsize=7.5,
            verticalalignment='top', bbox=props, family='monospace')
    
    fig.suptitle(
        'Photon-Ring Subrings as Discrete Ringdown Samplers\n'
        + r'$R^{\rm pix}_{n,q} \to e^{-q\gamma}\,e^{iq(\omega\tau-m\delta)}$'
        + '     '
        + r'$R^{\rm align}_{n,q} \to e^{-q\gamma}\,e^{iq\omega\tau}$',
        fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    os.makedirs('figures', exist_ok=True)
    plt.savefig('figures/subring_main.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[fig_main.py] Figure saved to figures/subring_main.png")

if __name__ == '__main__':
    make_figure()
