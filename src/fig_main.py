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
    critical_exponents_polar, qnm_l2m2_fundamental,
    subring_ratio, subring_polarization_field
)

rcParams.update({'font.size': 11, 'axes.labelsize': 12,
                 'legend.fontsize': 9, 'figure.figsize': (13, 9), 'figure.dpi': 150,
                 'axes.grid': True, 'grid.alpha': 0.3})

C_PIX = '#2166AC'     # blue: R^{pix}
C_ALIGN = '#B2182B'   # red: R^{align}
C_GEO = '#FF7F00'     # orange: geometric/demagnification
C_GW = '#4DAF4A'      # green: GW/dynamical

def make_figure():
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
    
    # ============================================================
    # Panel (a): Critical parameters vs spin
    # ============================================================
    ax = axes[0, 0]
    spins = np.linspace(0.0, 0.99, 200)
    gamma_arr, tau_arr = [], []
    for a in spins:
        c = critical_exponents_polar(a)
        gamma_arr.append(c['gamma'])
        tau_arr.append(c['tau'])
    
    ax.plot(spins, gamma_arr, color=C_GEO, linewidth=2, label=r'$\gamma$ (Lyapunov)')
    ax_twin = ax.twinx()
    ax_twin.plot(spins, tau_arr, '--', color=C_GEO, linewidth=2, alpha=0.6,
                 label=r'$\tau\,/\,M$ (time delay)')
    ax.set_xlabel('Black hole spin $a/M$')
    ax.set_ylabel(r'Lyapunov exponent $\gamma$', color=C_GEO)
    ax_twin.set_ylabel(r'Time delay $\tau\,/\,M$', color=C_GEO)
    ax.set_title('(a)  Photon-ring critical parameters (per half-orbit)', fontweight='bold')
    ax.legend(loc='upper right')
    ax_twin.legend(loc='center right')
    
    # ============================================================
    # Panel (b): Phase decomposition — R^{pix} vs R^{align}
    # ============================================================
    ax = axes[0, 1]
    a_scan = np.linspace(0.0, 0.98, 100)
    
    pix_phase, align_phase, wr_tau_deg, md_deg = [], [], [], []
    for a in a_scan:
        rp = subring_ratio(a, m=2, q=1, aligned=False)
        ra = subring_ratio(a, m=2, q=1, aligned=True)
        wr_tau_deg.append(rp['omega_R_tau_q_deg'])
        md_deg.append(rp['m_delta_q_deg'])
        pix_phase.append(rp['R_phase_deg'])
        align_phase.append(ra['R_phase_deg'])
    
    # Unwrap for smoothness
    pix_uw = np.unwrap(np.radians(pix_phase))
    align_uw = np.unwrap(np.radians(align_phase))
    wr_uw = np.unwrap(np.radians(wr_tau_deg))
    md_uw = np.unwrap(np.radians(md_deg))
    
    ax.plot(a_scan, np.degrees(wr_uw), '--', color=C_GW, linewidth=1.5, alpha=0.7,
            label=r'$\omega_R\tau$ (GW phase)')
    ax.plot(a_scan, np.degrees(md_uw), ':', color=C_GEO, linewidth=1.5, alpha=0.7,
            label=r'$m\delta$ (geometric)')
    ax.plot(a_scan, np.degrees(align_uw), '-', color=C_ALIGN, linewidth=2.5,
            label=r'$\arg R^{\rm align} = \omega_R\tau$')
    ax.plot(a_scan, np.degrees(pix_uw), '-', color=C_PIX, linewidth=2.5,
            label=r'$\arg R^{\rm pix} = \omega_R\tau - m\delta$')
    
    ax.set_xlabel('Black hole spin $a/M$')
    ax.set_ylabel('Phase advance per subring step (degrees)')
    ax.set_title('(b)  Two distinct observables: ' + r'$R^{\rm pix}$ vs $R^{\rm align}$', fontweight='bold')
    ax.legend(loc='upper left', ncol=2)
    
    # ============================================================
    # Panel (c): Magnitude |R| vs spin
    # ============================================================
    ax = axes[1, 0]
    mag_pix, mag_align, e_gamma_arr = [], [], []
    for a in a_scan:
        rp = subring_ratio(a, m=2, q=1, aligned=False)
        ra = subring_ratio(a, m=2, q=1, aligned=True)
        mag_pix.append(rp['R_mag'])
        mag_align.append(ra['R_mag'])
        e_gamma_arr.append(np.exp(rp['gamma_q']))
    
    ax.plot(a_scan, mag_pix, '-', color=C_PIX, linewidth=2,
            label=r'$|R| = e^{-\gamma + \omega_I\tau}$')
    ax.plot(a_scan, e_gamma_arr, '--', color='gray', linewidth=1.5, alpha=0.6,
            label=r'$e^\gamma$ (demagnification)')
    ax.axhline(y=1.0, color='gray', linestyle=':', alpha=0.3)
    ax.set_xlabel('Black hole spin $a/M$')
    ax.set_ylabel(r'$|R|$')
    ax.set_title('(c)  Subring ratio magnitude', fontweight='bold')
    ax.legend()
    ax.set_yscale('log')
    
    # ============================================================
    # Panel (d): R^{pix} and R^{align} convergence at a=0.7
    # ============================================================
    ax = axes[1, 1]
    a_demo = 0.7
    t_o = 20.0
    
    rp_pred = subring_ratio(a_demo, m=2, q=1, aligned=False)
    ra_pred = subring_ratio(a_demo, m=2, q=1, aligned=True)
    
    # Generate field and compute ratios
    ns, DeltaP = subring_polarization_field(t_o, a_demo, m=2, n_max=8)
    
    # R^{pix}: same φ
    ratios_pix = DeltaP[1:] / DeltaP[:-1]
    phases_pix = np.degrees(np.angle(ratios_pix)) % 360
    
    # R^{align}: φ → φ + δ 
    crit = critical_exponents_polar(a_demo)
    delta = crit['delta']
    ns_r = np.arange(1, 7)
    phases_align = []
    for n in ns_r:
        dp_n = DeltaP[n-1]
        dp_np1 = np.exp(-(n+1)*crit['gamma']) * np.exp(
            1j * (-(rp_pred['omega_R'] - 1j*rp_pred['omega_I']) * (t_o - (n+1)*crit['tau']) +
                  2 * ((0 + delta) - (n+1)*delta)))
        r = dp_np1 / dp_n
        phases_align.append(np.degrees(np.angle(r)) % 360)
    
    ax.plot(ns_r, phases_pix[:6], 'o', color=C_PIX, markersize=8, 
            label=r'$R^{\rm pix}$ numerical')
    ax.axhline(y=rp_pred['R_phase_deg'], color=C_PIX, linestyle='-', linewidth=2, alpha=0.7,
               label=f'predicted: {rp_pred["R_phase_deg"]:.1f}°')
    
    ax.plot(ns_r, phases_align, 's', color=C_ALIGN, markersize=8,
            label=r'$R^{\rm align}$ numerical')
    ax.axhline(y=ra_pred['R_phase_deg'], color=C_ALIGN, linestyle='-', linewidth=2, alpha=0.7,
               label=f'predicted: {ra_pred["R_phase_deg"]:.1f}°')
    
    ax.set_xlabel('Subring index $n$')
    ax.set_ylabel(r'$\arg R_{n,1}$ (degrees)')
    ax.set_title(f'(d)  Convergence: $R^{{\\rm pix}}$ vs $R^{{\\rm align}}$ ($a={a_demo}$)', fontweight='bold')
    ax.legend(fontsize=8, loc='center right')
    
    # Text box - precompute strings to avoid f-string backslash issues
    g_s = f'{crit["gamma"]:.2f}'
    t_s = f'{crit["tau"]:.1f}'
    d_s = f'{crit["delta_deg"]:.1f}'
    w_s = f'{rp_pred["omega_R"]:.4f}'
    wt_s = f'{rp_pred["omega_R_tau_q_deg"]:.1f}'
    md_s = f'{rp_pred["m_delta_q_deg"]:.1f}'
    
    textstr = (
        f'$a={a_demo}$  $m=2$  $q=1$\n'
        f'$\\gamma={g_s}$  $\\tau={t_s}M$\n'
        f'$\\delta={d_s}^\\circ$\n'
        f'$M\\omega_R={w_s}$\n'
        f'$\\omega_R\\tau={wt_s}^\\circ$\n'
        f'$m\\delta={md_s}^\\circ$'
    )
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.03, 0.97, textstr, transform=ax.transAxes, fontsize=7.5,
            verticalalignment='top', bbox=props, family='monospace')
    
    # ============================================================
    # Suptitle
    # ============================================================
    fig.suptitle(
        'Photon-Ring Subrings as Discrete Ringdown Samplers\n'
        + r'$R^{\rm pix}_{n,q} \to e^{-q\gamma}\,e^{iq(\omega\tau-m\delta)}$'
        + '     '
        + r'$R^{\rm align}_{n,q} \to e^{-q\gamma}\,e^{iq\omega\tau}$',
        fontsize=14, fontweight='bold', y=1.01
    )
    plt.tight_layout()
    
    outdir = '/home/zhangzelin/research/projects/subring-ringdown/figures'
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(f'{outdir}/subring_main.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[Figure saved to {outdir}/subring_main.png]")
    
    return fig


if __name__ == '__main__':
    make_figure()
