#!/usr/bin/env python3
"""
PPP JSONL Visualizer – Comprehensive Edition
Usage:  python3 ppp_viewer.py <path/to/_epochs.jsonl> [--show]

Outputs
-------
  <name>_viz.pdf          7-page diagnostic PDF
  <name>_viz_overview.png Quick-look overview PNG (page 1)

Pages
-----
  1. Position accuracy & geometry
  2. Pseudorange residuals & innovations  (+CP and Doppler histogram, residuals vs elev+CN0)
  3. Carrier phase innovations & ambiguities
  4. Doppler residuals & PVUW
  5. Atmosphere & clocks
  6. Signal quality  (C/N0 time series, C/N0 vs elevation, SV tracking strip, DOP zones)
  7. PPP diagnostics  (convergence curve, sigma-ratio, elevation-binned box plots)
"""

import sys, json, math
from collections import defaultdict

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse

matplotlib.rcParams.update({
    'font.size': 10, 'axes.titlesize': 10.5, 'axes.labelsize': 9.5,
    'legend.fontsize': 7.5, 'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5,
    'lines.linewidth': 1.3, 'axes.linewidth': 0.9, 'figure.dpi': 120,
})

# ─── I/O ────────────────────────────────────────────────────────────────────

def load(path):
    rows = []
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if s:
                rows.append(json.loads(s))
    return rows

def safe(rows, fn, default=float('nan')):
    out = []
    for r in rows:
        try:
            out.append(fn(r))
        except Exception:
            out.append(default)
    return np.array(out, dtype=float)

# ─── statistics ──────────────────────────────────────────────────────────────

def rmse(v):
    v = v[np.isfinite(v)]
    return float(np.sqrt(np.mean(v ** 2))) if len(v) else float('nan')

def bias(v):
    v = v[np.isfinite(v)]
    return float(np.mean(v)) if len(v) else float('nan')

def pct(v, p):
    v = np.abs(v[np.isfinite(v)])
    return float(np.percentile(v, p)) if len(v) else float('nan')

def ylim_pct(arrays, lo=1.5, hi=98.5, margin=0.12):
    vals = np.concatenate([a[np.isfinite(a)] for a in arrays if a is not None])
    if not len(vals):
        return -1.0, 1.0
    lo_v, hi_v = np.percentile(vals, lo), np.percentile(vals, hi)
    span = max(hi_v - lo_v, 1e-9)
    return lo_v - margin * span, hi_v + margin * span

# ─── plot helpers ─────────────────────────────────────────────────────────────

def grid(ax):
    ax.grid(True, alpha=0.25, lw=0.5)

def hline(ax, y=0.0):
    ax.axhline(y, color='k', lw=0.7, ls='--', alpha=0.5)

def cs_vlines(ax, sid, cs_dict):
    """Bright red solid vlines for cycle slips — always visible regardless of line colour."""
    for ep in cs_dict.get(sid, []):
        ax.axvline(ep, color='red', lw=2.0, ls='-', alpha=0.9, zorder=6)

def stat_box(ax, arr, unit='cm', loc='upper right'):
    r, b = rmse(arr), bias(arr)
    if math.isnan(r):
        return
    txt = f'RMSE = {r:.2f} {unit}\nbias  = {b:+.2f} {unit}'
    ha = 'right' if 'right' in loc else 'left'
    x  = 0.985 if ha == 'right' else 0.015
    ax.text(x, 0.975, txt, transform=ax.transAxes, va='top', ha=ha,
            fontsize=7.5, family='monospace',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.85, lw=0.5))

def pct_box(ax, arr, unit=''):
    q50, q75, q95 = pct(arr, 50), pct(arr, 75), pct(arr, 95)
    txt = f'|Q50|={q50:.2f}  |Q75|={q75:.2f}  |Q95|={q95:.2f} {unit}'
    ax.text(0.5, 0.01, txt, transform=ax.transAxes, va='bottom', ha='center',
            fontsize=7, family='monospace',
            bbox=dict(boxstyle='round,pad=0.25', fc='lightyellow', alpha=0.9, lw=0.5))

# ─── data builder ─────────────────────────────────────────────────────────────

def build(rows):
    n = len(rows)
    t = np.arange(n)
    D = {'t': t, 'n': n, 'has_err': 'err_enu' in rows[0].get('pos', {})}

    for key, fn, sc in [
        ('err_e',  lambda r: r['pos']['err_enu'][0],             100),
        ('err_n',  lambda r: r['pos']['err_enu'][1],             100),
        ('err_u',  lambda r: r['pos']['err_enu'][2],             100),
        ('sig_e',  lambda r: r['pos']['sigma_enu'][0],           100),
        ('sig_n',  lambda r: r['pos']['sigma_enu'][1],           100),
        ('sig_u',  lambda r: r['pos']['sigma_enu'][2],           100),
        ('n_sats', lambda r: r['n_sats'],                          1),
        ('hdop',   lambda r: r['geometry']['hdop'],                1),
        ('pdop',   lambda r: r['geometry']['pdop'],                1),
        ('vdop',   lambda r: r['geometry']['vdop'],                1),
        ('pvuw_pr',    lambda r: r['pvuw']['pr'],                  1),
        ('pvuw_phase', lambda r: r['pvuw']['phase'],               1),
        ('pvuw_dop',   lambda r: r['pvuw']['doppler'],             1),
        ('pvuw_gim',   lambda r: r['pvuw']['gim'],                 1),
        ('red_pr',     lambda r: r['redundancy']['pr'],            1),
        ('red_phase',  lambda r: r['redundancy']['phase'],         1),
        ('red_dop',    lambda r: r['redundancy']['doppler'],       1),
        ('red_gim',    lambda r: r['redundancy']['gim'],           1),
        ('tropo_zwd',  lambda r: r['tropo']['zwd_m'],             100),
        ('tropo_sig',  lambda r: r['tropo']['sigma_m'],           100),
        ('drift',      lambda r: r['clock']['drift_mps']['base'],   1),
    ]:
        D[key] = safe(rows, fn) * sc

    D['err_3d'] = np.sqrt(D['err_e']**2 + D['err_n']**2 + D['err_u']**2)

    clk_codes = sorted({k for r in rows for k in r.get('clock', {}).get('offsets_m', {})})
    D['clk_codes'] = clk_codes
    D['clk'] = {k: safe(rows, lambda r, k=k: r['clock']['offsets_m'][k]) for k in clk_codes}

    sat_ids  = sorted({s['id'] for r in rows for s in r.get('satellites', [])})
    prefixes = sorted({sid[:3] for sid in sat_ids})
    D['sat_ids'] = sat_ids
    D['prefixes'] = prefixes

    def sarr(key, sc=1.0):
        A = {sid: np.full(n, np.nan) for sid in sat_ids}
        for i, r in enumerate(rows):
            for s in r.get('satellites', []):
                v = s.get(key)
                if v is not None:
                    A[s['id']][i] = v * sc
        return A

    D['resid_pr']  = sarr('resid_pr_m',    100.0)
    D['resid_cp']  = sarr('resid_cp_m',   1000.0)
    D['resid_dop'] = sarr('resid_dop_mps',   1.0)
    D['innov_pr']  = sarr('innov_pr_m',    100.0)
    D['innov_cp']  = sarr('innov_cp_m',    100.0)
    D['iono']      = sarr('iono_tecu')
    D['iono_sig']  = sarr('iono_sigma_tecu')
    D['amb']       = sarr('amb_float_cyc')
    D['amb_sig']   = sarr('amb_sigma_cyc')
    D['elev']      = sarr('elev_deg')
    D['azim']      = sarr('azim_deg')
    D['cn0']       = sarr('cn0_dbhz')
    D['resid_gim'] = sarr('resid_gim_tecu')   # GIM post-fit residual [TECU]
    D['innov_gim'] = sarr('innov_gim_tecu')   # GIM pre-fit innovation [TECU]

    D['cs']  = defaultdict(list)
    D['rst'] = defaultdict(list)
    for i, r in enumerate(rows):
        for s in r.get('satellites', []):
            if s.get('cs'):    D['cs'][s['id']].append(i)
            if s.get('reset'): D['rst'][s['id']].append(i)

    prns   = sorted({sid[3:] for sid in sat_ids})
    cmap   = plt.cm.tab10
    pc     = {p: cmap(i / max(len(prns) - 1, 1)) for i, p in enumerate(prns)}
    # Use solid and dashed ONLY to distinguish same-PRN different-signal (e.g. G1C05 vs G5Q05)
    ls_map = {p: s for p, s in zip(prefixes, ['-', '--', ':', '-.'])}
    D['sat_col'] = lambda sid: pc[sid[3:]]
    D['sat_ls']  = lambda sid: ls_map.get(sid[:3], '-')
    D['ls_map']  = ls_map
    D['prefixes_ls_note'] = '  |  '.join(f'{p}={ls}' for p, ls in ls_map.items())
    return D

# ─── reusable panel builders ──────────────────────────────────────────────────

def sat_series(ax, t, arr_map, sat_ids, prefix, D,
               ylabel='', title='', skip_n=0, ylim_lo=1.5, ylim_hi=98.5):
    sats = [s for s in sat_ids if s.startswith(prefix)]
    for sid in sats:
        a = arr_map[sid].copy()
        if skip_n > 0:
            a[:skip_n] = np.nan
        c = D['sat_col'](sid)
        ax.plot(t, a, color=c, ls=D['sat_ls'](sid), lw=1.2, label=sid[3:], alpha=0.9)
        cs_vlines(ax, sid, D['cs'])
    hline(ax)
    ax.set_title(f'{title}  [{prefix}]')
    ax.set_xlabel('Epoch'); ax.set_ylabel(ylabel)
    ax.legend(ncol=5, loc='upper right')
    grid(ax)
    arrs = [arr_map[s] for s in sats]
    if arrs:
        lo, hi = ylim_pct(arrs, ylim_lo, ylim_hi)
        ax.set_ylim(lo, hi)
    flat = np.concatenate([arr_map[s][np.isfinite(arr_map[s])] for s in sats]) if sats else np.array([])
    if len(flat):
        pct_box(ax, flat, ylabel.split('[')[-1].rstrip(']') if '[' in ylabel else '')

def scatter_vs_2var(fig, ax_el, ax_cn, elev_map, val_map, cn0_map, sat_ids,
                    title_el='', title_cn='', ylabel='', skip_n=0,
                    ylim_lo=2.0, ylim_hi=98.0):
    """
    Two complementary scatter panels that together give a 3-variable view:
      Left:  X=elevation,  Y=residual,  colour=C/N0   → how residuals change with geometry
      Right: X=C/N0,       Y=residual,  colour=elev   → how residuals change with signal strength
    """
    has_cn0 = any(np.any(np.isfinite(cn0_map[s])) for s in sat_ids)

    all_vals = []
    sc1 = sc2 = None
    for sid in sat_ids:
        el = elev_map[sid]; v = val_map[sid]; cn = cn0_map[sid]
        if skip_n > 0:
            v = v.copy(); v[:skip_n] = np.nan
        mask = np.isfinite(el) & np.isfinite(v)
        if has_cn0:
            mask &= np.isfinite(cn)
        if not mask.any():
            continue
        all_vals.extend(v[mask].tolist())
        cn_vals = cn[mask] if has_cn0 else np.full(mask.sum(), 35.0)
        sc1 = ax_el.scatter(el[mask], v[mask], s=5, c=cn_vals,
                            cmap='viridis', vmin=20, vmax=52, alpha=0.55)
        if has_cn0:
            sc2 = ax_cn.scatter(cn[mask], v[mask], s=5, c=el[mask],
                                cmap='RdYlGn', vmin=0, vmax=90, alpha=0.55)

    if sc1:
        cb1 = fig.colorbar(sc1, ax=ax_el, pad=0.02, fraction=0.04)
        cb1.set_label('C/N₀ [dBHz]', fontsize=8)
    if sc2:
        cb2 = fig.colorbar(sc2, ax=ax_cn, pad=0.02, fraction=0.04)
        cb2.set_label('Elevation [°]', fontsize=8)

    hline(ax_el); hline(ax_cn)
    ax_el.set_xlabel('Elevation [°]'); ax_el.set_ylabel(ylabel)
    ax_el.set_title(title_el); ax_el.set_xlim(0, 90)
    ax_cn.set_xlabel('C/N₀ [dBHz]'); ax_cn.set_ylabel(ylabel)
    ax_cn.set_title(title_cn if has_cn0 else f'{title_cn}  [C/N₀ unavailable]')
    grid(ax_el); grid(ax_cn)

    if all_vals:
        arr = np.array(all_vals)
        lo = np.percentile(arr, ylim_lo); hi = np.percentile(arr, ylim_hi)
        span = max(hi - lo, 1e-9)
        for ax in [ax_el, ax_cn]:
            ax.set_ylim(lo - 0.12 * span, hi + 0.12 * span)

def residual_hist(ax, arr_map, sat_ids, title='', xlabel='', color='steelblue',
                  skip_n=0, bins=50):
    arrs = []
    for s in sat_ids:
        a = arr_map[s].copy()
        if skip_n > 0:
            a[:skip_n] = np.nan
        arrs.append(a[np.isfinite(a)])
    flat = np.concatenate(arrs) if arrs else np.array([])
    if not len(flat):
        return
    lo, hi = np.percentile(flat, 1), np.percentile(flat, 99)
    clipped = flat[(flat >= lo) & (flat <= hi)]
    ax.hist(clipped, bins=bins, color=color, edgecolor='none', alpha=0.75, density=True)
    mu, sg = np.mean(clipped), np.std(clipped)
    if sg > 0:
        xs = np.linspace(clipped.min(), clipped.max(), 400)
        ax.plot(xs, np.exp(-0.5 * ((xs - mu) / sg) ** 2) / (sg * np.sqrt(2 * np.pi)),
                color='crimson', lw=2.0, label=f'N(μ={mu:.2f}, σ={sg:.2f})')
    ax.axvline(mu, color='k', lw=1.2, ls='--', label=f'mean = {mu:.2f}')
    q50, q75, q95 = np.percentile(np.abs(clipped), [50, 75, 95])
    ax.text(0.97, 0.96,
            f'|Q50| = {q50:.2f}\n|Q75| = {q75:.2f}\n|Q95| = {q95:.2f}',
            transform=ax.transAxes, va='top', ha='right', fontsize=8, family='monospace',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.9, lw=0.5))
    ax.legend(); ax.set_title(title); ax.set_xlabel(xlabel)
    ax.set_ylabel('Density'); grid(ax)

# ─── PAGE 1 : Position accuracy & geometry ────────────────────────────────────

def page1(D):
    fig = plt.figure(figsize=(22, 12))
    fig.suptitle('Page 1 — Position Accuracy & Geometry', fontsize=12, fontweight='bold')
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.46, wspace=0.35)
    t = D['t']

    for j, (err, sig, lbl, col) in enumerate([
        (D['err_e'], D['sig_e'], 'East error [cm]',  'royalblue'),
        (D['err_n'], D['sig_n'], 'North error [cm]', 'darkorange'),
        (D['err_u'], D['sig_u'], 'Up error [cm]',    'seagreen'),
    ]):
        ax = fig.add_subplot(gs[0, j])
        if D['has_err']:
            ax.plot(t, err, color=col, lw=1.3)
            ax.fill_between(t, err - sig, err + sig, alpha=0.22, color=col, label='±1σ')
            hline(ax); stat_box(ax, err)
            ax.legend()
        else:
            ax.plot(t, sig, color=col, lw=1.3, label='±1σ'); ax.legend()
        ax.set_title(lbl); ax.set_xlabel('Epoch'); grid(ax)

    ax3 = fig.add_subplot(gs[0, 3])
    if D['has_err']:
        ax3.plot(t, D['err_3d'], color='crimson', lw=1.3)
        stat_box(ax3, D['err_3d'])
    ax3.set_title('3D Error [cm]'); ax3.set_xlabel('Epoch'); grid(ax3)

    ax_sc = fig.add_subplot(gs[1, 0])
    ax_dop = ax_sc.twinx()
    ax_sc.plot(t, D['n_sats'], color='black', lw=1.4, label='#sat')
    ax_dop.plot(t, D['pdop'], color='royalblue', lw=1.1,         label='PDOP')
    ax_dop.plot(t, D['hdop'], color='tomato',    lw=1.0, ls='--', label='HDOP')
    ax_dop.plot(t, D['vdop'], color='seagreen',  lw=1.0, ls=':',  label='VDOP')
    ax_sc.set_ylabel('#Satellites', fontsize=9)
    ax_dop.set_ylabel('DOP', color='royalblue', fontsize=9)
    ax_sc.set_title('Satellite Count & DOP'); ax_sc.set_xlabel('Epoch')
    h1, l1 = ax_sc.get_legend_handles_labels()
    h2, l2 = ax_dop.get_legend_handles_labels()
    ax_sc.legend(h1 + h2, l1 + l2, loc='upper right'); grid(ax_sc)

    ax_red = fig.add_subplot(gs[1, 1])
    for arr, lbl in [(D['red_pr'], 'PR'), (D['red_phase'], 'Phase'),
                     (D['red_dop'], 'Doppler'), (D['red_gim'], 'GIM')]:
        ax_red.plot(t, arr, lw=1.2, label=lbl)
    ax_red.set_title('Redundancy Number'); ax_red.set_xlabel('Epoch')
    ax_red.set_ylabel('Redundancy'); ax_red.legend(); grid(ax_red)

    ax_2d = fig.add_subplot(gs[1, 2])
    if D['has_err']:
        mask = np.isfinite(D['err_e']) & np.isfinite(D['err_n'])
        ee, nn = D['err_e'][mask], D['err_n'][mask]
        ax_2d.scatter(ee, nn, s=3, alpha=0.45, color='steelblue', label='epoch')
        if len(ee) > 4:
            cov = np.cov(ee, nn)
            vals, vecs = np.linalg.eigh(cov)
            idx = vals.argsort()[::-1]
            vals, vecs = vals[idx], vecs[:, idx]
            w, h = 2 * np.sqrt(5.991 * vals)
            angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            ell = Ellipse((ee.mean(), nn.mean()), w, h, angle=angle,
                          fc='none', ec='crimson', lw=2.0, ls='--', label='95% ellipse')
            ax_2d.add_patch(ell)
            cep50 = float(np.percentile(np.sqrt(ee**2 + nn**2), 50))
            ax_2d.text(0.03, 0.03, f'CEP₅₀ = {cep50:.1f} cm',
                       transform=ax_2d.transAxes, fontsize=8.5, va='bottom',
                       bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.85))
        ax_2d.legend()
    ax_2d.axhline(0, color='k', lw=0.6, ls='--', alpha=0.5)
    ax_2d.axvline(0, color='k', lw=0.6, ls='--', alpha=0.5)
    ax_2d.set_xlabel('East error [cm]'); ax_2d.set_ylabel('North error [cm]')
    ax_2d.set_title('Horizontal Error Scatter (with 95% error ellipse)')
    ax_2d.set_aspect('equal', 'datalim'); grid(ax_2d)

    ax_sky = fig.add_subplot(gs[1, 3], projection='polar')
    ax_sky.set_theta_zero_location('N'); ax_sky.set_theta_direction(-1)
    ax_sky.set_ylim(0, 90)
    ax_sky.set_yticks([15, 30, 45, 60, 75, 90])
    ax_sky.set_yticklabels(['75°', '60°', '45°', '30°', '15°', '0°'], fontsize=6.5)
    ax_sky.set_title('Sky Plot  (N=0°, clockwise)', pad=14)
    labeled = set()
    for sid in D['sat_ids']:
        el = D['elev'][sid]; az = D['azim'][sid]
        valid = np.isfinite(el) & np.isfinite(az)
        if not valid.any(): continue
        c = D['sat_col'](sid)
        ax_sky.scatter(np.deg2rad(az[valid]), 90 - el[valid], s=3, color=c, alpha=0.55)
        if sid not in labeled:
            first = np.argmax(valid)
            ax_sky.annotate(sid[3:], (np.deg2rad(az[first]), 90 - el[first]),
                            fontsize=6.5, color=c, ha='center', va='bottom', fontweight='bold')
            labeled.add(sid)

    fig.tight_layout()
    return fig

# ─── PAGE 2 : Pseudorange residuals & innovations ─────────────────────────────

def page2(D):
    fig = plt.figure(figsize=(22, 16))
    fig.suptitle('Page 2 — Pseudorange: Post-fit Residuals & Pre-fit Innovations',
                 fontsize=12, fontweight='bold')
    gs = gridspec.GridSpec(4, 4, figure=fig, hspace=0.50, wspace=0.35)
    t = D['t']; pfx = D['prefixes']

    for j, px in enumerate(pfx[:2]):
        ax = fig.add_subplot(gs[0, j * 2:(j + 1) * 2])
        sat_series(ax, t, D['resid_pr'], D['sat_ids'], px, D,
                   ylabel='Residual [cm]', title='PR post-fit residual')

    for j, px in enumerate(pfx[:2]):
        ax = fig.add_subplot(gs[1, j * 2:(j + 1) * 2])
        sat_series(ax, t, D['innov_pr'], D['sat_ids'], px, D,
                   ylabel='Innovation [cm]', title='PR pre-fit innovation',
                   skip_n=1, ylim_lo=5, ylim_hi=95)

    ax_el = fig.add_subplot(gs[2, 0:2])
    ax_cn = fig.add_subplot(gs[2, 2])
    scatter_vs_2var(fig, ax_el, ax_cn, D['elev'], D['resid_pr'], D['cn0'], D['sat_ids'],
                    title_el='PR residual vs Elevation  (colour = C/N₀)',
                    title_cn='PR residual vs C/N₀  (colour = elev)',
                    ylabel='Residual [cm]')

    ax_hi = fig.add_subplot(gs[2, 3])
    residual_hist(ax_hi, D['resid_pr'], D['sat_ids'], color='steelblue',
                  title='PR Residual Distribution', xlabel='Residual [cm]')

    # CP histogram + Doppler histogram in row 3
    ax_cp_hi = fig.add_subplot(gs[3, 0:2])
    residual_hist(ax_cp_hi, D['resid_cp'], D['sat_ids'], color='darkorange',
                  title='CP Post-fit Residual Distribution  (mm scale — filter absorbs most into ambiguity)',
                  xlabel='Residual [mm]')

    ax_dop_hi = fig.add_subplot(gs[3, 2:4])
    residual_hist(ax_dop_hi, D['resid_dop'], D['sat_ids'], color='mediumorchid',
                  title='Doppler Residual Distribution', xlabel='Residual [m/s]')

    fig.tight_layout()
    return fig

# ─── PAGE 3 : Carrier phase innovations & ambiguities ─────────────────────────

def page3(D):
    fig = plt.figure(figsize=(22, 16))
    fig.suptitle('Page 3 — Carrier Phase Innovations & Float Ambiguities',
                 fontsize=12, fontweight='bold')
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.50, wspace=0.35)
    t = D['t']; pfx = D['prefixes']

    for j, px in enumerate(pfx[:2]):
        ax = fig.add_subplot(gs[0, j * 2:(j + 1) * 2])
        sat_series(ax, t, D['innov_cp'], D['sat_ids'], px, D,
                   ylabel='Innovation [cm]', title='CP pre-fit innovation  (ep 0-2 excluded)',
                   skip_n=3, ylim_lo=5, ylim_hi=95)

    # Float ambiguity — all same-signal satellites share one colour; G5Q is dashed by design
    for j, px in enumerate(pfx[:2]):
        ax = fig.add_subplot(gs[1, j * 2:(j + 1) * 2])
        sats_px = [s for s in D['sat_ids'] if s.startswith(px)]
        for sid in sats_px:
            c = D['sat_col'](sid)
            ax.plot(t, D['amb'][sid], color=c, ls=D['sat_ls'](sid),
                    lw=1.2, label=sid[3:], alpha=0.9)
            cs_vlines(ax, sid, D['cs'])   # bright red vertical lines
            for ep in D['rst'].get(sid, []):
                ax.axvline(ep, color='darkorange', lw=2.0, ls='-', alpha=0.85, zorder=5)
        n_cs  = sum(len(D['cs'].get(s,  [])) for s in sats_px)
        n_rst = sum(len(D['rst'].get(s, [])) for s in sats_px)
        ls_note = D['ls_map'].get(px, '-')
        ax.set_title(f'Float Ambiguity [cyc]  [{px}]  —  line style "{ls_note}" = {px} signal\n'
                     f'Red vline = cycle slip ({n_cs} total)  |  Orange vline = filter reset ({n_rst} total)')
        ax.set_xlabel('Epoch'); ax.set_ylabel('Ambiguity [cyc]')
        ax.legend(ncol=5, loc='upper right'); grid(ax)

    ax_sig = fig.add_subplot(gs[2, 0:2])
    for sid in D['sat_ids']:
        c = D['sat_col'](sid)
        ax_sig.plot(t, D['amb_sig'][sid], color=c, ls=D['sat_ls'](sid),
                    lw=1.0, label=sid[3:], alpha=0.85)
    lo, hi = ylim_pct([D['amb_sig'][s] for s in D['sat_ids']], 0, 98)
    ax_sig.set_ylim(max(lo, 0), hi)
    ax_sig.set_title('Ambiguity σ [cyc]  — convergence: all sats should stabilise to <0.05 cyc')
    ax_sig.set_xlabel('Epoch'); ax_sig.set_ylabel('σ [cyc]')
    ax_sig.legend(ncol=5, loc='upper right'); grid(ax_sig)

    ax_cp_hi = fig.add_subplot(gs[2, 2:4])
    innov_clip = {sid: D['innov_cp'][sid].copy() for sid in D['sat_ids']}
    for sid in innov_clip:
        innov_clip[sid][:5] = np.nan
    residual_hist(ax_cp_hi, innov_clip, D['sat_ids'], color='darkorange',
                  title='CP Innovation Distribution  (ep 0-4 excluded)', xlabel='Innovation [cm]',
                  skip_n=5)

    fig.tight_layout()
    return fig

# ─── PAGE 4 : Doppler & PVUW ──────────────────────────────────────────────────

def page4(D):
    fig = plt.figure(figsize=(22, 14))
    fig.suptitle('Page 4 — Doppler Residuals & Post-fit Variance of Unit Weight',
                 fontsize=12, fontweight='bold')
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.50, wspace=0.35)
    t = D['t']; pfx = D['prefixes']

    for j, px in enumerate(pfx[:2]):
        ax = fig.add_subplot(gs[0, j * 2:(j + 1) * 2])
        sat_series(ax, t, D['resid_dop'], D['sat_ids'], px, D,
                   ylabel='Residual [m/s]', title='Doppler post-fit residual')

    ax_el = fig.add_subplot(gs[1, 0:2])
    ax_cn = fig.add_subplot(gs[1, 2])
    scatter_vs_2var(fig, ax_el, ax_cn, D['elev'], D['resid_dop'], D['cn0'], D['sat_ids'],
                    title_el='Doppler residual vs Elevation  (colour = C/N₀)',
                    title_cn='Doppler residual vs C/N₀  (colour = elev)',
                    ylabel='Residual [m/s]')

    ax_dop_hi = fig.add_subplot(gs[1, 3])
    residual_hist(ax_dop_hi, D['resid_dop'], D['sat_ids'], color='mediumorchid',
                  title='Doppler Residual Distribution', xlabel='Residual [m/s]')

    ax_pvuw = fig.add_subplot(gs[2, :])
    t1 = t[1:]
    phase_plot = np.where(np.abs(D['pvuw_phase']) > 50, np.nan, D['pvuw_phase'])
    series_map = {'PR': D['pvuw_pr'][1:], 'Phase': phase_plot[1:],
                  'Doppler': D['pvuw_dop'][1:], 'GIM': D['pvuw_gim'][1:]}
    colors_pvuw = {'PR': 'royalblue', 'Phase': 'seagreen', 'Doppler': 'darkorange', 'GIM': 'crimson'}
    for lbl, arr in series_map.items():
        ax_pvuw.plot(t1, arr, lw=1.5, color=colors_pvuw[lbl], label=lbl)
    ax_pvuw.axhline(1.0, color='k', lw=1.2, ls='--', label='ideal = 1')
    all_pvuw = np.concatenate([a[np.isfinite(a)] for a in series_map.values()])
    if len(all_pvuw):
        yhi = max(float(np.percentile(all_pvuw, 97)) * 1.3, 2.0)
        pr_arr = series_map['PR'][np.isfinite(series_map['PR'])]
        if len(pr_arr):
            for p, lbl2 in [(50, 'PR Q50'), (75, 'PR Q75')]:
                vp = float(np.percentile(pr_arr, p))
                ax_pvuw.axhline(vp, color='royalblue', lw=1.0, ls=':', alpha=0.8,
                                label=f'{lbl2} = {vp:.2f}')
        ax_pvuw.set_ylim(0, yhi)
    ax_pvuw.set_title('Post-Variance of Unit Weight  (epoch 0 excluded  |  ideal ≈ 1.0 means observations match noise model)')
    ax_pvuw.set_xlabel('Epoch'); ax_pvuw.legend(ncol=6); grid(ax_pvuw)

    fig.tight_layout()
    return fig

# ─── PAGE 5 : Atmosphere & clocks ─────────────────────────────────────────────

def page5(D):
    fig = plt.figure(figsize=(22, 18))
    fig.suptitle('Page 5 — Atmosphere & Clocks', fontsize=12, fontweight='bold')
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.50, wspace=0.35)
    t = D['t']

    ax_clk = fig.add_subplot(gs[0, 0:2])
    for k, v in D['clk'].items():
        ax_clk.plot(t, v, lw=1.3, label=k)
    ax_clk.set_title('Code Clock Offsets [m]'); ax_clk.set_xlabel('Epoch')
    ax_clk.set_ylabel('[m]'); ax_clk.legend(); grid(ax_clk)

    ax_drift = fig.add_subplot(gs[0, 2])
    ax_drift.plot(t, D['drift'], color='purple', lw=1.3)
    ax_drift.set_title('Clock Drift [m/s]'); ax_drift.set_xlabel('Epoch')
    ax_drift.set_ylabel('[m/s]'); grid(ax_drift)

    ax_tro = fig.add_subplot(gs[0, 3])
    ax_tro.plot(t, D['tropo_zwd'], color='teal', lw=1.3, label='ZWD')
    ax_tro.fill_between(t, D['tropo_zwd'] - D['tropo_sig'],
                        D['tropo_zwd'] + D['tropo_sig'],
                        alpha=0.3, color='teal', label='±1σ')
    ax_tro.set_title('Troposphere ZWD [cm]'); ax_tro.set_xlabel('Epoch')
    ax_tro.legend(); grid(ax_tro)

    first_pfx = D['prefixes'][0] if D['prefixes'] else ''
    ax_ion = fig.add_subplot(gs[1, 0:2])
    for sid in D['sat_ids']:
        if not sid.startswith(first_pfx): continue
        c = D['sat_col'](sid)
        ax_ion.plot(t, D['iono'][sid], color=c, lw=1.2, label=sid[3:], alpha=0.9)
        ax_ion.fill_between(t, D['iono'][sid] - D['iono_sig'][sid],
                            D['iono'][sid] + D['iono_sig'][sid],
                            alpha=0.12, color=c)
    ax_ion.set_title(f'Ionosphere TECU per PRN  [{first_pfx}]  (shaded band = ±1σ uncertainty)')
    ax_ion.set_xlabel('Epoch'); ax_ion.set_ylabel('TECU')
    ax_ion.legend(ncol=5, loc='upper right'); grid(ax_ion)

    ax_ie = fig.add_subplot(gs[1, 2:4])
    for sid in D['sat_ids']:
        if not sid.startswith(first_pfx): continue
        el = D['elev'][sid]; io = D['iono'][sid]
        mask = np.isfinite(el) & np.isfinite(io)
        if mask.any():
            ax_ie.scatter(el[mask], io[mask], s=5, color=D['sat_col'](sid),
                          alpha=0.65, label=sid[3:])
    ax_ie.set_xlabel('Elevation [°]'); ax_ie.set_ylabel('TECU')
    ax_ie.set_title(f'Ionosphere vs Elevation  [{first_pfx}]')
    ax_ie.set_xlim(0, 90); ax_ie.legend(ncol=5); grid(ax_ie)

    # ── Row 2: GIM pseudo-observable ─────────────────────────────────────────
    # Only one GIM state per physical PRN; all signals share it — plot per unique PRN
    # Use first-signal representative for each PRN
    prn_rep = {}   # PRN string → one satellite ID as representative
    for sid in D['sat_ids']:
        prn = sid[3:]
        if prn not in prn_rep:
            prn_rep[prn] = sid

    has_gim = any(np.any(np.isfinite(D['resid_gim'][s])) for s in D['sat_ids'])

    ax_gim_res = fig.add_subplot(gs[2, 0:2])
    ax_gim_inn = fig.add_subplot(gs[2, 2])
    ax_gim_hi  = fig.add_subplot(gs[2, 3])

    if has_gim:
        for prn, sid in prn_rep.items():
            c = D['sat_col'](sid)
            ax_gim_res.plot(t, D['resid_gim'][sid], color=c, lw=1.2, label=prn, alpha=0.9)
            ax_gim_inn.plot(t, D['innov_gim'][sid], color=c, lw=1.2, label=prn, alpha=0.9)
        residual_hist(ax_gim_hi, D['resid_gim'],
                      list(prn_rep.values()),   # one representative per PRN
                      color='mediumpurple',
                      title='GIM Residual Distribution', xlabel='Residual [TECU]')
        # add percentile ylim
        all_res = [D['resid_gim'][s] for s in prn_rep.values()]
        lo, hi = ylim_pct(all_res, 1, 99)
        ax_gim_res.set_ylim(lo, hi)
        all_inn = [D['innov_gim'][s] for s in prn_rep.values()]
        lo2, hi2 = ylim_pct(all_inn, 2, 98)
        ax_gim_inn.set_ylim(lo2, hi2)
    else:
        for ax in [ax_gim_res, ax_gim_inn, ax_gim_hi]:
            ax.text(0.5, 0.5, 'GIM residuals not in JSONL — re-run Java after EKF_PPP update',
                    transform=ax.transAxes, ha='center', va='center', fontsize=9, color='grey')

    hline(ax_gim_res); hline(ax_gim_inn)
    ax_gim_res.set_title('GIM Post-fit Residual per PRN [TECU]\n'
                         '(= estimated iono − GIM prior; small = GIM agrees with filter)')
    ax_gim_res.set_xlabel('Epoch'); ax_gim_res.set_ylabel('Residual [TECU]')
    ax_gim_res.legend(ncol=5, loc='upper right'); grid(ax_gim_res)

    ax_gim_inn.set_title('GIM Pre-fit Innovation per PRN [TECU]\n'
                         '(= GIM prior − predicted iono; large = GIM disagrees with filter prediction)')
    ax_gim_inn.set_xlabel('Epoch'); ax_gim_inn.set_ylabel('Innovation [TECU]')
    ax_gim_inn.legend(ncol=5, loc='upper right'); grid(ax_gim_inn)

    fig.tight_layout()
    return fig

# ─── PAGE 6 : Signal quality ──────────────────────────────────────────────────

def page6(D):
    fig = plt.figure(figsize=(22, 15))
    fig.suptitle('Page 6 — Signal Quality: C/N₀ & Satellite Tracking Visibility',
                 fontsize=12, fontweight='bold')
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.48, wspace=0.30)
    t = D['t']

    has_cn0 = any(np.any(np.isfinite(D['cn0'][s])) for s in D['sat_ids'])

    # ── Row 0: C/N0 time series (full width) ─────────────────────────────────
    ax_cn0 = fig.add_subplot(gs[0, :])
    if has_cn0:
        for sid in D['sat_ids']:
            c = D['sat_col'](sid)
            ax_cn0.plot(t, D['cn0'][sid], color=c, ls=D['sat_ls'](sid),
                        lw=1.2, label=sid[3:], alpha=0.9)
        ax_cn0.axhline(35, color='black',  lw=1.2, ls='--', alpha=0.7,
                       label='35 dBHz (good threshold)')
        ax_cn0.axhline(20, color='tomato', lw=1.2, ls='--', alpha=0.7,
                       label='20 dBHz (minimum tracking)')
        ax_cn0.set_ylim(0, 60)
        ax_cn0.fill_between(t, 0, 20, alpha=0.06, color='tomato')
        ax_cn0.fill_between(t, 20, 35, alpha=0.06, color='goldenrod')
        ax_cn0.fill_between(t, 35, 60, alpha=0.06, color='green')
    else:
        ax_cn0.text(0.5, 0.5, 'C/N₀ not present — re-run Java after the cn0_dbhz fix',
                    transform=ax_cn0.transAxes, ha='center', va='center',
                    fontsize=12, color='grey')
    ax_cn0.set_title('C/N₀ per Satellite [dBHz]  — higher = stronger signal\n'
                     f'Line styles: {D["prefixes_ls_note"]}')
    ax_cn0.set_xlabel('Epoch'); ax_cn0.set_ylabel('C/N₀ [dBHz]')
    ax_cn0.legend(ncol=8, loc='lower right'); grid(ax_cn0)

    # ── Row 1: C/N0 vs elevation scatter + 5°-bin median ─────────────────────
    ax_ce = fig.add_subplot(gs[1, 0])
    if has_cn0:
        all_el, all_cn = [], []
        for sid in D['sat_ids']:
            el = D['elev'][sid]; cn = D['cn0'][sid]
            mask = np.isfinite(el) & np.isfinite(cn) & (cn > 0)
            if mask.any():
                ax_ce.scatter(el[mask], cn[mask], s=4, color=D['sat_col'](sid), alpha=0.4)
                all_el.extend(el[mask].tolist()); all_cn.extend(cn[mask].tolist())
        all_el = np.array(all_el); all_cn = np.array(all_cn)
        bin_mids, bin_meds = [], []
        for lo in range(0, 90, 5):
            m = (all_el >= lo) & (all_el < lo + 5) & np.isfinite(all_cn)
            if m.sum() >= 3:
                bin_mids.append(lo + 2.5)
                bin_meds.append(np.median(all_cn[m]))
        if bin_mids:
            ax_ce.plot(bin_mids, bin_meds, 'k-o', lw=2.0, ms=5, zorder=5,
                       label='5° bin median')
        ax_ce.axhline(35, color='black',  lw=1.2, ls='--', alpha=0.6)
        ax_ce.axhline(20, color='tomato', lw=1.2, ls='--', alpha=0.6)
        ax_ce.set_ylim(0, 60)
        ax_ce.legend()
    ax_ce.set_title('C/N₀ vs Elevation\n(good antennas: C/N₀ rises smoothly with elevation)')
    ax_ce.set_xlabel('Elevation [°]'); ax_ce.set_ylabel('C/N₀ [dBHz]')
    ax_ce.set_xlim(0, 90); grid(ax_ce)

    # ── Row 1: DOP with colour-zone shading ───────────────────────────────────
    ax_dop = fig.add_subplot(gs[1, 1])
    ymax = max(float(np.nanmax(D['pdop'])), 6.0) * 1.1
    for lo_z, hi_z, col, lbl in [(0, 2, '#d4edda', 'Excellent (<2)'),
                                  (2, 4, '#fff3cd', 'Good (2–4)'),
                                  (4, 6, '#fde8c8', 'Moderate (4–6)'),
                                  (6, 20, '#f8d7da', 'Poor (>6)')]:
        if lo_z < ymax:
            ax_dop.axhspan(lo_z, min(hi_z, ymax), color=col, alpha=0.8, label=lbl)
    ax_dop.plot(t, D['pdop'], color='royalblue', lw=1.5,         label='PDOP', zorder=4)
    ax_dop.plot(t, D['hdop'], color='tomato',    lw=1.2, ls='--', label='HDOP', zorder=4)
    ax_dop.plot(t, D['vdop'], color='seagreen',  lw=1.2, ls=':',  label='VDOP', zorder=4)
    ax_dop.set_ylim(0, ymax)
    ax_dop.set_title('DOP with Quality Zones')
    ax_dop.set_xlabel('Epoch'); ax_dop.set_ylabel('DOP')
    ax_dop.legend(ncol=2, fontsize=7.5); grid(ax_dop)

    # ── Row 2: SV visibility / tracking strip (full width) ────────────────────
    ax_sv = fig.add_subplot(gs[2, :])
    sat_ids_sorted = sorted(D['sat_ids'])
    cmap_elev = plt.cm.RdYlGn

    for sid in sat_ids_sorted:
        y = sat_ids_sorted.index(sid)
        el = D['elev'][sid]
        tracked = np.isfinite(el)
        seg_start = None
        for ep in range(len(t)):
            at_end = ep == len(t) - 1
            if tracked[ep] and seg_start is None:
                seg_start = ep
            elif (not tracked[ep] or at_end) and seg_start is not None:
                seg_end = ep + 1 if (at_end and tracked[ep]) else ep
                seg_el = np.nanmean(el[seg_start:seg_end])
                color = cmap_elev(np.clip(seg_el / 90.0, 0, 1))
                ax_sv.broken_barh([(seg_start, seg_end - seg_start)], (y - 0.42, 0.84),
                                  facecolors=[color], edgecolors='none', alpha=0.9)
                seg_start = None
        # Cycle slip markers: tall red vertical lines, clearly visible
        for ep in D['cs'].get(sid, []):
            ax_sv.vlines(ep, y - 0.48, y + 0.48, colors='crimson', lw=3.0, zorder=6)

    ax_sv.set_yticks(range(len(sat_ids_sorted)))
    ax_sv.set_yticklabels(sat_ids_sorted, fontsize=8)
    ax_sv.set_xlim(0, len(t) - 1)
    ax_sv.set_xlabel('Epoch', fontsize=9)
    ax_sv.set_title('SV Tracking Timeline\n'
                    'Colour = mean elevation during segment  (red = low elevation, green = high elevation)\n'
                    'Crimson vertical bar = cycle slip detected', fontsize=10)

    sm = plt.cm.ScalarMappable(cmap=cmap_elev, norm=plt.Normalize(0, 90))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_sv, orientation='vertical', pad=0.01, fraction=0.015)
    cbar.set_label('Elevation [°]', fontsize=8)
    cbar.ax.tick_params(labelsize=7.5)
    grid(ax_sv)

    fig.tight_layout()
    return fig

# ─── PAGE 7 : PPP diagnostics ─────────────────────────────────────────────────

def page7(D):
    fig = plt.figure(figsize=(22, 15))
    fig.suptitle('Page 7 — PPP Diagnostics: Convergence, σ-Ratio & Elevation-Binned Residuals',
                 fontsize=12, fontweight='bold')
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.34)
    t = D['t']

    # ── Row 0: Convergence curve ──────────────────────────────────────────────
    ax_conv = fig.add_subplot(gs[0, :])
    if D['has_err']:
        err3d_m = D['err_3d'] / 100.0
        err_h_m = np.sqrt(D['err_e']**2 + D['err_n']**2) / 100.0
        err_u_m = np.abs(D['err_u']) / 100.0

        # Threshold bands (coloured background zones, log-scale friendly)
        ax_conv.axhspan(0.001, 0.1,  color='#d4edda', alpha=0.5, zorder=0)   # green  < 10 cm
        ax_conv.axhspan(0.1,   0.3,  color='#fff3cd', alpha=0.5, zorder=0)   # yellow 10-30 cm
        ax_conv.axhspan(0.3,   1.0,  color='#fde8c8', alpha=0.5, zorder=0)   # orange 30-100 cm
        ax_conv.axhspan(1.0,   1000, color='#f8d7da', alpha=0.4, zorder=0)   # red    > 1 m

        # Threshold lines with direct text labels (no legend entry)
        for val, lbl in [(0.1, '0.1 m'), (0.3, '0.3 m'), (1.0, '1.0 m')]:
            ax_conv.axhline(val, color='grey', lw=1.0, ls='--', zorder=1)
            ax_conv.text(len(t) * 0.01, val * 1.08, lbl, fontsize=8.5, color='dimgrey', va='bottom')

        # Data lines — 3D prominent, H and U secondary
        ax_conv.semilogy(t, err_u_m,  color='tomato',   lw=1.5, ls=':', alpha=0.75, label='Up')
        ax_conv.semilogy(t, err_h_m,  color='steelblue',lw=1.5, ls='--', alpha=0.75, label='Horizontal')
        ax_conv.semilogy(t, err3d_m,  color='navy',     lw=2.5, alpha=1.0, label='3D (primary)')

        ax_conv.set_ylim(1e-3, max(float(np.nanmax(err3d_m[np.isfinite(err3d_m)])) * 2, 2.0))
    else:
        ax_conv.text(0.5, 0.5, 'No reference position — err_enu not available',
                     transform=ax_conv.transAxes, ha='center', va='center',
                     fontsize=12, color='grey')
    ax_conv.set_title('Convergence Curve — position error vs epoch  (log scale)\n'
                      'Green zone < 0.1 m  |  Yellow 0.1–0.3 m  |  Orange 0.3–1.0 m  |  Red > 1.0 m',
                      fontsize=10)
    ax_conv.set_xlabel('Epoch'); ax_conv.set_ylabel('Error [m]  (log scale)')
    ax_conv.legend(fontsize=9, loc='upper right'); grid(ax_conv)

    # ── Row 1: σ-ratio per component ─────────────────────────────────────────
    if D['has_err']:
        for j, (err, sig, lbl, col) in enumerate([
            (D['err_e'], D['sig_e'], 'East',  'royalblue'),
            (D['err_n'], D['sig_n'], 'North', 'darkorange'),
            (D['err_u'], D['sig_u'], 'Up',    'seagreen'),
        ]):
            ax = fig.add_subplot(gs[1, j])
            ratio = np.abs(err) / np.clip(sig, 1e-9, None)
            ax.plot(t, ratio, color=col, lw=1.3, alpha=0.9)
            ax.axhline(1.0, color='black',  lw=1.5, ls='--', label='ideal = 1  (well-calibrated)')
            ax.axhline(2.0, color='tomato', lw=1.2, ls=':',  label='2σ bound')
            lo, hi = ylim_pct([ratio], 0, 98)
            ax.set_ylim(0, max(hi, 3.0))
            within = np.sum(ratio[np.isfinite(ratio)] <= 2.0) / max(np.sum(np.isfinite(ratio)), 1)
            ax.text(0.03, 0.97, f'{within:.1%} of epochs within 2σ',
                    transform=ax.transAxes, va='top', fontsize=8.5, family='monospace',
                    bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.88, lw=0.5))
            ax.set_title(f'σ-ratio  |{lbl} error| / σ\n'
                         f'> 1 → filter overconfident  |  < 1 → underconfident')
            ax.set_xlabel('Epoch'); ax.set_ylabel('|error| / σ')
            ax.legend(fontsize=8); grid(ax)
    else:
        for j in range(3):
            ax = fig.add_subplot(gs[1, j])
            ax.text(0.5, 0.5, 'No reference position available',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=11, color='grey')
            ax.set_title('σ-ratio  (unavailable)')

    # ── Row 2: Elevation-binned residual box plots ─────────────────────────────
    bin_edges = list(range(0, 95, 10))
    bin_labels = [f'{lo}–{lo+10}°' for lo in bin_edges[:-1]]

    for j, (arr_map, ylabel, col, title) in enumerate([
        (D['resid_pr'],  'Residual [cm]',  'steelblue',    'PR post-fit  vs elevation bin'),
        (D['innov_pr'],  'Innovation [cm]','darkorange',   'PR pre-fit innovation  vs elevation bin'),
        (D['resid_dop'], 'Residual [m/s]', 'mediumorchid', 'Doppler post-fit  vs elevation bin'),
    ]):
        ax = fig.add_subplot(gs[2, j])
        boxdata = []
        for lo in bin_edges[:-1]:
            vals = []
            for sid in D['sat_ids']:
                el = D['elev'][sid]; v = arr_map[sid]
                mask = np.isfinite(el) & np.isfinite(v) & (el >= lo) & (el < lo + 10)
                vals.extend(v[mask].tolist())
            if len(vals) > 4:
                p1, p99 = np.percentile(vals, 1), np.percentile(vals, 99)
                clipped = [x for x in vals if p1 <= x <= p99]
            else:
                clipped = vals
            boxdata.append(clipped if clipped else [0])
        ax.boxplot(boxdata, positions=range(len(bin_labels)), widths=0.65,
                   patch_artist=True, showfliers=False,
                   boxprops=dict(facecolor=col, alpha=0.5),
                   medianprops=dict(color='crimson', lw=2.0),
                   whiskerprops=dict(lw=1.0, color='dimgrey'),
                   capprops=dict(lw=1.0, color='dimgrey'))
        ax.set_xticks(range(len(bin_labels)))
        ax.set_xticklabels(bin_labels, rotation=38, ha='right', fontsize=7.5)
        ax.axhline(0, color='k', lw=1.0, ls='--', alpha=0.6)
        ax.set_title(f'{title}\n(box = IQR, red line = median, whiskers = 5th–95th pct)')
        ax.set_xlabel('Elevation bin')
        ax.set_ylabel(ylabel); grid(ax)

    fig.tight_layout()
    return fig

# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    args = sys.argv[1:]
    path = next((a for a in args if not a.startswith('--')), None)
    show = '--show' in args
    if not path:
        print("Usage: python3 ppp_viewer.py <epochs.jsonl> [--show]")
        sys.exit(1)

    rows = load(path)
    if not rows:
        print("No data found.")
        sys.exit(1)

    n_sats_unique = len({s['id'] for r in rows for s in r.get('satellites', [])})
    print(f"Loaded {len(rows)} epochs, {n_sats_unique} unique sat/signal IDs")
    D = build(rows)

    stem     = path.replace('.jsonl', '')
    pdf_path = stem + '_viz.pdf'
    png_path = stem + '_viz_overview.png'

    pages = [page1, page2, page3, page4, page5, page6, page7]

    page_figs = []
    with PdfPages(pdf_path) as pdf:
        for i, fn in enumerate(pages, 1):
            print(f"  Rendering page {i}/{len(pages)} …")
            fig = fn(D)
            pdf.savefig(fig, bbox_inches='tight')
            page_figs.append(fig)

    page_figs[0].savefig(png_path, dpi=150, bbox_inches='tight')
    print(f"Saved PDF     → {pdf_path}")
    print(f"Saved PNG     → {png_path}")

    if show:
        plt.show()
    else:
        for fig in page_figs:
            plt.close(fig)


if __name__ == '__main__':
    main()
