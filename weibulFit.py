import json
import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def weibull(t, M_inf, tau, beta):
    return M_inf * (1 - np.exp(-(t / tau) ** beta))


def load_and_prepare(filepath):
    """Load one JSON result, exclude trapped drug from the denominator,
    and return t, M(t)_effective, plus metadata."""
    with open(filepath, 'r') as f:
        data = json.load(f)

    info = data['dataInfo']
    ndrug = float(info['ndrug'])
    trapped = float(info.get('trapped', 0))
    effective_ndrug = max(ndrug - trapped, 1)  # avoid divide-by-zero

    t = np.array(data['releaseData']['MCstep'], dtype=float)
    cumRelease = np.array(data['releaseData']['cumDrugRelease'], dtype=float)

    # Re-normalize release % against only the drugs that are NOT trapped
    M_effective = (cumRelease / effective_ndrug) * 100
    M_effective = np.clip(M_effective, 0, 100)

    return {
        't': t,
        'M_effective': M_effective,
        'ndrug': ndrug,
        'trapped': trapped,
        'effective_ndrug': effective_ndrug,
        'polymerCon': info.get('polymerCon'),
        'polymerSize': info.get('polymerSize'),
        'drugCon': info.get('drugCon'),
    }


def fit_weibull_to(t, M):
    p0 = [max(M.max(), 1.0), t[len(t) // 4] or 10, 0.7]
    bounds = ([0, 0.1, 0.01], [100, t.max() * 5, 3])
    popt, _ = curve_fit(weibull, t, M, p0=p0, bounds=bounds, maxfev=20000)
    predicted = weibull(t, *popt)
    ss_res = np.sum((M - predicted) ** 2)
    ss_tot = np.sum((M - np.mean(M)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return popt, r_squared


def process_folder(json_folder, out_folder='./weibull_results/'):
    os.makedirs(out_folder, exist_ok=True)
    files = glob.glob(f'{json_folder}/*.json')
    print(f"Found {len(files)} JSON files")

    records = []
    loaded = []
    for fp in files:
        try:
            d = load_and_prepare(fp)
            d['file'] = fp
            loaded.append(d)
        except Exception as e:
            print(f"  Skipping {fp}: {e}")

    # Group by (polymerCon, drugCon) -> one chart per group, one line per polymerSize
    groups = {}
    for d in loaded:
        key = (d['polymerCon'], d['drugCon'])
        groups.setdefault(key, []).append(d)

    print(f"Found {len(groups)} (polymerCon, drugCon) combinations")

    for (pcon, dcon), items in groups.items():
        items = sorted(items, key=lambda x: x['polymerSize'])

        plt.figure(figsize=(10, 7))
        for d in items:
            t, M = d['t'], d['M_effective']

            try:
                popt, r2 = fit_weibull_to(t, M)
                M_inf, tau, beta = popt
                fit_ok = True
            except (RuntimeError, ValueError) as e:
                print(f"  Fit failed: polymerSize={d['polymerSize']}, "
                      f"polymerCon={pcon}, drugCon={dcon}: {e}")
                M_inf, tau, beta, r2 = np.nan, np.nan, np.nan, np.nan
                fit_ok = False

            records.append({
                'polymerCon': pcon, 'drugCon': dcon,
                'polymerSize': d['polymerSize'],
                'ndrug': d['ndrug'], 'trapped': d['trapped'],
                'effective_ndrug': d['effective_ndrug'],
                'M_inf': M_inf, 'tau': tau, 'beta': beta,
                'r_squared': r2, 'fit_success': fit_ok,
            })

            # Plot: thinned raw data as dots, fitted curve as solid line
            thin = slice(0, len(t), max(1, len(t) // 200))
            line, = plt.plot(t[thin], M[thin], '.', markersize=3, alpha=0.35)
            if fit_ok:
                t_smooth = np.linspace(t.min(), t.max(), 300)
                plt.plot(t_smooth, weibull(t_smooth, *popt), '-',
                         color=line.get_color(),
                         label=f'L={d["polymerSize"]}  '
                               f'(M∞={M_inf:.1f}, τ={tau:.1f}, β={beta:.2f})')
            else:
                plt.plot([], [], '-', color=line.get_color(),
                         label=f'L={d["polymerSize"]} (fit failed)')

        plt.xlabel('MC Steps (time)')
        plt.ylabel('Effective Drug Released (%)  [trapped drug excluded]')
        plt.title(f'Weibull Fit — Varying Polymer Chain Length\n'
                  f'polymerCon={pcon}, drugCon={dcon}')
        plt.legend(fontsize=9)
        plt.grid(True)
        plt.tight_layout()

        fname = f'{out_folder}weibull_polymerCon{pcon}_drugCon{dcon}.png'
        plt.savefig(fname, dpi=150)
        plt.close()
        print(f"  Saved {fname}")

    df = pd.DataFrame(records)
    df.to_csv(f'{out_folder}weibull_fit_results.csv', index=False)
    print(f"\nSaved fit results table: {out_folder}weibull_fit_results.csv")
    return df


if __name__ == '__main__':
    json_folder = './_data/varrypolySize/'
    df = process_folder(json_folder, out_folder='./_weibull_out/')

    failed = df[~df['fit_success']]
    if len(failed):
        print(f"\n{len(failed)} fits failed:")
        print(failed[['polymerCon', 'drugCon', 'polymerSize']].to_string(index=False))