import glob
import re
import os
import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# -----------------------------
# Configuration
# -----------------------------
PWO_REGION = 0
SCI_REGION = 1

# Trigger security cuts (in mm)
PWO_RMIN = 130.0
PWO_RMAX = 220.0

SCI_RMIN = 320.0
SCI_RMAX = 560.0

OUTPUT_DIR = "histogram_plots_crystalball"
os.makedirs(OUTPUT_DIR, exist_ok=True)

MIN_BINS = 25
MAX_BINS = 60

USE_LOG_Y = True


# -----------------------------
# Helper functions
# -----------------------------
# Bringing in the output, merged root files from the job
def extract_energy_from_filename(fname):
    m = re.search(r'ecal_(\d+(?:\.\d+)?)GeV_.*Apr23.*\.root', fname)
    if not m:
        return None
    return float(m.group(1))


def rms_stats(arr):
    arr = np.asarray(arr)
    arr = arr[np.isfinite(arr)]

    if len(arr) < 2:
        return np.nan, np.nan, np.nan

    mu = np.mean(arr)
    sigma = np.std(arr, ddof=1)
    res = sigma / mu if mu != 0 else np.nan
    return mu, sigma, res


def crystal_ball_left(x, A, mu, sigma, alpha, n):
    x = np.asarray(x)
    sigma = np.abs(sigma)
    t = (x - mu) / sigma

    result = np.empty_like(t, dtype=float)

    core = t > -alpha
    result[core] = A * np.exp(-0.5 * t[core] ** 2)

    a = (n / alpha) ** n * np.exp(-0.5 * alpha ** 2)
    b = n / alpha - alpha

    tail = ~core
    result[tail] = A * a * (b - t[tail]) ** (-n)

    return result


def crystal_ball_fit(arr, nbins=None):
    arr = np.asarray(arr)
    arr = arr[np.isfinite(arr)]

    if len(arr) < 50:
        return np.nan, np.nan, np.nan, None, None, None, None

    if nbins is None:
        nbins = max(MIN_BINS, min(MAX_BINS, int(1.25 * np.sqrt(len(arr)))))

    counts, edges = np.histogram(arr, bins=nbins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    if np.sum(counts) == 0:
        return np.nan, np.nan, np.nan, None, None, counts, edges

    peak_idx = np.argmax(counts)
    mu0 = centers[peak_idx]

    sigma0 = np.std(arr, ddof=1)
    if not np.isfinite(sigma0) or sigma0 <= 0:
        return np.nan, np.nan, np.nan, None, None, counts, edges

    mask = (centers > 0.60 * mu0) & (centers < 1.08 * mu0)
    xfit = centers[mask]
    yfit = counts[mask]

    populated = yfit > 0
    xfit = xfit[populated]
    yfit = yfit[populated]

    if len(xfit) < 6:
        return np.nan, np.nan, np.nan, None, None, counts, edges

    p0 = [np.max(yfit), mu0, sigma0, 1.5, 3.0]
    bounds = (
        [0.0, 0.0, 1e-6, 0.2, 1.01],
        [np.inf, np.inf, np.inf, 10.0, 50.0]
    )

    try:
        popt, _ = curve_fit(
            crystal_ball_left,
            xfit,
            yfit,
            p0=p0,
            bounds=bounds,
            maxfev=50000
        )

        _, mu_fit, sigma_fit, _, _ = popt
        sigma_fit = abs(sigma_fit)

        if not np.isfinite(mu_fit) or not np.isfinite(sigma_fit) or mu_fit <= 0:
            return np.nan, np.nan, np.nan, None, None, counts, edges

        xcurve = np.linspace(np.min(xfit), np.max(xfit), 500)
        ycurve = crystal_ball_left(xcurve, *popt)

        resolution = sigma_fit / mu_fit
        return mu_fit, sigma_fit, resolution, xcurve, ycurve, counts, edges

    except Exception:
        return np.nan, np.nan, np.nan, None, None, counts, edges


def select_all_fiducial(total, region, r):
    mask = (
        ((region == PWO_REGION) & (r > PWO_RMIN) & (r < PWO_RMAX)) |
        ((region == SCI_REGION) & (r > SCI_RMIN) & (r < SCI_RMAX))
    )
    return np.asarray(total)[mask]


def plot_one_hist(ax, arr, title, color):
    arr = np.asarray(arr)
    arr = arr[np.isfinite(arr)]

    if len(arr) == 0:
        ax.text(0.5, 0.5, "No events", ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        ax.set_xlabel("Deposited Energy (MeV)", fontsize=18)
        ax.set_ylabel("Counts", fontsize=14)
        return

    mu_rms, sigma_rms, res_rms = rms_stats(arr)
    mu_fit, sigma_fit, res_fit, xcurve, ycurve, counts, edges = crystal_ball_fit(arr)

    nbins = max(MIN_BINS, min(MAX_BINS, int(1.25 * np.sqrt(len(arr)))))
    ax.hist(arr, bins=nbins, alpha=0.75, color=color)

    if xcurve is not None and ycurve is not None:
        ax.plot(xcurve, ycurve, 'k-', lw=2, label="Crystal Ball fit")

    if USE_LOG_Y:
        ax.set_yscale('log')
        ax.set_ylim(bottom=0.5)

    text_lines = [
        rf"$\mathrm{{N}}_{{events}} = {len(arr)}$",
        f"RMS/mean = {100*res_rms:.2f}%" if np.isfinite(res_rms) else "RMS/mean = nan"
    ]

    if np.isfinite(res_fit):
        text_lines.append(f"CB fit = {100*res_fit:.2f}%")
        text_lines.append(rf"$\mu_{{\mathrm{{CB}}}} = {mu_fit:.1f} MeV$")
        text_lines.append(rf"$\sigma_{{\mathrm{{CB}}}} = {sigma_fit:.1f} MeV$")
    else:
        text_lines.append("CB fit = failed")

    ax.text(
        0.67, 0.97,
        "\n".join(text_lines),
        transform=ax.transAxes,
        ha='right', va='top',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
    )

    ax.set_title(title, fontsize=16)
    ax.set_xlabel("Deposited Energy (MeV)", fontsize=14)
    ax.set_ylabel("Counts", fontsize=14)
    #ax.grid(True, which='both', alpha=0.3)

    if xcurve is not None and ycurve is not None:
        ax.legend(loc="upper left")


all_files = glob.glob("ecal_*GeV_Apr23.root")
files = [f for f in all_files if extract_energy_from_filename(f) is not None]
files = sorted(files, key=lambda f: extract_energy_from_filename(f))

print("Matched files:", files)

if not files:
    raise RuntimeError("No matching ROOT files found in the current directory.")

for fname in files:
    E = extract_energy_from_filename(fname)
    if E is None:
        continue

    with uproot.open(fname) as f:
        tree = f["events"]
        total = tree["totalEdep_MeV"].array(library="np")
        region = tree["initialRegion"].array(library="np")
        r = tree["impactR_mm"].array(library="np")

    all_events = select_all_fiducial(total, region, r)

    pwo_events = total[
        (region == PWO_REGION) &
        (r > PWO_RMIN) &
        (r < PWO_RMAX)
    ]

    sci_events = total[
        (region == SCI_REGION) &
        (r > SCI_RMIN) &
        (r < SCI_RMAX)
    ]

    print(f"\nProcessing {fname}")
    print(f"  Beam energy: {E} GeV")
    print(f"  All events fiducial: {len(all_events)}")
    print(f"  PWO fiducial: {len(pwo_events)}")
    print(f"  SciGlass fiducial: {len(sci_events)}")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    plot_one_hist(
        axes[0],
        all_events,
        f"All Triggered Events, {E} GeV",
        color="tab:blue"
    )

    plot_one_hist(
        axes[1],
        pwo_events,
        rf"PbWO$_{{4}}$ Triggered, {E} GeV",
        color="tab:orange"
    )

    plot_one_hist(
        axes[2],
        sci_events,
        f"SciGlass Triggered, {E} GeV",
        color="tab:green"
    )

    fig.suptitle(f"EEEMCal Deposited Energy Responses Fit with Crystal Ball ({E} GeV)", fontsize=20)
    plt.tight_layout()

    outfile = os.path.join(OUTPUT_DIR, f"histograms_crystalball_{E:g}GeV.png")
    plt.savefig(outfile, dpi=200)
    plt.show()
    plt.close(fig)

    print(f"  Saved plot to: {outfile}")

