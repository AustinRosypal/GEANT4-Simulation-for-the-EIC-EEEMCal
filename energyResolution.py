import glob
import re
import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# -----------------------------
# Configuration
# -----------------------------
PWO_REGION = 0
SCI_REGION = 1

# Fiducial cuts (in mm) to ensure purity of material selection
PWO_RMIN = 130.0
PWO_RMAX = 220.0

SCI_RMIN = 320.0
SCI_RMAX = 560.0

# Trigger threshold for RMS/mean comparison plot only:
# discard events with measured energy <= 30% of beam energy
TRIGGER_FRACTION = 0.30

# -----------------------------
# Helper functions
# -----------------------------
# Bringing in the output, merged root files from the job
def extract_energy_from_filename(fname):
    m = re.search(r'ecal_(\d+(?:\.\d+)?)GeV_.*Apr23.*\.root', fname)
    if not m:
        return None
    return float(m.group(1))


def rms_resolution(arr):
    arr = np.asarray(arr)
    arr = arr[np.isfinite(arr)]

    if len(arr) < 2:
        return np.nan, np.nan, np.nan

    mu = np.mean(arr)
    sigma = np.std(arr, ddof=1)
    res = sigma / mu if mu != 0 else np.nan
    return mu, sigma, res


def apply_energy_trigger(arr, beam_energy_gev, trigger_fraction=TRIGGER_FRACTION):
    arr = np.asarray(arr)
    arr = arr[np.isfinite(arr)]

    threshold_mev = trigger_fraction * beam_energy_gev * 1000.0
    return arr[arr > threshold_mev]


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


def crystal_ball_resolution(arr, nbins=None):
    arr = np.asarray(arr)
    arr = arr[np.isfinite(arr)]

    if len(arr) < 50:
        return np.nan, np.nan, np.nan, None, None, None, None

    if nbins is None:
        nbins = max(20, min(50, int(np.sqrt(len(arr)))))

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

        A_fit, mu_fit, sigma_fit, alpha_fit, n_fit = popt
        sigma_fit = abs(sigma_fit)

        if not np.isfinite(mu_fit) or not np.isfinite(sigma_fit) or mu_fit <= 0:
            return np.nan, np.nan, np.nan, None, None, counts, edges

        xcurve = np.linspace(np.min(xfit), np.max(xfit), 500)
        ycurve = crystal_ball_left(xcurve, *popt)

        resolution = sigma_fit / mu_fit
        return mu_fit, sigma_fit, resolution, xcurve, ycurve, counts, edges

    except Exception:
        return np.nan, np.nan, np.nan, None, None, counts, edges


def finite_percent(arr):
    arr = np.array(arr, dtype=float)
    return 100.0 * arr


def print_stats(label, arr, rms_res, fit_res):
    n = len(arr)
    left = f"  {label:<18} N = {n:6d}"

    if np.isfinite(rms_res):
        left += f"   RMS/mean = {100*rms_res:.2f}%"
    else:
        left += "   RMS/mean = nan%"

    if np.isfinite(fit_res):
        left += f"   Crystal Ball = {100*fit_res:.2f}%"
    else:
        left += "   Crystal Ball = nan%"

    print(left)


def select_all_fiducial(total, region, r):
    """
    Apply the radial cuts for the two materials,
    but combined into a single 'all events' sample:
      - PbWO4-like events:   region == 0 and 130 < r < 220 mm
      - SciGlass-like events: region == 1 and 320 < r < 560 mm
    """
    total = np.asarray(total)
    region = np.asarray(region)
    r = np.asarray(r)

    mask = (
        ((region == PWO_REGION) & (r > PWO_RMIN) & (r < PWO_RMAX)) |
        ((region == SCI_REGION) & (r > SCI_RMIN) & (r < SCI_RMAX))
    )
    return total[mask]


# -----------------------------
# Bring in root files
# -----------------------------
all_files = glob.glob("ecal_*GeV_Apr23.root")
files = [f for f in all_files if extract_energy_from_filename(f) is not None]
files = sorted(files, key=lambda f: extract_energy_from_filename(f))

print("Matched files:", files)

if not files:
    raise RuntimeError("No matching ROOT files found in the current directory.")

# -----------------------------
# Storage for scan results
# -----------------------------
energies = []

# RMS/mean after trigger
res_all_rms = []
res_pwo_rms = []
res_sci_rms = []

# Crystal Ball fit resolutions
res_all_cb = []
res_pwo_cb = []
res_sci_cb = []

# -----------------------------
# Main scan loop
# -----------------------------
for fname in files:
    E = extract_energy_from_filename(fname)
    if E is None:
        continue

    with uproot.open(fname) as f:
        tree = f["events"]
        total = tree["totalEdep_MeV"].array(library="np")
        region = tree["initialRegion"].array(library="np")
        r = tree["impactR_mm"].array(library="np")

    print(f"\nProcessing file: {fname}")
    print("  initialRegion unique:", np.unique(region))
    print("  impactR_mm min/max:", np.min(r), np.max(r))
    print("  Total events in file:", len(total))

    # All events use the same fiducial region cuts as the two material samples
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

    print("  All-events fiducial count:", len(all_events))
    print("  PWO subset count:", len(pwo_events))
    print("  SciGlass subset count:", len(sci_events))

    # Apply trigger ONLY for RMS/mean comparison plot
    all_events_triggered = apply_energy_trigger(all_events, E)
    pwo_events_triggered = apply_energy_trigger(pwo_events, E)
    sci_events_triggered = apply_energy_trigger(sci_events, E)

    print("  Triggered all-events count:", len(all_events_triggered))
    print("  Triggered PWO count:", len(pwo_events_triggered))
    print("  Triggered SciGlass count:", len(sci_events_triggered))

    # RMS-based resolutions after trigger
    _, _, rall_rms = rms_resolution(all_events_triggered)
    _, _, rpwo_rms = rms_resolution(pwo_events_triggered)
    _, _, rsci_rms = rms_resolution(sci_events_triggered)

    # Crystal Ball fit-based resolutions
    all_events_cb = apply_energy_trigger(all_events, E)
    pwo_events_cb = apply_energy_trigger(pwo_events, E)
    sci_events_cb = apply_energy_trigger(sci_events, E)
    print("  CB-fit All:", len(all_events_cb))
    print("  CB-fit PWO:", len(pwo_events_cb))
    print("  CB-fit Sci:", len(sci_events_cb))
    _, _, rall_cb, _, _, _, _ = crystal_ball_resolution(all_events)
    _, _, rpwo_cb, _, _, _, _ = crystal_ball_resolution(pwo_events)
    _, _, rsci_cb, _, _, _, _ = crystal_ball_resolution(sci_events)

    # Store
    energies.append(E)

    res_all_rms.append(rall_rms)
    res_pwo_rms.append(rpwo_rms)
    res_sci_rms.append(rsci_rms)

    res_all_cb.append(rall_cb)
    res_pwo_cb.append(rpwo_cb)
    res_sci_cb.append(rsci_cb)

    print(fname)
    print_stats("All events:", all_events_triggered, rall_rms, rall_cb)
    print_stats("PWO fiducial:", pwo_events_triggered, rpwo_rms, rpwo_cb)
    print_stats("SciGlass fiducial:", sci_events_triggered, rsci_rms, rsci_cb)

# -----------------------------
# Convert to arrays
# -----------------------------
energies = np.array(energies, dtype=float)

res_all_rms = np.array(res_all_rms, dtype=float)
res_pwo_rms = np.array(res_pwo_rms, dtype=float)
res_sci_rms = np.array(res_sci_rms, dtype=float)

res_all_cb = np.array(res_all_cb, dtype=float)
res_pwo_cb = np.array(res_pwo_cb, dtype=float)
res_sci_cb = np.array(res_sci_cb, dtype=float)

print("\nenergies =", energies)
print("res_all_cb =", res_all_cb)
print("res_pwo_cb =", res_pwo_cb)
print("res_sci_cb =", res_sci_cb)
print("res_all_rms =", res_all_rms)
print("res_pwo_rms =", res_pwo_rms)
print("res_sci_rms =", res_sci_rms)

# -----------------------------
# Plot 1: Crystal Ball resolutions
# -----------------------------
plt.figure(figsize=(9, 6))

plt.plot(
    energies,
    finite_percent(res_pwo_cb),
    'o-',
    markersize=9,
    label=rf"PbWO$_{{4}}$ - Triggered, Crystal Ball"
)

plt.plot(
    energies,
    finite_percent(res_sci_cb),
    's-',
    markersize=9,
    label='SciGlass - Triggered, Crystal Ball'
)

plt.plot(
    energies,
    finite_percent(res_all_cb),
    '^-',
    markersize=9,
    label='All Events - Triggered, Crystal Ball'
)

plt.xlabel("Beam Energy (GeV)", fontsize=14)
plt.ylabel("Energy Resolution (%)", fontsize=14)
plt.title("EEEMCal Resolution vs. Energy Fit with Crystal Ball", fontsize=20)
plt.legend()
plt.tight_layout()
plt.show()

# -----------------------------
# Plot 2: RMS/mean resolutions after 30% trigger
# -----------------------------
plt.figure(figsize=(9, 6))

plt.plot(energies, finite_percent(res_pwo_rms), 'o-', markersize=9, label=rf"PbWO$_{{4}}$ Triggered, RMS/mean")
plt.plot(energies, finite_percent(res_sci_rms), 's-', markersize=9, label='SciGlass Triggered, RMS/mean')
plt.plot(energies, finite_percent(res_all_rms), '^-', markersize=9, label='All Events Triggered, RMS/mean')

plt.xlabel("Beam Energy (GeV)", fontsize=14)
plt.ylabel("Energy Resolution (%)", fontsize=14)
plt.title("EEEMCal Resolution vs. Energy (RMS/mean)", fontsize=20)
plt.legend()
plt.tight_layout()
plt.show()

