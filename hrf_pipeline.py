import argparse
import glob
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import loadmat, savemat
from scipy.optimize import fmin
from scipy.signal import find_peaks
from scipy.stats import gamma
from tqdm import tqdm


TIME_POINTS = np.arange(0, 26, 2)
TIME_FINE = np.arange(0, 24.01, 0.01)
GEM_COLORS = np.array(
    [
        [0, 0.4470, 0.7410],
        [0.8500, 0.3250, 0.0980],
        [0.9290, 0.6940, 0.1250],
        [0.4940, 0.1840, 0.5560],
        [0.4660, 0.6740, 0.1880],
        [0.3010, 0.7450, 0.9330],
        [0.6350, 0.0780, 0.1840],
    ]
)


@dataclass
class HrfMetrics:
    peak: float
    t_peak: float
    fwhm: float
    trough: float
    t_trough: float


def hrf_model(params: np.ndarray, t: np.ndarray) -> np.ndarray:
    return (
        params[4] * gamma.pdf(t, params[0], scale=params[1])
        - params[5] * gamma.pdf(t, params[2], scale=params[3])
        + params[6]
    )


def _penalty(params: np.ndarray, lb: np.ndarray, ub: np.ndarray) -> float:
    if np.any(params < lb) or np.any(params > ub):
        return 1e6
    return 0.0


def fit_hrf(y: np.ndarray, verbose: bool) -> np.ndarray:
    init_p = np.array([6, 1, 12, 1.2, 30, 15, y[0]], dtype=float)
    lb = np.array([4, 0.5, 4, 0.5, 10, 8, -0.2], dtype=float)
    ub = np.array([7, 1.5, 16, 2.5, 80, 50, 0.5], dtype=float)

    def obj(params: np.ndarray) -> float:
        return float(np.sum((hrf_model(params, TIME_POINTS) - y) ** 2) + _penalty(params, lb, ub))

    disp = bool(verbose)
    best_params = fmin(obj, init_p, disp=disp)
    return best_params


def _weighted_average(data: np.ndarray, cnt: np.ndarray, rows: Iterable[int]) -> np.ndarray:
    weighted = np.zeros(data.shape[2], dtype=float)
    total = 0.0
    for row in rows:
        for hemi in range(2):
            weight = float(cnt[row, hemi])
            weighted += data[row, hemi, :] * weight
            total += weight
    if total == 0:
        return weighted
    return weighted / total


def extract_regions(mat_path: str) -> Dict[str, np.ndarray]:
    mat = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    hrf_roi = mat["hrf_ROI"]
    data = np.asarray(hrf_roi.data, dtype=float)
    cnt = np.asarray(hrf_roi.cnt, dtype=float)

    v1 = _weighted_average(data, cnt, rows=[0, 1])
    v2 = _weighted_average(data, cnt, rows=[2, 3])
    v3 = _weighted_average(data, cnt, rows=[4, 5])
    v4 = _weighted_average(data, cnt, rows=[6])
    all_roi = np.asarray(hrf_roi.all, dtype=float).reshape(-1)

    return {"V1": v1, "V2": v2, "V3": v3, "V4": v4, "All": all_roi}


def compute_metrics(y_fine: np.ndarray, interval: float = 0.01) -> HrfMetrics:
    peaks, locs = find_peaks(y_fine)
    peak = float("nan")
    t_peak = float("nan")
    if len(peaks) > 0:
        peak = float(y_fine[peaks[0]])
        t_peak = float(peaks[0] * interval)
        if t_peak < 4 and len(peaks) > 1:
            peak = float(y_fine[peaks[1]])
            t_peak = float(peaks[1] * interval)

    fwhm = float(np.sum(y_fine > (peak / 2.0)) * interval) if np.isfinite(peak) else float("nan")

    trough = float("nan")
    t_trough = float("nan")
    troughs, tlocs = find_peaks(-y_fine)
    if len(troughs) > 0:
        trough = float(-y_fine[troughs[0]])
        t_trough = float(troughs[0] * interval)
        if trough < 6 and len(troughs) > 1:
            trough = float(-y_fine[troughs[1]])
            t_trough = float(troughs[1] * interval)

    return HrfMetrics(peak=peak, t_peak=t_peak, fwhm=fwhm, trough=trough, t_trough=t_trough)


def plot_subject(subject: str, regions: Dict[str, np.ndarray], params: Dict[str, np.ndarray], out_dir: str) -> None:
    plt.figure(figsize=(8, 6))
    for idx, key in enumerate(["V1", "V2", "V3", "V4", "All"]):
        color = GEM_COLORS[idx]
        y = regions[key]
        p = params[key]
        y_fine = hrf_model(p, TIME_FINE)
        
        # R² 계산
        y_fitted = hrf_model(p, TIME_POINTS)
        ss_res = np.sum((y - y_fitted) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        
        plt.plot(TIME_POINTS, y, "o", linewidth=2, color=color)
        plt.plot(TIME_FINE, y_fine, "-", linewidth=2, color=color, label=f"{key} (R²={r2:.3f})")

    plt.xlabel("Time (s)")
    plt.ylabel("Signal")
    plt.title(subject)
    plt.legend(loc="best")
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{subject}_hrf_fit.png")
    plt.savefig(out_path, dpi=160)
    plt.close()


def subject_from_path(path: str) -> str:
    base = os.path.basename(path)
    return base.replace("_hrf.mat", "")


def collect_subjects(input_dir: str, subjects: List[str]) -> List[str]:
    if subjects:
        return subjects

    mat_files = glob.glob(os.path.join(input_dir, "SUB*_hrf.mat"))
    return sorted(subject_from_path(path) for path in mat_files)


def save_params(subject: str, params: Dict[str, np.ndarray], out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    
    # MAT 파일 저장
    mat_path = os.path.join(out_dir, f"{subject}_params.mat")
    savemat(mat_path, {"params": params})
    
    # CSV 파일 저장
    csv_path = os.path.join(out_dir, f"{subject}_params.csv")
    rows = []
    param_names = ["p1", "q1", "p2", "q2", "a1", "a2", "c"]
    for roi, p in params.items():
        row = {"subject": subject, "roi": roi}
        row.update({param_names[i]: p[i] for i in range(7)})
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)


def process_subject(mat_path: str, verbose: bool) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    regions = extract_regions(mat_path)
    fitted = {name: fit_hrf(y, verbose=verbose) for name, y in regions.items()}
    return regions, fitted


def write_metrics_tsv(rows: List[List[str]], out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("subject\troi\tpeak\tt_peak\tfwhm\ttrough\tt_trough\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit HRF curves and extract metrics.")
    parser.add_argument("--input-dir", default=".", help="Directory with SUB*_hrf.mat files.")
    parser.add_argument("--output-dir", default="./output", help="Directory for params, plots, and stats.")
    parser.add_argument("--subjects", nargs="*", default=[], help="Optional subject IDs like SUB701.")
    parser.add_argument("--verbose", action="store_true", help="Print optimizer progress.")
    args = parser.parse_args()

    subjects = collect_subjects(args.input_dir, args.subjects)
    if not subjects:
        raise SystemExit("No subjects found. Provide --subjects or place SUB*_hrf.mat in input-dir.")

    stats_rows: List[List[str]] = []
    for subject in tqdm(subjects, desc="Processing subjects"):
        mat_path = os.path.join(args.input_dir, f"{subject}_hrf.mat")
        regions, fitted = process_subject(mat_path, verbose=args.verbose)
        save_params(subject, fitted, os.path.join(args.output_dir, "params"))

        for roi_name, params in fitted.items():
            y_fine = hrf_model(params, TIME_FINE)
            metrics = compute_metrics(y_fine)
            stats_rows.append(
                [
                    subject,
                    roi_name,
                    f"{metrics.peak:.6f}",
                    f"{metrics.t_peak:.6f}",
                    f"{metrics.fwhm:.6f}",
                    f"{metrics.trough:.6f}",
                    f"{metrics.t_trough:.6f}",
                ]
            )

        plot_subject(subject, regions, fitted, os.path.join(args.output_dir, "plots"))

    write_metrics_tsv(stats_rows, os.path.join(args.output_dir, "stats.tsv"))


if __name__ == "__main__":
    main()
