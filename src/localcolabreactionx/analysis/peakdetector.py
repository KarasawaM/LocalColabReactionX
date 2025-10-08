from __future__ import annotations

# --- general imports ---
from pathlib import Path

# --- third party imports ---
import pandas as pd
from scipy.signal import find_peaks
from loguru import logger


def _findpeaks_in_energy(df: pd.DataFrame, prominence: float = 0.01, distance: int = None):
    """
    output: peak_indices (0-based index) np.ndarray
    """
    x = df["index"].astype(int).to_numpy()
    y = df["Delta E [kcal/mol]"].astype(float).to_numpy()

    peaks, _ = find_peaks(y, prominence=prominence, distance=distance)  # add other params if needed.

    if len(peaks) == 0:
        logger.info("No peaks found in the data.")
        return None, None
    else:
        df_table = df.loc[peaks, ["index", "Delta E [kcal/mol]"]].reset_index(drop=True)
        peak_indices_1based = df_table["index"].astype(int).to_numpy()
        peak_indices = peak_indices_1based - 1 # convert to 0-based index
        
    logger.info(f"Detected peaks as follows:\n{df_table.to_string(index=False)}")

    return df_table, peak_indices


def detect_peaks_from_df(
    df: pd.DataFrame,
    prominence: float = 0.01,
    distance: int | None = None,
    outdir: str | Path = "maxima",
):
    """
    - Detect peaks (energy maxima) from a scan/DMF energy CSV.
    - Write `peaks.csv` under outdir.
    - Return (df_table: pd.DataFrame, 0-based indices for traj).
    """
    df_table, peak_indices = _findpeaks_in_energy(df, prominence=prominence, distance=distance)
    if df_table is None or peak_indices is None:
        return [], []

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    # Save a copy of peaks to CSV for auditability
    df_table.to_csv(outdir / "peaks.csv", index=False)

    # highest energy peak index (0-based)
    max_peak_index = df_table["Delta E [kcal/mol]"].idxmax()
    max_peak_index = peak_indices[max_peak_index]

    if peak_indices.size == 0:
        logger.warning("No energy maxima found. Skipping vibrational analysis.")
    else:
        logger.info(f"Found {peak_indices.size} peak(s). Saved list to: {outdir}/peaks.csv")

    return [max_peak_index], peak_indices


