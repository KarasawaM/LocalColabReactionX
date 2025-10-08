from __future__ import annotations

import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import plotly.graph_objects as go
import plotly.io as pio
from loguru import logger
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt


def _pick_col(df, cands):
    df.columns = [c.strip() for c in df.columns]
    for c in cands:
        if c in df.columns:
            return c
    return None


def plot_1d_scan(df, x_label, y_label, out_prefix="scan_final"):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df[x_label], df[y_label], marker='o', linestyle='-')
    ax.set_xlabel("Index", fontsize=13)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel(r"$\Delta E$ [kcal/mol]", fontsize=13)
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.pdf", bbox_inches='tight')

    logger.info(f"1D scan plot saved to: {out_prefix}.pdf")


def plot_2d_scan(df, x_label, y_label, out_prefix="scan_final"):
    # Auto-pick columns
    x_col = _pick_col(df, [x_label])
    y_col = _pick_col(df, [y_label])
    z_col = _pick_col(df, ["Delta E [kcal/mol]"])
    if not x_col or not y_col or not z_col:
        raise ValueError(f"Could not determine columns. Found: {list(df.columns)}")

    # Ensure numeric
    for c in [x_col, y_col, z_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=[x_col, y_col, z_col]).copy()

    x = df[x_col].to_numpy()
    y = df[y_col].to_numpy()
    z = df[z_col].to_numpy()

    # === 2) Grid & interpolation for surface ===
    nx, ny = 100, 100
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    XI, YI = np.meshgrid(xi, yi)
    ZI = griddata((x, y), z, (XI, YI), method="linear")
    if np.all(np.isnan(ZI)):
        ZI = griddata((x, y), z, (XI, YI), method="nearest")

    # === 3) Build Plotly figure ===
    fig = go.Figure()

    # Surface trace (interpolated)
    fig.add_trace(go.Surface(
        x=XI, y=YI, z=ZI,
        colorbar=dict(title=z_col),
        name="Interpolated surface",
        opacity=0.8,
        # hovertemplate for surface
        hovertemplate=(
            f"{x_col}=%{{x:.4f}}<br>"
            f"{y_col}=%{{y:.4f}}<br>"
            f"{z_col}=%{{z:.4f}}"
            "<extra></extra>"
        ),
    ))

    # === Points (original data) with index + ΔE in hover ===
    indices = df.index.to_numpy()  # 0-based
    indices = indices + 1   # Convert to 1-based for display
    customdata = np.stack((indices, z), axis=-1)

    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode="markers",
        marker=dict(size=2, color="red"),
        name="Data points",
        customdata=customdata,
        # show index + ΔE（= z_col）, x_col, y_col in hover
        hovertemplate=(
            "index=%{customdata[0]}<br>"
            f"{z_col}=%{{customdata[1]:.4f}}<br>"
            f"{x_col}=%{{x:.4f}}<br>"
            f"{y_col}=%{{y:.4f}}"
            "<extra></extra>"
        ),
    ))

    fig.update_layout(
        title=f"3D Plot: {x_col} vs {y_col} vs {z_col}",
        scene=dict(
            xaxis_title=x_col,
            yaxis_title=y_col,
            zaxis_title=z_col,
        ),
        margin=dict(l=0, r=0, t=40, b=0),
        width=900,
        height=600,
    )

    # save to HTML
    out_html = "scan_final.html"
    pio.write_html(fig, file=out_html, include_plotlyjs="cdn", auto_open=False)
    logger.info(f"Interactive 2D-scan plot saved to: {out_html}")

