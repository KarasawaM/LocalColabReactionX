from __future__ import annotations

# --- general modules ---
from pathlib import Path
import sys

# --- third-party modules ---
from ase.optimize.minimahopping import MinimaHopping, MHPlot
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.units import kcal, mol
from loguru import logger

# --- LocalColabReactionX modules ---
from ..builders.toml2ase import select_builder
from ..builders.optimizerselect import get_optimizer
from ..analysis.traj2csv import traj_to_csv
from ..analysis.optlog_parser import check_convergence_from_log
from ..formats.outpath_format import make_outpath


def run_minhop(toml_path: str | Path):
    builder = select_builder(toml_path)
    atoms = builder.build()

    # MD settings
    calcdata = builder.data.get("calculation", {})
    T0 = float(calcdata.get("T0", 1000.0))     # K, initial MD ‘temperature’
    beta1 = float(calcdata.get("beta1", 1.1))  # temperature adjustment parameter
    beta2 = float(calcdata.get("beta2", 1.1))  # temperature adjustment parameter
    beta3 = float(calcdata.get("beta3", 1 / 1.1))  # temperature adjustment parameter
    Ediff0_kcal = float(calcdata.get("Ediff0", 12))  # kcal/mol. initial energy acceptance threshold
    alpha1 = float(calcdata.get("alpha1", 0.98))  # temperature adjustment parameter
    alpha2 = float(calcdata.get("alpha2", 1 / 0.98))  # temperature adjustment parameter
    mdmin = int(calcdata.get("mdmin", 2))  # criteria to stop MD simulation (no. of minima)
    minima_threshold = float(calcdata.get("minima_threshold", 0.5))  # A, threshold for identical configs
    timestep = float(calcdata.get("timestep", 1.0))  # fs, MD time step
    totalsteps = int(calcdata.get("totalstep", 10))  # total steps for MD and optimization

    # optimizer settings
    method = calcdata.get("method", "FIRE")
    optimizer_class = get_optimizer(method)
    fmax = float(calcdata.get("fmax", 0.05))

    # output filenames
    out_prefix = "minima"
    minima_traj = f"{out_prefix}.traj"
    plot_out = "hop_summary.pdf"

    Ediff0 = Ediff0_kcal * (kcal / mol)  # convert to eV

    # run MH
    opt = MinimaHopping(
        atoms,
        Ediff0=Ediff0,
        T0=T0,
        beta1=beta1,
        beta2=beta2,
        beta3=beta3,
        alpha1=alpha1,
        alpha2=alpha2,
        optimizer=optimizer_class,
        fmax=fmax,
        timestep=timestep,
        minima_threshold=minima_threshold,
        mdmin=mdmin,
        minima_traj=minima_traj,
        logfile="hop.log"
    )
    opt(totalsteps=totalsteps)

    # plot result
    MHPlot().save_figure(plot_out)

    # convert .traj to coords and .csv
    base, suffix = make_outpath(builder, multi_model=True)  # get basename and output suffix (e.g., .xyz, .gjf)
    result_traj = Trajectory(minima_traj)
    try: 
        write(f"{out_prefix}{suffix}", result_traj)
    except Exception as e:
        logger.error(f"Failed to save final structures: {e}")

    traj_to_csv(result_traj, f"{out_prefix}.csv", return_df=False)

    # check convergence
    for i in range(totalsteps):
        opt_log = f"qn{i:05}.log"
        _, _ = check_convergence_from_log(
            logfile_path=opt_log, fmax_thresh=fmax, maxstep=1000000,
            header=f"Step [{i+1}/{totalsteps}]:"
        )

    logger.info(f"Minima-Hopping terminated: {totalsteps} steps.")
    logger.info(f"Results written to: {out_prefix}.traj, {suffix}, .csv")
    logger.info(f"Detailed simulation logs are available in hop.log.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python minhop_runner.py input.toml")
        sys.exit(1)
    run_minhop(sys.argv[1])
