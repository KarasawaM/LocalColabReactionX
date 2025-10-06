from __future__ import annotations

# ---- general modules ----
from pathlib import Path
import sys

# ---- third party modules ----
from loguru import logger
from ase.io import write
from ase.io.trajectory import Trajectory
# from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG

# ---- LocalColabReactionX modules ----
from ..builders.toml2ase import select_builder
from ..builders.optimizerselect import get_optimizer
from ..analysis.traj2csv import traj_to_csv
from ..analysis.optlog_parser import check_convergence_from_log
from ..formats.outpath_format import make_outpath, do_ase_write


def run_opt(
    toml_path: str | Path,
    history_prefix: str = "opt_history",
    logfile: str = "opt.log"
) -> None:
    """
    Run geometry optimization using parameters from TOML input.
    ASE's optimizer handles its own logging via logfile or trajectory.
    """
    builder = select_builder(toml_path)
    atoms = builder.build()

    calcdata = builder.data.get("calculation", {})
    fmax = float(calcdata.get("fmax", 0.05))
    maxstep = int(calcdata.get("maxstep", 1000))
    method = calcdata.get("method", "LBFGS")

    optimizer_class = get_optimizer(method)
    optimizer = optimizer_class(atoms, logfile=logfile, trajectory=f"{history_prefix}.traj")
    optimizer.run(fmax=fmax, steps=maxstep)

    # --- Write final structure ---
    # output format setting
    base, suffix = make_outpath(builder, multi_model=False)  # get basename and output suffix (e.g., .xyz, .gjf)
    final_path = f"{base}_opt{suffix}"
    do_ase_write(final_path, atoms, builder)

    # write history
    traj = Trajectory(f"{history_prefix}.traj")
    write(f"{history_prefix}.xyz", traj)    # write XYZ format for history
    traj_to_csv(traj, out_csv_path=f"{history_prefix}.csv")

    # check convergence
    _, _ = check_convergence_from_log(
        logfile_path=logfile, fmax_thresh=fmax, maxstep=maxstep,
        header="Optimization result:",
    )

    logger.info(f"Optimization terminated.")
    logger.info(f"Final structure written to: {final_path}")
    logger.info(f"Optimization process written to: {history_prefix}.traj, .xyz, .csv")


if __name__ == "__main__":
    # --- Debugging / Example usage ---
    # run with a TOML file path:  python3 opt_runner.py input.toml
    if len(sys.argv) != 2:
        print("Usage: python opt_runner.py input.toml")
        sys.exit(1)

    run_opt(sys.argv[1])
