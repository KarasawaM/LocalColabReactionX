# ---- general modules ----
from pathlib import Path
import sys
import os
from functools import partial

# ---- third party modules ----
import numpy as np
from loguru import logger
from ase.io import write
from ase import Atoms
from ase.vibrations import Vibrations
from ase.io.trajectory import Trajectory
from sella import Sella, Constraints, IRC

# ---- LocalColabReactionX modules ----
from ..builders.toml2ase import select_builder
from ..analysis.ircanalyzer import concat_irc_traj, make_df_oneside_irc, merge_df_irc
from ..analysis.traj2csv import traj_to_csv
from ..analysis.optlog_parser import check_convergence_from_log
from ..formats.outpath_format import make_outpath, do_ase_write
from ..formats.gv_irc_format import write_gaussian_irc_log
from ..utils.plot_utlis import make_2Dplot_from_df


def _check_calculator(builder):
    """
    Gaussian and ORCA have built-in TS optimization algorithms.
    This function raises an error if such calculators are detected.
    """
    calculator = builder.data.get("calculator", {})
    if calculator.get("name", "").lower() in ["gaussian", "orca"]:
        raise ValueError(f"Calculator {calculator.get('name')} has built-in TS optimization. Please consider using it directly.")


def _autograd_hessian(atoms, calc):
    atoms.get_potential_energy()
    atoms.get_forces()
    H = calc.get_hessian(atoms, vmap=False)
    H = 0.5 * (H + H.T)

    return H


def _vibrations_hessian(atoms):
    vib = Vibrations(atoms, name="vib_tmp")
    vib.run()
    vibdata = vib.get_vibrations()
    H = vibdata.get_hessian_2d()
    H = 0.5 * (H + H.T)
    vib.clean()
    return H
    

def _get_initial_hessian(H0_type="autograd", atoms=None, calc=None):
    """
    Initialize Hessian for Sella TS optimization.
    input: H0
        none, null, false -> None. Sella will estimate default hessian.
        calc -> try autograd Hessian, if fails, numerical Hessian by ASE Vibrations.
        filepath (.csv or .npy) -> load from file. hessian_2d() format.

    return: H0 as np.ndarray or None
        Cartesian Hessian matrix for Sella.
    """
    if H0_type.lower() in ["sella", "none", "null", "false"]:
        H0 = None

    elif Path(H0_type).is_file():
        # accept .csv or npy
        H0 = np.loadtxt(H0_type, delimiter=',') if H0_type.endswith('.csv') else np.load(H0_type)
        H0 = 0.5 * (H0 + H0.T)  # make sure the Hessian is symmetric

    elif H0_type.lower() in ["autograd"]:
        try:
            H0 = _autograd_hessian(atoms, calc)
        except Exception as e:
            logger.info(f"Failed to compute Hessian by autograd: {e}. Will try numerical Hessian.")
            H0 = _vibrations_hessian(atoms)

    elif H0_type.lower() in ["vib", "vibrations", "numerical", "freq", "frequency"]:
        H0 = _vibrations_hessian(atoms)

    return H0


def _set_hessian_function(update_hess, calc):
    """
    Set hessian_function for Sella.
    This function returns a callable that computes the Hessian matrix.
    
    input: update_hess
        None, false, no -> None. No hessian update during optimization.
        autograd -> use calc.get_hessian with vmap=False.
        vib -> use ASE Vibrations to compute numerical Hessian.
    """
    if update_hess.lower() in ["autograd", "analytical"]:
        hessian_function = partial(calc.get_hessian, vmap=False)

    elif update_hess.lower() in ["vib", "vibrations", "numerical", "freq", "frequency"]:
        hessian_function = _vibrations_hessian

    else:
        hessian_function = None

    return hessian_function


def _convert_hessian_to_coord_system(tsopt, H0_cart, internal: bool):
    """
    Convert Hessian to the desired coordinate system.
    If internal=True, convert cartesian Hessian to internal coordinates.
    If internal=False, return cartesian Hessian as is.

    Note: This function assumes H0 is in cartesian coordinates if internal=False.
    """
    if internal and H0_cart is not None:
        # Convert cartesian Hessian to internal coordinates
        # Using Sella's internal function
        H0 = tsopt.pes._convert_cartesian_hessian_to_internal(H0_cart)
    else:
        H0 = H0_cart
    return H0


def run_ts(
    toml_path: str | Path,
    history_prefix: str = "ts_opt",
) -> None:
    """
    transition state optimization using Sella.
    """
    # --- prepare system ---
    builder = select_builder(toml_path)
    _check_calculator(builder)
    atoms = builder.build()

    calcdata = builder.data.get("calculation", {})
    fmax = float(calcdata.get("fmax", 0.001))
    maxstep = int(calcdata.get("maxstep", 5000))
    internal = bool(calcdata.get("internal", True))
    H0_type = calcdata.get("initial_hess", "autograd")
    update_hess = str(calcdata.get("update_hess", "vib"))

    # irc settings
    irc_blocks = calcdata.get("irc", [])
    ircdata = irc_blocks[0] if len(irc_blocks) > 0 else None
    irc_flg = False
    endpt_opt = False

    if ircdata is not None:
        irc_flg = True
        dx = float(ircdata.get("irc_dx", 0.1))
        irc_fmax = float(ircdata.get("irc_fmax", 0.05))
        irc_steps = int(ircdata.get("irc_steps", 100))
        endpt_opt = ircdata.get("endpoint_opt", "LBFGS")

    # --- run Sella ---
    H0_cart = _get_initial_hessian(H0_type=H0_type, atoms=atoms, calc=atoms.calc)
    hessian_function = _set_hessian_function(update_hess=update_hess, calc=atoms.calc)

    os.makedirs("ts_opt_history", exist_ok=True)
    history_prefix = os.path.join("ts_opt_history", history_prefix)

    tsopt = Sella(
        atoms,
        order=1,
        internal=internal,
        trajectory=f'{history_prefix}.traj',
        logfile=f'{history_prefix}.log',
        hessian_function=hessian_function,
        H0=None,  # will set below: pes.set_H()
    )

    # in internal coordinates, need to convert the Hessian
    H0 = _convert_hessian_to_coord_system(tsopt, H0_cart, internal=internal)
    if H0 is not None:
        tsopt.pes.set_H(H0, initialized=True)

    tsopt.run(fmax, maxstep)

    # --- output results ---
    # write final TS structure
    base, suffix = make_outpath(builder)  # get basename and output suffix (e.g., .xyz, .gjf)
    final_path = f"{base}_TS{suffix}"
    do_ase_write(final_path, atoms, builder)

    # write history
    traj = Trajectory(f"{history_prefix}.traj")
    write(f"{history_prefix}.xyz", traj, plain=True)    # write XYZ format for history
    traj_to_csv(traj, out_csv_path=f"{history_prefix}.csv")

    # check convergence
    _, _ = check_convergence_from_log(
        logfile_path=f"{history_prefix}.log", fmax_thresh=fmax, maxstep=maxstep,
        header="TS optimization:")

    logger.info(f"Final TS structure written to: {final_path}")
    logger.info(f"TS Optimization process written to: {history_prefix}.traj, .xyz, .csv")
    logger.info("TS Optimization terminated.")

    # --- run IRC calculations ---
    if irc_flg:
        irc_endpoints = run_irc_both(atoms_ts=atoms, builder=builder, 
                                     update_hess=update_hess, dx=dx, 
                                     fmax=irc_fmax, steps=irc_steps)

        # write IRC endpoints
        do_ase_write(f"IRC-F/IRC-F_endpt{suffix}", irc_endpoints[0], builder)
        do_ase_write(f"IRC-R/IRC-R_endpt{suffix}", irc_endpoints[1], builder)
        logger.info(f"IRC endpoints written to: IRC-F/IRC-F_endpt{suffix}, IRC-R/IRC-R_endpt{suffix}")       
    else:
        logger.info("No IRC block found in input. Skipping IRC calculations.")

    if endpt_opt and irc_flg:
        run_endopt(irc_endpoints=irc_endpoints, 
                   builder=builder, method=endpt_opt, fmax=fmax, maxstep=maxstep)


def run_irc_both(atoms_ts: Atoms, builder, update_hess="vib", dx=0.1, fmax=0.05, steps=100) -> list[Atoms]:
    irc_endpoints = []
    # --- IRC forward/reverse---
    for direction, tag in (("forward", "F"), ("reverse", "R")):
        logger.info(f"Running {direction} IRC (IRC-{tag})...")
        os.makedirs(f"IRC-{tag}", exist_ok=True)
        atoms = atoms_ts.copy()
        calc = builder.make_calculator()
        atoms.calc = calc
        hess_func = _set_hessian_function(update_hess=update_hess, calc=atoms.calc)
        
        irc = IRC(
            atoms,
            trajectory=f"IRC-{tag}/IRC-{tag}.traj",
            logfile=f"IRC-{tag}/IRC-{tag}.log",
            hessian_function=hess_func,
            dx=dx,
        )

        irc.run(fmax=fmax, steps=steps, direction=direction)
        irc_endpoints.append(irc.atoms)

    # --- final outputs ---
    # merge IRC F/R and analyze
    os.makedirs("IRC_merged", exist_ok=True)
    traj_F = "IRC-F/IRC-F.traj"
    traj_R = "IRC-R/IRC-R.traj"

    concat_irc_traj(irc_f=traj_F, irc_r=traj_R, prefix="IRC")
    df_F = make_df_oneside_irc(traj_F)
    df_R = make_df_oneside_irc(traj_R)
    df_FtoR = merge_df_irc(df_F, df_R)
    df_RtoF = merge_df_irc(df_R, df_F)

    for df, direction in [(df_FtoR, "FtoR"), (df_RtoF, "RtoF")]:
        df.to_csv(f"IRC_merged/IRC_{direction}.csv", index=False)
        make_2Dplot_from_df(df=df,
                            df_xlabel="IRC [amu^1/2 Angstrom]",
                            df_ylabel="DeltaE [kcal/mol]",
                            plot_xlabel=r"IRC [amu$^{1/2}$ Å]",
                            plot_ylabel=r"$\Delta E$ vs. TS [kcal mol$^{-1}$]",
                            outpath=f"IRC_merged/IRC_{direction}.pdf")

        # gaussian irc log format
        rxcoords = df["IRC [amu^1/2 Bohr]"].to_numpy()
        write_gaussian_irc_log(f"IRC_merged/IRC_{direction}.traj", calclevel=builder.data["calculator"],
                               rxcoords=rxcoords, output_log=f"IRC_merged/IRC_{direction}_gv.log")

    logger.info("IRC calculations terminated.")
    
    return irc_endpoints


def run_endopt(irc_endpoints: list[Atoms], builder, method="LBFGS", fmax=0.01, maxstep=1000) -> Atoms:
    from ..builders.optimizerselect import get_optimizer
    optimizer_class = get_optimizer(method)

    for i, atoms in enumerate(irc_endpoints):
        direction = ["F", "R"][i]
        outdir = f"IRC-{direction}"
        history_dir = f"{outdir}/endpt_opt_history"
        base = f"IRC-{direction}_endpt_opt"
        os.makedirs(f"{history_dir}", exist_ok=True)
        logger.info(f"Running endpoint optimization for IRC-{direction}...")

        calc = builder.make_calculator()
        atoms.calc = calc

        optimizer = optimizer_class(atoms, logfile=f"{history_dir}/{base}.log",
                                    trajectory=f"{history_dir}/{base}.traj")
        optimizer.run(fmax=fmax, steps=maxstep)

        # write final structure
        _, suffix = make_outpath(builder, multi_model=False)  # get basename and output suffix (e.g., .xyz, .gjf)
        final_path = f"{outdir}/{base}{suffix}"
        do_ase_write(final_path, atoms, builder)

        # save history        
        traj = Trajectory(f"{history_dir}/{base}.traj")
        write(f"{history_dir}/{base}.xyz", traj)    # write XYZ format for history
        traj_to_csv(traj, out_csv_path=f"{history_dir}/{base}.csv")

        # check convergence
        _, _ = check_convergence_from_log(
            logfile_path=f"{history_dir}/{base}.log", fmax_thresh=fmax, maxstep=maxstep,
            header=f"[IRC-{direction}] Endpoint optimization:")

        logger.info(f"[IRC-{direction}] Optimized structure written to: {final_path}")
        logger.info(f"[IRC-{direction}] Endpoint optimization process written to: {history_dir}/{base}.traj, .xyz, .csv")
        logger.info(f"[IRC-{direction}] Endpoint optimization terminated.")


if __name__ == "__main__":
    # --- Debugging / Example usage ---
    # run with a TOML file path:  python3 opt_runner.py input.toml
    if len(sys.argv) != 2:
        print("Usage: python opt_runner.py input.toml")
        sys.exit(1)

    run_ts(sys.argv[1])
