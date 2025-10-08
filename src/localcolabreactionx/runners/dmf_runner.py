from __future__ import annotations

# --- general imports ---
import os
import threading
import sys
from pathlib import Path
from typing import TYPE_CHECKING

# ---- third party modules ----
import numpy as np
import ase.parallel
from ase.io import write
from loguru import logger
from types import MethodType
from dmf import DirectMaxFlux, interpolate_fbenm
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---- LocalColabReactionX modules ----
from ..builders.toml2ase import select_builder
from ..analysis.traj2csv import traj_to_csv
from ..analysis.peakdetector import detect_peaks_from_df
from ..analysis.vibration_runner import run_vibrations_for_traj
from ..formats.gv_irc_format import write_gaussian_irc_log
from ..formats.outpath_format import make_outpath
if TYPE_CHECKING:
    from ..builders.toml2ase import SingleASEBuilder, DoubleASEBuilder


# --- helper functions ---
def parallel_settings(builder: SingleASEBuilder | DoubleASEBuilder, n_parallel: int = None) -> dict:
    """
    Returns the parallel settings for the specified program.
    output:
        outpaths: dict. Working dirs for gaussian and orca. idx+1
        overrides_template: dict. Overrides for make_calclator function.
    """
    # parallel settings
    parallel = builder.data["calculation"].get("parallel", False)

    if parallel is False:
        return {}

    # if parallel. cpu and memory settings from builder data.
    totalcores = int(builder.data["calculator"].get("totalcores", 4))
    totalmem = int(builder.data["calculator"].get("totalmem", 8))  # in GB
    programm = builder.data["calculator"].get("type", "").lower()
    dmf_nmove = int(builder.data.get("calculation", {}).get("dmf_nmove", 5))

    if n_parallel is None or n_parallel <= 0:
        n_parallel = dmf_nmove

    # nproc and memory per job
    nproc_per_job = totalcores // n_parallel if totalcores is not None else 2
    mem_per_job_GB = totalmem // n_parallel if totalmem is not None else 2  # in GB

    if programm == "gaussian" or programm == "g16":
        mem_per_job = mem_per_job_GB  # in GB for gaussian
        overrides_template = {"outpath": "g16_idx{idx}/g16", "mem": f"{mem_per_job}GB", "nprocshared": str(nproc_per_job)}
        return overrides_template

    elif programm == "orca":
        orcablocks_in = builder.data["calculator"].get("orcablocks", "")
        mem_per_job = (mem_per_job_GB * 10 // nproc_per_job) * 100  # in MB per job * cores.
        overrides_template = {"outpath": "orca_idx{idx}", "orcablocks": f"{orcablocks_in}\n%pal nprocs {nproc_per_job} end\n%maxcore {mem_per_job}\n"}
        return overrides_template

    else:
        return {}


def wrap_calculator_with_logging(calc, idx: int, world=None):
    """
    Wrap calc.calculate to log start/end with point idx (0-based), rank, tid, pid.
    This version is robust to different binding/call patterns in ASE.
    """
    world = world or ase.parallel.world
    orig_calculate = calc.calculate  # this is usually a bound method

    def logged_calculate(self, *args, **kwargs):
        # self is the calculator instance (bound via MethodType below)
        rank, size = world.rank, world.size
        tid = threading.get_ident()
        pid = os.getpid()
        # try to peek "properties" if passed positionally or by keyword (for readability)
        props = None
        if "properties" in kwargs:
            props = kwargs.get("properties")
        elif len(args) >= 2:
            props = args[1]
        tag = f"[index={idx+1:03d} | rank={rank}/{size} | tid={tid} | pid={pid}]"
        logger.debug(f"START {tag} props={list(props) if props is not None else 'unknown'}")
        try:
            # orig_calculate is a bound method; do not pass self here
            return orig_calculate(*args, **kwargs)
        finally:
            e = getattr(self, "results", {}).get("energy", None)
            logger.debug(f"END   {tag} E={e}")

    # Bind wrapper to this instance so ASE always calls it as a bound method.
    calc.calculate = MethodType(logged_calculate, calc)
    return calc


def write_energy_and_force_history(mxflx: DirectMaxFlux):
    with open("energy_history.txt", "w") as f:
        for step, energies in enumerate(mxflx.history.energies):
            f.write(f"# Iteration {step}\n")
            for i, energy in enumerate(energies):
                f.write(f"Image {i}: {energy:.8f} eV\n")
            f.write("\n")

    with open("force_history.txt", "w") as f:
        for step, forces in enumerate(mxflx.history.forces):
            f.write(f"# Iteration {step}\n")
            for i, force_array in enumerate(forces):
                f.write(f"Image {i}:\n")
                for vec in force_array:
                    f.write(f"{vec[0]:.6f} {vec[1]:.6f} {vec[2]:.6f}\n")
                f.write("\n")


def _check_parallel(parallel: bool, world) -> None:
    if parallel and world.size > 1:
        mode = "mpi"
    elif parallel and world.size == 1:
        mode = "threads"
    else:
        mode = "serial"
    logger.info(f"Parallel mode: {mode} (world.size={world.size}, parallel={parallel})")


# --- main DMF runner function ---
def generate_initial_path(react, prod, correlated, nmove) -> np.ndarray:
    ref_images = [react, prod]
    mxflx_fbenm = interpolate_fbenm(ref_images, correlated=correlated, nmove=nmove)
    write("DMF_init.traj", mxflx_fbenm.images)
    write("DMF_init.xyz", mxflx_fbenm.images)
    np.save("DMF_init_coefs.npy", mxflx_fbenm.coefs.copy())
    return mxflx_fbenm.coefs.copy()


def solve_dmf(react, prod, coefs, charge, mult, calc_factory, nmove, update_teval, tol, parallel, **overrides_template) -> DirectMaxFlux:
    ref_images = [react, prod]
    mxflx = DirectMaxFlux(ref_images, coefs=coefs, nmove=nmove,
                          update_teval=update_teval, parallel=parallel)
    mxflx.add_ipopt_options({"output_file": "DMF_ipopt.out"})

    # parallel mode detection
    world = ase.parallel.world
    _check_parallel(parallel, world)

    # attach calculators to images
    for i, image in enumerate(mxflx.images, start=0):
        image.info["charge"] = charge
        image.info["spin"]   = mult

        if parallel and len(overrides_template) > 0:
            overrides = {**overrides_template, "outpath": overrides_template["outpath"].format(idx=i+1)}
        else:
            overrides = {}

        base_calc = calc_factory(**overrides)
        image.calc = wrap_calculator_with_logging(base_calc, idx=i, world=world)

    try:
        mxflx.solve(tol=tol)
    except Exception as e:
        logger.warning(f"DMF.solve failed: {e}")
        write("DMF_last_before_error.traj", mxflx.images)
        write("DMF_last_before_error.xyz", mxflx.images)

    return mxflx


def _evaluate_one_image(i, img, charge, mult, calc_factory, world, **overrides_template):
    """Compute a single image's energy (and forces if必要) and return (idx, atoms)."""
    try:
        atoms = img.copy()
        atoms.info["charge"] = img.info.get("charge", charge)
        atoms.info["spin"] = img.info.get("spin", mult)

        if len(overrides_template) > 0:
            overrides = {**overrides_template, "outpath": overrides_template["outpath"].format(idx=i+1)}
        else:
            overrides = {}

        # attach calculator and wrap for logging
        base_calc = calc_factory(**overrides)
        atoms.calc = wrap_calculator_with_logging(base_calc, idx=i, world=world)

        # one-point evaluation (energy only; forcesが必要なら get_forces に切替／併用)
        _ = atoms.get_potential_energy()
        # _ = atoms.get_forces()

        return i, atoms, None
    except Exception as e:
        logger.warning(f"[final SP] idx={i+1} failed: {e}")
        return i, img, e  # fall back to original image


def evaluate_final_images(images, charge, mult, calc_factory, parallel: bool = False, **overrides_template):
    """
    Re-evaluate final images (one-point) possibly in parallel, mirroring DMF's strategy:
      - MPI: world.size > 1 -> split indices across ranks, gather on root, broadcast
      - Threads: world.size == 1 and parallel -> thread pool
      - Serial: otherwise
    Returns: list[Atoms] with the same ordering as input.
    """
    world = ase.parallel.world
    n_points = len(images)
    _check_parallel(parallel, world)

    # --- MPI path ---
    if parallel and world.size > 1:
        logger.info("MPI parallelization not tested yet. Serial fallback will be used.")
        """
        my_idxs = [i for i in range(n_points) if (i % world.size) == world.rank]
        local_results = []
        for i in my_idxs:
            local_results.append(_evaluate_one_image(i, images[i], charge, mult, calc_factory, world, **overrides_template))

        # gather lists of tuples (i, atoms, err) to rank 0
        gathered = world.gather(local_results, destination=0)

        if world.rank == 0:
            # flatten and restore original order
            flat = [item for sub in gathered for item in sub]
            out = [None] * n_points
            for i, atoms, err in flat:
                out[i] = atoms
            # broadcast back to all ranks
            world.broadcast(out, destination=0)
            final_images = out
        else:
            # receive the reconstructed list
            final_images = world.broadcast(None, destination=0)

        return final_images
    """

    # --- threading path ---
    if parallel and world.size == 1:
        logger.debug("Using ThreadPoolExecutor for parallel SPC.")
        final_images = [None] * n_points
        max_workers = min(n_points - 2, 32)  # (npoints - 2) = dmf_nmove. max 32 threads.
        with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="final-spc") as ex:
            futs = {}
            for i in range(n_points):
                future = ex.submit(_evaluate_one_image, i, images[i], charge, mult, calc_factory, world, **overrides_template)
                futs[future] = i

            for fut in as_completed(futs):
                i, atoms, err = fut.result()
                final_images[i] = atoms

        return final_images

    # --- serial path ---
    final_images = []
    for i, img in enumerate(images):
        _, atoms, _ = _evaluate_one_image(i, img, charge, mult, calc_factory, world)
        final_images.append(atoms)

    return final_images


def run_dmf(toml_path: str | Path):
    # structure settings
    builder = select_builder(toml_path)
    react, prod = builder.build()    
    calcdata = builder.data.get("calculation", {})
    charge = react.info.get("charge", 0)
    mult = react.info.get("spin", 1)

    # dmf settings
    dmf_nmove = int(calcdata.get("dmf_nmove", 5))
    update_teval = bool(calcdata.get("update_teval", False))
    tol = calcdata.get("dmf_convergence", "tight")
    parallel = bool(calcdata.get("parallel", False))
    assert isinstance(update_teval, bool), "update_teval must be boolean"
    assert isinstance(parallel, bool), "parallel must be boolean"

    # parallel control. if not parallel, overrides_template = {}
    dmf_overrides_templt = parallel_settings(builder, n_parallel=dmf_nmove)

    # === Initial path ===
    # fbenm settings
    fbenm = calcdata.get("fbenm", "cfbenm")
    fbenm_nmove = int(calcdata.get("fbenm_nmove", 10))
    if fbenm.upper() == "FBENM":
        correlated = False
        logger.info("Using FBENM (uncorrelated) for initial path estimation.")
    else:
        correlated = True
        logger.info("Using CFBENM (correlated) for initial path estimation.")

    # Generate initial path
    logger.info(f"Generating initial path: fbenm_nmove = {fbenm_nmove}.")
    coefs = generate_initial_path(react=react, prod=prod,
                                  correlated=correlated, nmove=fbenm_nmove)
    logger.info("Initial path generation terminated normally.")

    # === Path optimization ===
    # Solve DMF
    logger.info(f"Optimizing the reaction path: dmf_nmove = {dmf_nmove}, update_teval = {update_teval}, tol = {tol}, parallel = {parallel}.")

    mxflx = solve_dmf(
        react=react, prod=prod, coefs=coefs,
        charge=charge, mult=mult,
        calc_factory=builder.make_calculator,
        nmove=dmf_nmove, update_teval=update_teval, tol=tol, 
        parallel=parallel, **dmf_overrides_templt
    )

    logger.info("Path optimization terminated normally.")

    # === Final Energy Evaluation ===
    logger.info("Final single-point energy calculation for all points.")
    final_images = evaluate_final_images(images=mxflx.images,
                                         charge=charge, mult=mult, 
                                         calc_factory=builder.make_calculator,
                                         parallel=parallel,
                                         **dmf_overrides_templt
                                         )

    # === Output ===
    # output final path
    base, suffix = make_outpath(builder)  # get basename and output suffix (e.g., .xyz, .gjf)
    write("DMF_final.traj", final_images)

    try:
        write(f"DMF_final{suffix}", final_images)
    except Exception as e:
        logger.warning(f"Failed to convert DMF_final.traj to {suffix}: {e}. Skipping...")

    e_df = traj_to_csv(final_images, out_csv_path="DMF_final.csv", return_df=True)

    # output tmax
    write("DMF_tmax.traj", mxflx.history.images_tmax)
    write(f"DMF_tmax{suffix}", mxflx.history.images_tmax)

    # output history
    write_energy_and_force_history(mxflx)
    logger.info("All DMF calculations terminated normally.")

    # output gview irc style log
    calclevel = builder.data["calculator"]
    write_gaussian_irc_log("DMF_final.traj", calclevel=calclevel,
                           output_log="DMF_final_gv.log")

    # === find peaks and vibrtional analysis ===
    max_peak_index, peak_indices = detect_peaks_from_df(e_df, prominence=0.01, distance=None, outdir="maxima")

    # vibrational analysis
    peak_vib = calcdata.get("peak_vibration", "highest")
    if peak_vib.lower() == "highest":
        traj_indices = max_peak_index
    elif peak_vib.lower() == "all":
        traj_indices = peak_indices
    else:
        logger.info("Vibration analysis for peaks skipped...") 
        traj_indices = []
 
    if len(traj_indices) > 0:
        run_vibrations_for_traj(traj=final_images, traj_indices=traj_indices, builder=builder, 
                                outdir="maxima", out_prefix="DMF_final", suffix=suffix)


if __name__ == "__main__":
    # --- Debugging / Example usage ---
    # run with a TOML file path:  python3 dmf_runner.py input.toml
    import sys
    if len(sys.argv) != 2:
        print("Usage: python dmf_runner.py input.toml")
        sys.exit(1)

    toml_path = sys.argv[1]

    logger.add("dmf_runner.log", level="INFO", rotation="1 MB", backtrace=True, diagnose=True)
    logger.info(f"Running DMF with input: {toml_path}")
    run_dmf(toml_path)
    logger.info("DMF calculation terminated normally.")

