from __future__ import annotations

# --- general imports ---
from pathlib import Path
import sys

# --- third-party modules ---
import ase.parallel
from ase.io import write
from ase.mep import NEB
from loguru import logger

# --- LocalColabReactionX modules ---
from ..builders.toml2ase import select_builder
from ..builders.optimizerselect import get_optimizer
from ..analysis.traj2csv import traj_to_csv
from ..runners.dmf_runner import parallel_settings, wrap_calculator_with_logging
from ..formats.gv_irc_format import write_gaussian_irc_log
from ..formats.outpath_format import make_outpath
from ..analysis.vibration_runner import run_vibrations_for_traj
from ..analysis.peakdetector import detect_peaks_from_df


def opt_neb_mep(images, charge, mult, calc_factory,
                idpp, k_spring, climb, parallel, neb_method,
                opt_method, fmax, maxstep, **overrides_template):
    
    # parallel settings
    world = ase.parallel.world
    mode = ("serial", "threads")[parallel and world.size == 1] if not (parallel and world.size > 1) else "mpi"
    logger.info(f"NEB parallel mode: {mode} (world.size={world.size}, parallel={parallel})")

    for i, image in enumerate(images):
        image.info["charge"] = charge
        image.info["spin"] = mult

        if parallel and len(overrides_template) > 0:
            overrides = {**overrides_template, "outpath": overrides_template["outpath"].format(idx=i+1)}
        else:
            overrides = {}

        base_calc = calc_factory(**overrides)
        image.calc = wrap_calculator_with_logging(base_calc, idx=i, world=world)

    # initial interpolation
    neb = NEB(images, climb=climb, k=k_spring,
              remove_rotation_and_translation=True, parallel=parallel, method=neb_method)

    neb.interpolate(method="idpp" if idpp else "linear")
    write("NEB_init.xyz", images)
    write("NEB_init.traj", images)

    # optimize mep
    optimizer_class = get_optimizer(opt_method)
    opt = optimizer_class(neb, trajectory="NEB.traj", logfile="NEB.log")
    opt.run(fmax=fmax, steps=maxstep)

    # aseneb: images[0] and images[-1] are not calculated. Get their energies from initial/final states.
    if neb_method == "aseneb":
        images[0].get_potential_energy()
        images[-1].get_potential_energy()

    return images


def run_neb(toml_path: str | Path):
    """
    Run a NEB calculation configured by a TOML file.

    Required TOML sections (same as DMF):
        [calculation]
        type = "NEB"
        neb_method = "aseneb"   # or "improvedtangent",  "eb", "spline", "string"
        n_images = 7            # total images including endpoints
        climb = true            # climbing image on/off
        idpp = true             # use IDPP interpolation
        method = "FIRE"         # optimizer ("FIRE", "BFGS", "LBFGS", "BFGSLineSearch")
        fmax = 0.05             # eV/Å
        maxstep = 500           # optimizer max steps
        k = 0.1                 # (optional) NEB spring constant (eV/Å^2)
        parallel = false        # (optional) run images in parallel. only for CPU calculators.
    """
    builder = select_builder(toml_path)
    react, prod = builder.build()
    calcdata = builder.data.get("calculation", {})
    charge = react.info.get("charge", 0)
    mult = react.info.get("spin", 1)

    # ---- Read NEB controls ----
    n_images: int = int(calcdata.get("n_images", 7))
    if n_images < 2:
        raise ValueError("n_images must be >= 2 (includes endpoints).")

    # neb settings
    neb_method: str = calcdata.get("neb_method", "aseneb").lower()
    climbing_image: bool = bool(calcdata.get("climb", calcdata.get("climbing_image", True)))
    use_idpp: bool = bool(calcdata.get("idpp", True))
    opt_method: str = calcdata.get("method", "FIRE")
    fmax: float = float(calcdata.get("fmax", 0.05))
    maxstep: int = int(calcdata.get("maxstep", 500))
    k_spring: float = float(calcdata.get("k", 0.1))  
    parallel = bool(calcdata.get("parallel", False))
    assert isinstance(use_idpp, bool), "use_idpp must be boolean"
    assert isinstance(parallel, bool), "parallel must be boolean"
    
    logger.info(
        f"NEB settings: n_images={n_images}, climb={climbing_image}, "
        f"idpp={use_idpp}, opt_method={opt_method}, fmax={fmax}, maxstep={maxstep}, k={k_spring}, "
        f"parallel={parallel}, neb_method={neb_method}"
    )

    # ---- Prepare images ----
    # set initial images
    r = react.copy()
    p = prod.copy()

    images = [r]
    images += [r.copy() for _ in range(n_images - 2)]
    images += [p]

    # parallel settings. same as DMF runner
    n_parallel = n_images - 2
    neb_overrides_templt = parallel_settings(builder, n_parallel=n_parallel)

    # opt neb image
    images = opt_neb_mep(
        images=images, charge=charge, mult=mult, calc_factory=builder.make_calculator,
        idpp=use_idpp, k_spring=k_spring, climb=climbing_image, parallel=parallel, opt_method=opt_method,
        fmax=fmax, maxstep=maxstep, neb_method=neb_method, **neb_overrides_templt
    )

    logger.info("NEB calculations terminated normally.")

    # == output ==
    write("NEB_final.traj", images)

    _, suffix = make_outpath(builder)  # get basename and output suffix (e.g., .xyz, .gjf)
    try:
        write(f"NEB_final{suffix}", images)
    except Exception as e:
        logger.warning(f"Failed to convert NEB_final.traj to {suffix}: {e}. Skipping...")

    e_df = traj_to_csv(images, out_csv_path="NEB_final.csv", return_df=True)

    # write Gaussian-style IRC log
    calclevel = builder.data["calculator"]
    write_gaussian_irc_log("NEB_final.traj", calclevel=calclevel,
                           output_log="NEB_final_gv.log")

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
        run_vibrations_for_traj(traj=images, traj_indices=traj_indices, builder=builder, 
                                outdir="maxima", out_prefix="NEB_final", suffix=suffix)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python neb_runner.py input.toml")
        sys.exit(1)

    # Light per-file logger (main CLI sets up a richer one)
    logger.add("neb_runner.log", level="INFO", rotation="1 MB", backtrace=True, diagnose=True)
    logger.info(f"Running NEB with input: {sys.argv[1]}")
    run_neb(sys.argv[1])
    logger.info("NEB calculation terminated normally.")
