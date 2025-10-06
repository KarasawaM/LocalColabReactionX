from __future__ import annotations

# ---- general modules ----
from pathlib import Path
import sys

# ---- third party modules ----
from loguru import logger
from ase import units
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md import MDLogger

# ---- LocalColabReactionX modules ----
from ..builders.toml2ase import select_builder
from ..analysis.traj2csv import traj_to_csv
from ..formats.outpath_format import make_outpath, do_ase_write


def run_nvemd(
    toml_path: str | Path,
    trajfile: str = "md.traj",
    logfile: str = "md.log"
) -> None:
    """
    Run geometry optimization using parameters from TOML input.
    ASE's optimizer handles its own logging via logfile or trajectory.
    """
    builder = select_builder(toml_path)
    atoms = builder.build()

    calcdata = builder.data.get("calculation", {})
    timestep = float(calcdata.get("timestep", 1.0))
    nsteps = int(calcdata.get("nsteps", 1000))
    gen_temp = float(calcdata.get("gen_temp", 300.0))
    nstlog = int(calcdata.get("nstlog", 10))
    logfile = calcdata.get("logfile", logfile)
    trajfile = calcdata.get("trajfile", trajfile)

    # logging function
    def print_dyn():
        nstep = dyn.get_number_of_steps() + 1
        etot  = atoms.get_total_energy() / (units.kcal / units.mol)
        ekin  = atoms.get_kinetic_energy() / (units.kcal / units.mol)
        epot  = atoms.get_potential_energy() / (units.kcal / units.mol)
        temp_K = atoms.get_temperature()
        logger.info(f"   {nstep:>8d}     {etot:.9f}     {ekin:.9f}    {epot:.9f}   {temp_K:.2f}")

    # --- velocity verlet dynamics ---
    MaxwellBoltzmannDistribution(atoms, temperature_K=gen_temp)
    dyn = VelocityVerlet(atoms=atoms,
                         timestep=timestep * units.fs,
                         trajectory=trajfile,
                         loginterval=nstlog
                         )

    # Print statements
    dyn.attach(print_dyn, interval=nstlog)
    logger.info("nstep     Etot(kcal/mol)    Ekin(kcal/mol)    Epot(kcal/mol)    T(K)")

    dyn.attach(MDLogger(dyn, atoms, logfile, header=True, stress=False, peratom=False, mode="w"), interval=nstlog)
    dyn.run(nsteps)

    # final output
    base, suffix = make_outpath(builder, multi_model=True)  # get basename and output suffix (e.g., .xyz, .gjf)

    try:
        traj = Trajectory(f"{trajfile}")
        write(f"{base}_md{suffix}", traj)    # write final structure
        write(f"{base}_final{suffix}", traj[-1])
    except Exception as e:
        logger.error(f"Error writing trajectory: {e}")

    # final log
    logger.info("NVE MD simulation terminated.")
    logger.info(f"Trajectory written to: {trajfile}, {base}_md{suffix}")
    logger.info(f"Final snapshot written to: {base}_final{suffix}")


if __name__ == "__main__":
    # --- Debugging / Example usage ---
    # run with a TOML file path:  python3 nve_md_runner.py input.toml
    if len(sys.argv) != 2:
        print("Usage: python nve_md_runner.py input.toml")
        sys.exit(1)

    run_nvemd(sys.argv[1])
