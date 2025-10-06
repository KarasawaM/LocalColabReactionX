# ---- general imports -----
from pathlib import Path
import json
import pprint
import sys

# ----- third-party imports -----
from loguru import logger
from ase.constraints import FixInternals


def setup_logger(verbose: bool = False, quiet: bool = False):
    """
    Configure global loguru logger.
    Call this at the beginning of any entry point (e.g., CLI, script).
    """    
    logger.remove()

    if quiet:
        level = "WARNING"
    elif verbose:
        level = "DEBUG"
    else:
        level = "INFO"

    # Set log format
    if verbose:
        log_format = ("<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
                      "<level>{level: <8}</level> | "
                      "<cyan>{module}:{function}:{line}</cyan> - "
                      "<level>{message}</level>"
                      )
    else:
        log_format = ("<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
                      "<level>{level: <8}</level> | "
                      "<cyan>{module}:{line}</cyan> - "
                      "<level>{message}</level>"
                      )

    logger.add(sys.stderr, level=level, format=log_format)
    logger.add("lcrx_{time:YYYY-MM-DD}.log", level=level, format=log_format, rotation="10 MB")


def log_launch_info():
    """
    Log basic information.
    """
    logger.info("LocalColabReactionX: Starting calculation...")
    logger.info(f"Command line arguments: {sys.argv}")


def log_toml_data(toml_data: dict, label="Input TOML"):
    logger.info(f"{label}:\n{pprint.pformat(toml_data)}")


def log_atoms_info(atoms, label="Atoms"):
    info_lines = [
        f"{label} structure loaded as follows:",
        f"Formula: {atoms.get_chemical_formula()}",
        f"Atom count: {len(atoms)}",
        f"Charge: {atoms.info.get('charge')}, Spin (mult): {atoms.info.get('spin')}",
        f"Calculator: {atoms.calc.__class__.__name__ if atoms.calc else 'None'}",
    ]

    logger.info("\n".join(info_lines))


def log_constraints(atoms, label="Atoms"):
    if not atoms.constraints:
        logger.info(f"No constraints applied on the {label}.")
        return

    constr_log = [f"Constraints on the {label} applied as follows (0-based index):"]
    for constr in atoms.constraints:
        if isinstance(constr, FixInternals):
            constr.initialize(atoms)

            for b in constr.bonds:
                indices = b[1]
                length = b[0]
                constr_log.append(f"Distance: indices={indices}, value={length:.2f} Å")

            for a in constr.angles:
                indices = a[1]
                angle = a[0]
                constr_log.append(f"Angle: indices={indices}, value={angle:.1f} deg")

            for d in constr.dihedrals:
                indices = d[1]
                dihedral = d[0]
                constr_log.append(f"Dihedral: indices={indices}, value={dihedral:.1f} deg")

        else:
            constr_log.append(f"{type(constr).__name__}: {constr}")

    if constr_log:
        logger.info("\n".join(constr_log))


