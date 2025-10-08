from __future__ import annotations

import sys
import textwrap
from ase.io import read
from loguru import logger
from ase.units import Hartree
from ..formats.logo import get_logo


def _colabreaction_header(level: dict) -> str:
    logo = get_logo()

    level_txt = ", ".join(f"{k}={v}" for k, v in level.items())
    level_txt = f"# IRC compatible format, Energy({level_txt})"

    level_wrap_list = textwrap.wrap(
        level_txt,
        width=70,
        break_long_words=True,
        break_on_hyphens=False
    )

    level_wrap = "\n".join(level_wrap_list)

    return f"""{logo}
-----------------------------------------------------------------------
{level_wrap}
-----------------------------------------------------------------------

IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC
IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC"""


def _format_structure(atoms):
    lines = []
    lines.append("                         Standard orientation:                       ")
    lines.append("---------------------------------------------------------------------")
    lines.append("Center     Atomic      Atomic             Coordinates (Angstroms)")
    lines.append("Number     Number       Type             X           Y           Z")
    lines.append("---------------------------------------------------------------------")
    for i, atom in enumerate(atoms):
        atomic_num = atom.number
        x, y, z = atom.position
        lines.append(f"{i+1:5d}\t{atomic_num:<2d}\t0\t{x: .10f}\t{y: .10f}\t{z: .10f}")
    lines.append("---------------------------------------------------------------------")
    return "\n".join(lines)


def write_gaussian_irc_log(traj_file, calclevel, output_log, rxcoords: list[float] = None):
    traj = read(traj_file, ":")
    empty_indices = []
    last_energy = 0.0

    if rxcoords is None:
        rxcoords = [i for i in range(len(traj))]

    with open(output_log, "w") as f:
        f.write(_colabreaction_header(calclevel) + "\n\n")

        for i, atoms in enumerate(traj):
            coord = _format_structure(atoms)

            try:
                energy_ev = atoms.get_potential_energy()
                if energy_ev is None:
                    raise ValueError("Energy is None")
                energy_hartree = energy_ev / Hartree  # Convert eV to Hartree
                last_energy = energy_hartree
            except Exception:
                energy_hartree = last_energy
                empty_indices.append(i)

            f.write(coord + "\n")
            f.write(f"SCF Done:  E(scf) =  {energy_hartree: .10f}     A.U.\n\n")

            f.write("IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC\n")
            f.write(f"Pt {i} Step number   1 out of a maximum of  1\n")
            f.write(f"NET REACTION COORDINATE UP TO THIS POINT = {float(rxcoords[i]):20.10f}\n")
            f.write("IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC\n\n")

        f.write("Normal termination of Gaussian. Formatted by ColabReaction.\n")

    # warning if empty energy
    if empty_indices:
        logger.warning("Energy not found at the following structure(s):")
        logger.warning("  " + ", ".join(f"#{i}" for i in empty_indices))

    logger.info(f"Log file compatible with GaussView is saved to: {output_log}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python traj_to_gaussian_log.py input.traj output.log")
        sys.exit(1)

    traj_file = sys.argv[1]
    output_log = sys.argv[2]
    write_gaussian_irc_log(traj_file, output_log)
