from __future__ import annotations

import sys
from pathlib import Path
import textwrap
from ase.io import read
from loguru import logger
from ase.units import Hartree
from ..formats.logo import get_logo


def _colabreaction_header(level: dict) -> str:
    logo = get_logo()

    level_txt = ", ".join(f"{k}={v}" for k, v in level.items())
    level_txt = f"# Scan compatible format, Energy({level_txt})"

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
"""


def modredandant_2d_txt(scan_blocks) -> str:
    """
    Example of scan blocks: scan_type, atom_id (0-based index), values (list of float)
        scan blocks: [
                    ('distance', [0, 5], array([3. , 2.5, 2. ])), 
                    ('angle', [1, 9, 0], array([110. , 107.5, 105. ]))
                    ]

    gaussian ModRedundant format for 2D scan:
        The following ModRedundant input section has been read:
        D       1       5       8      11 S   4 15.000                                
        A       8      11      12 S   4 5.0000 
    """
    # only for 2D scan
    if scan_blocks is None or len(scan_blocks) != 2:
        return ""

    if len(scan_blocks) == 2:
        modredundant_txt = ["The following ModRedundant input section has been read:"]
        for stype, atomid, values in (scan_blocks or []):
            if stype == "distance":
                line = f"B {atomid[0]+1:>7} {atomid[1]+1:>7} S {len(values)-1:>3} {values[1]-values[0]:.4f}"
            elif stype == "angle":
                line = f"A {atomid[0]+1:>7} {atomid[1]+1:>7} {atomid[2]+1:>7} S {len(values)-1:>3} {values[1]-values[0]:.4f}"
            elif stype == "dihedral":
                line = f"D {atomid[0]+1:>7} {atomid[1]+1:>7} {atomid[2]+1:>7} {atomid[3]+1:>7} S {len(values)-1:>3} {values[1]-values[0]:.4f}"
            else:
                raise ValueError(f"Unknown scan type: {stype}")

            modredundant_txt.append(line)
        modredundant_txt.append("\n")
        return "\n".join(modredundant_txt)


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


def write_gaussian_scan_log(traj_file, calclevel, output_log, scan_blocks=[]):
    traj = read(traj_file, ":")

    empty_indices = []
    last_energy = 0.0

    with open(output_log, "w") as f:
        f.write(_colabreaction_header(calclevel) + "\n")
        f.write(modredandant_2d_txt(scan_blocks))
        f.write("GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n")
        f.write("GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n\n")

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

            f.write("GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n")
            f.write(f"Step number   1 out of a maximum of   1 on scan point {i+1:>5} out of {len(traj):>5}\n")
            f.write("GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n\n")

        f.write("Normal termination of Gaussian. Formatted by ColabReaction.\n")

    # warning if empty energy
    if empty_indices:
        logger.warning("Energy not found at the following structure(s):")
        logger.warning("  " + ", ".join(f"#{i}" for i in empty_indices))

    logger.info(f"Log file compatible with GaussView is saved to: {output_log}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python traj_to_gaussian_log.py input.traj output.log scan_blocks")
        sys.exit(1)

    traj_file = sys.argv[1]
    output_log = sys.argv[2]
    scan_blocks = sys.argv[3] 

#    scan_blocks =  [
#        ('distance', [0, 5], [3. , 2.5, 2. ]), 
#        ('angle', [1, 9, 0], [110. , 107.5, 105. ])
#        ]
    write_gaussian_scan_log(traj_file=traj_file, output_log=output_log, calclevel={"level": "test"}, scan_blocks=scan_blocks)
