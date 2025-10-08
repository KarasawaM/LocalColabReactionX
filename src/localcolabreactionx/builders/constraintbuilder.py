from __future__ import annotations

from ase import Atoms
from ase.constraints import (ExternalForce, FixAtoms, FixBondLength,
                             FixInternals, Hookean)
from ase.units import kcal, mol

"""constraintbuilder.py
This module defines the ConstraintBuilder class, which parses constraint
settings from a TOML-loaded dictionary and applies them to an ASE Atoms object.

Supported constraint types:
- FixAtoms ("freeze"): fix atomic positions based on atom indices.
- FixInternals ("internals"): fix internal coordinates (distance, angle, dihedral).
- Hookean ("hookean"): apply harmonic restraints between atom pairs.

TOML format example:
[constraints]
[[constraints.freeze]]
atomid = "1, 2, 3"  # or [1, 2, 3]. 1-based index.

[[constraints.distance]]
atomid = "4, 5"  # or [4, 5]. 1-based index.
value = 1.5  # angstrom

[[constraints.angle]]
atomid = "6, 7, 8"  # or [6, 7, 8].  1-based index.
value = 90  # degrees

[[constraints.dihedral]]
atomid = "9, 10, 11, 12" # or [9, 10, 11, 12]. 1-based index.
value = 120  # degrees

[[constraints.hookean]]
atomid = "13, 14"  # or [13, 14]. 1-based index.
k = 1.0   # force constant (eV/Å²)
r0 = 1.5  # angstrom. When distance > r0, apply harmonic restraint.

Usage:
    builder = ConstraintBuilder(constraint_dict)
    builder.apply(atoms)
"""


class ConstraintBuilder:
    def __init__(self, constraint_data: dict):
        self.constraint_data = constraint_data

    def _parse_atomid(self, atomid: str | list[int]) -> list[int]:
        if isinstance(atomid, str):
            return [int(x.strip()) - 1 for x in atomid.split(",")]
        elif isinstance(atomid, list):
            return [int(x) - 1 for x in atomid]
        raise ValueError("atomid must be string or list of integers")

    def _apply_fix_atoms(self) -> list:
        constraints = []
        # ("freeze", []: default)
        for block in self.constraint_data.get("freeze", []):
            indices = self._parse_atomid(block["atomid"])
            constraints.append(FixAtoms(indices=indices))
        return constraints

    def _apply_fix_distances(self) -> list:
        constraints = []
        for block in self.constraint_data.get("distance", []):
            value = float(block["value"])
            atoms = self._parse_atomid(block["atomid"])
            fixbond = [value, atoms]
            constraints.append(FixInternals(bonds=[fixbond]))
        return constraints

    def _apply_fix_angles(self) -> list:
        constraints = []
        for block in self.constraint_data.get("angle", []):
            value = float(block["value"])
            atoms = self._parse_atomid(block["atomid"])
            fixangle = [value, atoms]
            constraints.append(FixInternals(angles_deg=[fixangle]))
        return constraints

    def _apply_fix_dihedrals(self) -> list:
        constraints = []
        for block in self.constraint_data.get("dihedral", []):
            value = float(block["value"])
            atoms = self._parse_atomid(block["atomid"])
            fixdihedral = [value, atoms]
            constraints.append(FixInternals(dihedrals_deg=[fixdihedral]))
        return constraints

    def _apply_hookean(self) -> list:
        constraints = []
        for block in self.constraint_data.get("hookean", []):
            atomid = self._parse_atomid(block["atomid"])
            k_kcalmol = float(block["k"])  # kcal mol-1 angstroam-2
            k = 2.0 * k_kcalmol * (kcal / mol)   # convert to eV angstroam-2. Hookean uses 1/2 * k(r-r0)^2.
            r0 = float(block["r0"])
            constraints.append(Hookean(atomid[0], atomid[1], k, r0))
        return constraints

    def apply(self, atoms: Atoms):
        all_constraints = []
        all_constraints += self._apply_fix_atoms()
        all_constraints += self._apply_fix_distances()
        all_constraints += self._apply_fix_angles()
        all_constraints += self._apply_fix_dihedrals()
        all_constraints += self._apply_hookean()
        if all_constraints:
            atoms.set_constraint(all_constraints)




