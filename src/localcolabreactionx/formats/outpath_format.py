from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional, Any
from ase.io import write


def _resolve_input_path(bulder: Any) -> tuple[Path, Literal["single", "double"]]:
    """
    input: builder (SingleASEBuilder or DoubleASEBuilder)
    output: Path to the input (or reactant) structure file.
    """
    struct_data = bulder.struct_data
    # SingleASEBuilder
    if "input" in struct_data:
        return Path(struct_data["input"]), "single"  # single structure case

    # DoubleASEBuilder
    if "reactant" in struct_data and "product" in struct_data:
        reactant_path = Path(struct_data["reactant"])
        product_path = Path(struct_data["product"])
        if reactant_path.exists() and product_path.exists():
            return reactant_path, "double"  # double structure case
        else:
            raise FileNotFoundError("Reactant or product file not found.")
    raise KeyError("No valid input path found in struct_data.")


def make_outpath(builder: Any, multi_model=False) -> str:
    """
    input: builder (SingleASEBuilder or DoubleASEBuilder)
    multi_model: True if multiple models are required (e.g. MD, MinHop)
    """
    input_path, builder_type = _resolve_input_path(builder)
    base, suffix_in = input_path.stem, input_path.suffix

    # multi_model: only xyz or pdb supported.
    if multi_model:
        if suffix_in.lower() in [".pdb", ".sdf"]:
            suffix = ".pdb"
        else:
            suffix = ".xyz"
        return base, suffix

    # SingleASEBuilder case
    if builder_type == "single":
        if suffix_in.lower() in [".pdb", ".sdf"]:
            suffix = ".pdb"
        elif suffix_in.lower() in [".gjf", ".com"]:
            suffix = suffix_in
        else:
            suffix = ".xyz"

    # DoubleASEBuilder case
    elif builder_type == "double":
        if suffix_in.lower() in [".pdb", ".sdf"]:
            suffix = ".pdb"
        else:
            suffix = ".xyz"

    return base, suffix


def do_ase_write(final_path, atoms, builder):
    """
    Write final structure using ASE's write function.
    """
    suffix = Path(final_path).suffix

    if suffix.lower() in [".gjf", ".com"]:
        # write Gaussian input file
        write(final_path, atoms, 
              charge=builder.struct_data.get("charge", 0), 
              mult=builder.struct_data.get("mult", 1),
              )
    elif suffix.lower() in [".xyz"]:
        write(final_path, atoms, plain=True)
    else:
        # write other formats (e.g., .pdb)
        write(final_path, atoms)
