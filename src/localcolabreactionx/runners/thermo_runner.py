from __future__ import annotations

# ---- general modules ----
from pathlib import Path
from typing import TYPE_CHECKING

# ---- third party modules ----
from loguru import logger

# ---- LocalColabReactionX modules ----
from ..builders.toml2ase import select_builder, SingleASEBuilder
from ..formats.outpath_format import make_outpath
from ..analysis.vibration_runner import pick_vibration_runner, ASEVibrationRunner
if TYPE_CHECKING:
    from typing import Union, Optional


def run_thermo(
    target: Union[str, Path, SingleASEBuilder],
    outdir: str | Path = "vibrations",
    out_prefix: Optional[str] = None,
) -> None:
    """
    Run vibrational analysis for either a TOML file or a SingleASEBuilder instance.

    Examples
    --------
    >>> run_vibration("input.toml")
    >>> builder = SingleASEBuilder("input.toml")
    >>> run_vibration(builder)
    """
    # --- Builder preparation ---
    if isinstance(target, (str, Path)):
        toml_path = Path(target)
        builder = select_builder(str(toml_path))
        if not isinstance(builder, SingleASEBuilder):
            raise TypeError(f"run_vibration supports only SingleASEBuilder, got {type(builder).__name__}")
    elif isinstance(target, SingleASEBuilder):
        builder = target
    else:
        raise TypeError("target must be a TOML path or a SingleASEBuilder")

    # --- Atoms and Runner setup ---
    atoms = builder.build()
    base, _ = make_outpath(builder, multi_model=False)
    out_prefix = base
    runner = pick_vibration_runner(builder, Path(outdir), out_prefix)
    logger.info(f"Output prefix set to: {out_prefix}")

    if isinstance(runner, ASEVibrationRunner):
        atoms.calc = builder.make_calculator()
        logger.info("Calculator attached to Atoms for ASE numerical vibrations.")

    # --- Execution ---
    logger.info(f"Running vibrational analysis ({runner.__class__.__name__}): prefix={out_prefix}")
    runner.run(atoms)
    logger.info("Vibrational analysis completed successfully.")

