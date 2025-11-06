from __future__ import annotations

# --- general imports ---
import os
import subprocess
from abc import ABC, abstractmethod
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING

# --- third party imports ---
import numpy as np
from ase.vibrations import Vibrations
from ase.io import Trajectory, write
from loguru import logger

# --- LocalColabReactionX modules ---
from ..formats.gv_freq_format import write_gaussian_freq_log
if TYPE_CHECKING:
    from ..builders.toml2ase import SingleASEBuilder, DoubleASEBuilder


@contextmanager
def _temp_chdir(path: Path):
    """Temporarily change current working directory."""
    old = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(old)


# === Vibration Class ===

class BaseVibrationRunner(ABC):
    """
    Base class for vibration backends. Concrete classes implement
    how vibrational analysis is performed for each backend.
    """

    def __init__(
        self,
        builder: "SingleASEBuilder | DoubleASEBuilder",
        outdir: str | Path = "vibrations",
        out_prefix: str = "",
        suffix: str = ".xyz",
    ):
        self.builder = builder
        self.outdir = Path(outdir)
        self.out_prefix = out_prefix
        self.suffix = suffix
        self.outdir.mkdir(exist_ok=True)

    @abstractmethod
    def run(self, atoms, idx=None) -> None:
        """Run vibrational analysis for the given geometry."""
        ...

    # -------------------------
    # helpers (shared methods)
    # -------------------------
    def _write_geometry(self, atoms, outname: str) -> Path:
        """
        Write geometry for the current index into outdir.
        Returns: path to the written file (without extension).
        """
        atoms.info = {str(k): v for k, v in atoms.info.items()}
        path_noext = self.outdir / outname
        outfile = path_noext.with_suffix(self.suffix)

        if self.suffix == ".xyz":
            write(str(outfile), atoms, plain=True)
        else:
            write(str(outfile), atoms)

        return path_noext

    @staticmethod
    def _get_charge_mult(atoms) -> tuple[int, int]:
        """Fetch charge and multiplicity from atoms.info with safe defaults."""
        charge = int(atoms.info.get("charge", 0))
        mult = int(atoms.info.get("spin", 1))  # multiplicity (2S+1)
        return charge, mult


class ASEVibrationRunner(BaseVibrationRunner):
    """
    ASE numerical vibrations for ML/semi-empirical/backends (UMA, TBLite, EMT, ...).
    Also writes a GaussView-like frequency log using your existing utility.
    """

    def run(self, atoms, idx=None) -> None:
        if idx is not None:
            outname = f"{self.out_prefix}_idx{idx + 1}"
        else:
            outname = f"{self.out_prefix}"
        
        self._write_geometry(atoms, outname)

        logger.info(f"[ASE] Running Vibrations analysis for {outname} in {self.outdir}.\n"
                    f"Vibration table: {outname}_vib.txt\n"
                    f"Each mode: {outname}_vib/{outname}_mode.N.xyz, .traj\n"
                    f"Coordinates for sequential use: {outname}{self.suffix}, .gjf"
                    )

        with _temp_chdir(self.outdir):
            os.makedirs(f"{outname}_vib", exist_ok=True)  # ./maxima/idx{i+1}_vib
            os.chdir(f"{outname}_vib")

            # vibrational analysis
            vib = Vibrations(atoms, name=f"{outname}_mode")
            vib.run()
            vdata = vib.get_vibrations()
            vib.clean()  # remove temp files

            # Save table
            with open(f"{outname}_vib.txt", "w") as f:
                f.write(vdata.tabulate())

            # save Hessian
            H_2d = vdata.get_hessian_2d()
            logger.debug(f"Hessian shape: {H_2d.shape}, atoms: {len(atoms)}")
            np.save(f"{outname}_hessian.npy", H_2d)

            # Extract and save imaginary modes
            freqs = vdata.get_frequencies()
            imag_modes = np.where(freqs.imag > 0)[0]  # np.complex
            for m in imag_modes:
                vib.write_mode(m)
                trj = Trajectory(f"{outname}_mode.{m}.traj")
                write(f"{outname}_mode.{m}.xyz", trj[:], plain=True)
                logger.info(f"[ASE] Imaginary mode {m}: {freqs[m].imag:.4f}i cm^-1")

            os.chdir("..")  # back to outdir
            # Write GaussView-like log at outdir level
            charge, mult = self._get_charge_mult(atoms)
            write_gaussian_freq_log(
                vib=vib,
                atoms=atoms,
                filename=str(f"{outname}_gv.log"),
                charge=charge,
                mult=mult,
            )

            # write gaussian input file (for sebuqential use)
            write(
                f"{outname}.gjf", 
                images=atoms,
                charge=charge, mult=mult,
                extra="opt=(calcfc,tight,ts,noeigentest) freq=(noraman) nosymm"
            )

        logger.info(f"[ASE] Vibrational analysis for {outname} terminated.")


class GaussianVibrationRunner(BaseVibrationRunner):
    """
    Run Gaussian as an external process and compute analytic frequencies.
    This avoids expensive numerical finite-difference in ASE for G/ORCA cases.
    """

    def run(self, atoms, idx=None) -> None:
        if idx is not None:
            outname = f"{self.out_prefix}_idx{idx + 1}"
        else:
            outname = f"{self.out_prefix}"

        path_noext = self._write_geometry(atoms, outname)

        # manual freq input
        vib_data = (self.builder.data.get("vibration") or {})
        if len(vib_data) != 0: 
            cmd: str = vib_data.get("gaussian_cmd", "g16")
            mem: str = vib_data.get("mem", "8GB")
            nproc: int = int(vib_data.get("nprocshared", 4))
            route: str = vib_data.get("route", "B3LYP/6-31G(d) freq=noraman nosymm")

        # automatic freq input
        cmd: str = "g16"
        method = self.builder.data.get("calculator", {}).get("method", "")
        basis = self.builder.data.get("calculator", {}).get("basis", "")
        extra = self.builder.data.get("calculator", {}).get("extra", "")
        scf = self.builder.data.get("calculator", {}).get("scf", [])
        mem = self.builder.data.get("calculator", {}).get("totalmem", [])
        nproc = self.builder.data.get("calculator", {}).get("totalcores", [])

        scf = f"scf=({','.join(scf)})" if len(scf) > 0 else ""
        route = f"{method} {basis} freq=noraman {scf} {extra}"

        gjf = path_noext.with_suffix(".gjf")
        chk = path_noext.with_suffix(".chk")
        log = path_noext.with_suffix(".log")
        charge, mult = self._get_charge_mult(atoms)

        # Write Gaussian input (minimal but robust)
        self._write_gaussian_input(
            atoms=atoms,
            gjf_path=gjf,
            chk_path=chk,
            mem=mem,
            nprocshared=nproc,
            route=route,
            title=f"Frequency calculation for {outname}",
            charge=charge,
            mult=mult,
        )

        # Launch Gaussian
        logger.info(f"[Gaussian] Running Frequency calculation: {cmd} {gjf.name} (cwd={self.outdir})")
        with _temp_chdir(self.outdir):
            subprocess.run(
                [cmd, gjf.name],
                check=True,
                shell=False
            )
 
        logger.info(f"[Gaussian] Frequency calculation terminated. Log file: {log.name}")

    # -------------------------
    # Gaussian helpers
    # -------------------------
    @staticmethod
    def _write_gaussian_input(
        atoms,
        gjf_path: Path,
        chk_path: Path,
        mem: str,
        nprocshared: int,
        route: str,
        title: str,
        charge: int,
        mult: int,
    ) -> None:
        """Create a standard Gaussian input file."""
        lines = []
        lines.append(f"%chk={chk_path.name}")
        lines.append(f"%mem={mem}GB")
        lines.append(f"%nprocshared={nprocshared}")
        lines.append(f"#P {route}")
        lines.append("")                # blank
        lines.append(title)
        lines.append("")                # blank
        lines.append(f"{charge} {mult}")
        for a in atoms:
            x, y, z = a.position
            lines.append(f"{a.symbol:2s}  {x: .10f}  {y: .10f}  {z: .10f}")
        lines.append("\n")                # trailing blank

        gjf_path.write_text("\n".join(lines))


class ORCAVibrationRunner(BaseVibrationRunner):
    """
    Run ORCA as an external process and compute analytic (or numerical) frequencies.
    The route line is controlled via [vibration] options in TOML.
    """

    def run(self, atoms, idx=None) -> None:
        if idx is not None:
            outname = f"{self.out_prefix}_idx{idx + 1}"
        else:
            outname = f"{self.out_prefix}"

        path_noext = self._write_geometry(atoms, outname)

        # manual freq input
        vib_data = (self.builder.data.get("vibration") or {})
        if len(vib_data) != 0: 
            cmd: str = vib_data.get("orca_cmd", "orca")

        # automatic freq input
        cmd = self.builder.data.get("calculator", {}).get("orcapath", None)
        orcasimpleinput_in = self.builder.data.get("calculator", {}).get("orcasimpleinput", "")
        orcablocks_in = self.builder.data.get("calculator", {}).get("orcablocks", "")
        totalcores = self.builder.data.get("calculator", {}).get("totalcores", 16)
        totalmem = self.builder.data.get("calculator", {}).get("totalmem", 60)  # in GB

        # Maxcore: in MB / nprocs
        def _make_orca_block(orcablocks, totalcores, totalmem):
            orcablocks += f"\n%pal nprocs {totalcores} end\n"
            maxcore = (totalmem * 10 // totalcores) * 100  # 60 GB, 16procs -> maxcore 3700 MB.
            orcablocks += f"\n%maxcore {maxcore}\n"
            return orcablocks
        
        def _mod_orca_simpleinput(orcasimpleinput):
            import re
            pattern = re.compile(r'grad', re.IGNORECASE)
            simplein_list = orcasimpleinput.split()
            simplein_list = [s for s in simplein_list if not pattern.search(s)]
            simplein_list.append("Freq")
            simplein_list.insert(0, "!")
            return " ".join(simplein_list)

        orcablocks = _make_orca_block(orcablocks_in, totalcores, totalmem) if (totalcores and totalmem) else orcablocks_in
        orcasimpleinput = _mod_orca_simpleinput(orcasimpleinput_in)

        inp = path_noext.with_suffix(".inp")
        out = path_noext.with_suffix(".out")
        charge, mult = self._get_charge_mult(atoms)

        self._write_orca_input(atoms, inp, orcasimpleinput, charge, mult, orcablocks)

        logger.info(f"[ORCA] Running Frequency calculation: {cmd} {inp.name} (cwd={self.outdir})")

        # Launch ORCA
        with _temp_chdir(self.outdir):
            with open(out.name, "w") as fh:
                # Write both stdout and stderr into the same file
                subprocess.run(
                    [cmd, inp.name],
                    stdout=fh,
                    stderr=subprocess.STDOUT,
                    check=True,
                    shell=False
                )

        logger.info(f"[ORCA] Frequency calculation terminated. Log file: {out.name}")

    # -------------------------
    # ORCA helpers
    # -------------------------
    @staticmethod
    def _write_orca_input(atoms, inp_path: Path, orcasimpleinput: str, charge: int, mult: int, orcablocks) -> None:
        """Create a minimal ORCA input file."""
        coords = []
        for a in atoms:
            x, y, z = a.position
            coords.append(f"{a.symbol:2s}  {x: .10f}  {y: .10f}  {z: .10f}")

        lines = [
            orcasimpleinput,
            orcablocks,
            "",
            f"* xyz {charge} {mult}",
            *coords,
            "*",
            "",
        ]
        inp_path.write_text("\n".join(lines))


# === main functions ===

def pick_vibration_runner(
    builder: "SingleASEBuilder | DoubleASEBuilder",
    outdir: str | Path = "vibrations",
    out_prefix: str = "",
    suffix: str = ".xyz",
):
    """
    Factory that selects the best backend based on calculator type and
    optional [vibration] settings in the TOML.
    """
    calc_type = str(builder.data.get("calculator", {}).get("type", "")).lower()
    vib_cfg = (builder.data.get("vibration") or {})
    backend = str(vib_cfg.get("backend", "auto")).lower()

    # manual selection
    if backend == "ase":
        return ASEVibrationRunner(builder, outdir, out_prefix, suffix)
    if backend == "gaussian":
        return GaussianVibrationRunner(builder, outdir, out_prefix, suffix)
    if backend == "orca":
        return ORCAVibrationRunner(builder, outdir, out_prefix, suffix)

    # auto selection:
    if calc_type in {"uma", "uma-s-1p1", "uma-m-1p1", "tblite", "xtb", "emt"}:
        return ASEVibrationRunner(builder, outdir, out_prefix, suffix)
    if calc_type == "gaussian":
        return GaussianVibrationRunner(builder, outdir, out_prefix, suffix)
    if calc_type == "orca":
        return ORCAVibrationRunner(builder, outdir, out_prefix, suffix)

    # default fallback
    return ASEVibrationRunner(builder, outdir, out_prefix, suffix)


def run_vibrations_for_traj(
    traj,
    traj_indices: np.ndarray,
    builder: SingleASEBuilder | DoubleASEBuilder,
    outdir: str | Path = "vibrations",
    out_prefix: str = "",
    suffix: str = ".xyz",
) -> None:
    """
    Run vibrational analysis for each peak index using a pluggable backend:
      - UMA/TBLite/EMT -> ASE numerical vibrations
      - Gaussian/ORCA  -> external executable for analytic vibrations

    Notes:
      * For ASE backend we attach the calculator from builder.
      * For Gaussian/ORCA backend we just write inputs and call subprocess.
    """
    runner = pick_vibration_runner(builder, outdir=outdir, out_prefix=out_prefix, suffix=suffix)

    for i in traj_indices:
        atoms = traj[i]
        atoms.info = {str(k): v for k, v in atoms.info.items()}
        atoms.info["index"] = int(i + 1)

        # Attach calculator only for ASE backend
        if isinstance(runner, ASEVibrationRunner):
            atoms.calc = builder.make_calculator()
            logger.info("Calculator attached for ASE vibrations.")

        logger.info(f"Vibrational analysis at index {i+1} (1-based) via {runner.__class__.__name__}")
        runner.run(atoms=atoms, idx=i)


