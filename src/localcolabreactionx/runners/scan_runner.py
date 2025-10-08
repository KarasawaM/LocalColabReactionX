from __future__ import annotations

# --- general imports ---
import os
import tomllib
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
from typing import List, Tuple

# --- third-party modules ---
from ase import Atoms
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.constraints import FixInternals
from ase.units import Hartree, kcal, mol
from loguru import logger

# --- LocalColabReactionX modules ---
from ..builders.toml2ase import select_builder
from ..builders.optimizerselect import get_optimizer
from ..analysis.optlog_parser import check_convergence_from_log
from ..analysis.peakdetector import detect_peaks_from_df
from ..analysis.vibration_runner import run_vibrations_for_traj
from ..analysis.scan_visualizer import plot_1d_scan, plot_2d_scan
from ..formats.outpath_format import make_outpath, do_ase_write
from ..formats.gv_scan_format import write_gaussian_scan_log


class BaseScanRunner(ABC):
    """
    Abstract base class for multi-dimensional scan runners.

    Subclasses must implement how to generate initial structures 
    and how to manage dependencies between scan points.
    """

    def __init__(self, toml_path: str):
        """
        Parameters:
            toml_path (str): The path to the TOML file. 
        """
        self.builder = select_builder(toml_path)  # from toml2ase module
        self.atoms = self.builder.build()
        self.calcdata = self.builder.data.get("calculation", {})
        self.calc_factory = self.builder.make_calculator

        if self.builder.data.get("scan"):
            self.scandata = self.builder.data.get("scan", {})
        else:
            raise ValueError("No [scan] block found in TOML file.")

        # Optimization parameters
        method = self.calcdata.get("method", "FIRE")
        self.optimizer_class = self.get_optimizer(method)
        self.fmax = float(self.calcdata.get("fmax", 0.05))
        self.maxstep = int(self.calcdata.get("maxstep", 500))
        self.peak_vib = self.calcdata.get("peak_vibration", "highest")

        # Scan grid setup
        self.scan_blocks = self.generate_scan_block()
        self.scan_grids = [v for _, _, v in self.scan_blocks]
        self.all_points = self.generate_scan_points()

    @staticmethod
    def get_optimizer(method: str):
        return get_optimizer(method)

    def _parse_scan_block(self, block):
        atomid = [int(x) - 1 for x in block["atomid"].split(",")]
        start, stop, num = map(float, block["value"].split(","))
        values = np.linspace(start, stop, int(num))
        return atomid, values

    def generate_scan_block(self):
        """
        Parse scan blocks from TOML data.
        Returns a list of tuples (scan_type, [atomid], array(values)). 0-based index.
        """
        scan_blocks = []
        for stype, blocks in self.scandata.items():
            if not isinstance(blocks, list):
                continue  # Skip boolean or other non-list keys (e.g., "sequential")
            for block in blocks:
                atomid, values = self._parse_scan_block(block)
                scan_blocks.append((stype, atomid, values))
        return scan_blocks

    def generate_scan_points(self) -> list[Tuple[float]]:
        """
        Generate scan points with serpentine order for 2D scans.
        Scan order Example:
            scan.distance -> outer loop
            scan.angle -> inner loop
        """
        if len(self.scan_grids) > 3:
            raise ValueError("Relaxed scan not supported for N > 3 dimensions.")
        
        elif len(self.scan_grids) == 1:
            scan_points = [(v,) for v in self.scan_grids[0]]

        # serpentine order for 2D scan (zigzag)
        elif len(self.scan_grids) == 2:
            inner = self.scan_grids[0]  # e.g. distance inner loop 
            outer = self.scan_grids[1]  # e.g. angle outer loop 

            scan_points = []
            # inner: (m, inn), outer: (n, out)
            for n, out in enumerate(outer):
                if n % 2 == 0:
                    for inn in inner:
                        scan_points.append((inn, out))
                else:
                    for inn in inner[::-1]:
                        scan_points.append((inn, out))

        elif len(self.scan_grids) == 3:
            inner = self.scan_grids[0]  # e.g. distance
            middle = self.scan_grids[1]  # e.g. angle
            outer = self.scan_grids[2]  # e.g. dihedral

            scan_points = []
            for o in outer:
                for m in middle:
                    for i in inner:
                        scan_points.append((i, m, o))

        return scan_points

    def generate_constraints(self, point: Tuple[float]) -> List[FixInternals]:
        """
        Generate FixInternals constraints for a given scan point.

        Args:
            point (tuple): Tuple of scan values for THIS point. e.g. (1.0, 90)

        Returns:
            List of FixInternals constraints to apply.
        """
        constraints = []
        for (stype, atomid, _), val in zip(self.scan_blocks, point):
            match stype:
                case "distance":
                    constraints.append(FixInternals(bonds=[[val, atomid]]))
                case "angle":
                    constraints.append(FixInternals(angles_deg=[[val, atomid]]))
                case "dihedral":
                    constraints.append(FixInternals(dihedrals_deg=[[val, atomid]]))
                case _:
                    raise ValueError(f"Unsupported scan type: {stype}")
        return constraints

    @abstractmethod
    def generate_initial_structure(self, i: int, point: Tuple[float]) -> Atoms:
        """
        Abstract method to generate the initial structure for a given scan point.
        """
        pass

    def _make_df_from_results(self, traj_path: str, 
                              convergence_status: List[Tuple[float, bool]]) -> pd.DataFrame:
        """
        Create a pandas DataFrame from the trajectory and convergence status.
        Parameters:
            traj_path (str): Path to the trajectory file.
            convergence_status (list of tuple): List of convergence status for each point. fmax and bool.
        Returns:
            pd.DataFrame: DataFrame containing scan results.
        """
        # Prepare scan info
        actual_headers = []
        for s in self.scan_blocks:
            if s[0] == "distance":
                actual_headers.append(f"{s[0]} [angstrom]")
            elif s[0] == "angle" or s[0] == "dihedral":
                actual_headers.append(f"{s[0]} [degree]")

        # Compute actual internal values
        actual_values = []
        energies = []
        for atoms_i in Trajectory(traj_path):
            try:
                energy = atoms_i.get_potential_energy()
            except Exception as e:
                logger.warning(f"Failed to get energy: {e}")
                energy = None
            energies.append(energy)

            values = []
            for stype, atomid, _ in self.scan_blocks:
                match stype:
                    case "distance":
                        val = atoms_i.get_distance(*atomid)
                    case "angle":
                        val = atoms_i.get_angle(*atomid)
                    case "dihedral":
                        val = atoms_i.get_dihedral(*atomid)
                    case _:
                        val = None
                values.append(val)
            actual_values.append(values)

        # Prepare DataFrame for CSV output
        df_data = {
            "index": range(1, len(self.all_points) + 1),
            **dict(zip(actual_headers, zip(*actual_values))),  # transpose list of rows to columns -> dict
            "energy [eV]": energies,
            "energy [Hartree]": [e / Hartree if e is not None else None for e in energies],
            "energy [kcal/mol]": [e / (kcal / mol) if e is not None else None for e in energies],
            "Delta E [kcal/mol]": [
                (e - energies[0]) / (kcal / mol) if e is not None else None for e in energies
            ],
            "fmax [eV/angstrom]": [fmax_calc for fmax_calc, _ in convergence_status],
            "converged": [converged for _, converged in convergence_status],
        }

        df = pd.DataFrame(df_data)
        traj_prefix, _ = os.path.splitext(traj_path)
        df.to_csv(f"{traj_prefix}.csv", index=False)
        return df

    # --- Main run function ---
    def run(self,
            out_prefix: str = "scan_final",
            log_prefix: str = "scan",
            log_dir: str = "opt_logs",
            fmax: float = None,
            maxstep: int = None,
            optimizer_class=None,
            ):
        """
        Run optimization for all scan points.

        Parameters:
            out_prefix (str): Prefix for output files.
            log_dir (str): Directory to store per-step logs.
            fmax (float): Optional override of force threshold.
            maxstep (int): Optional override of max steps.
            optimizer_class: Optional override of optimizer class (e.g., FIRE).
        """
        os.makedirs(log_dir, exist_ok=True)
        traj = Trajectory(f"{out_prefix}.traj", mode="w", properties=["energy", "forces"])
        energies = []
        scan_values = []
        convergence_status = []

        fmax = fmax if fmax is not None else self.fmax
        maxstep = maxstep if maxstep is not None else self.maxstep
        optimizer_class = optimizer_class if optimizer_class else self.optimizer_class

        logger.debug(f"Scan data: {self.scandata}")
        logger.debug(f"Scan blocks: {self.scan_blocks}")
        logger.debug(f"Scan grids: {self.scan_grids}")
        logger.info(f"Scan points: {self.all_points}")

        for i, point in enumerate(self.all_points):
            atoms_i = self.generate_initial_structure(i, point)

            logfile = os.path.join(log_dir, f"{log_prefix}_{i+1:04}.log")
            optimizer = optimizer_class(atoms_i, logfile=logfile, trajectory=None)
            optimizer.run(fmax=fmax, steps=maxstep)

            try:
                energy = atoms_i.get_potential_energy()
            except Exception as e:
                logger.warning(f"Failed to get energy at step {i}: {e}")
                energy = None

            traj.write(atoms_i)
            energies.append(energy)
            scan_values.append(point)
            # convergence check
            if energy is not None:
                fmax_calc, converged = check_convergence_from_log(
                    logfile_path=logfile, fmax_thresh=fmax, maxstep=maxstep,
                    header=f"[{i + 1}/{len(self.all_points)}] scan{point}:",
                )
                convergence_status.append((fmax_calc, converged))
            
            else:
                logger.warning(f"[{i + 1}/{len(self.all_points)}] scan{point} failed")
            
        traj.close()
        traj = Trajectory(f"{out_prefix}.traj")

        # --- Output results ---
        _, suffix = make_outpath(self.builder, multi_model=True)  # get basename and output suffix (e.g., .xyz, .gjf)
        try:
            do_ase_write(f"{out_prefix}{suffix}", traj, self.builder)
        except Exception as e:
            logger.error(f"Failed to save final structures: {e}")

        df = self._make_df_from_results(traj_path=f"{out_prefix}.traj",
                                        convergence_status=convergence_status)
       
        # save Gaussian log file compatible to GaussView
        write_gaussian_scan_log(traj_file=f"{out_prefix}.traj",
                                output_log=f"{out_prefix}_gv.log",
                                calclevel=self.builder.data.get("calculator", {}),
                                scan_blocks=self.scan_blocks)

        # --- plot results--- 
        # 1D scan
        if len(self.scan_blocks) == 1:
            plot_1d_scan(df, x_label="index", y_label="Delta E [kcal/mol]", out_prefix=out_prefix)
        # 2D scan
        elif len(self.scan_blocks) == 2:
            plot_2d_scan(df, x_label=df.columns[1], y_label=df.columns[2], out_prefix=out_prefix)
        else:
            logger.info("Scan visualization for N-dimensional scans (N > 2) is not implemented. Skipped.")

        logger.info(f"Scan terminated: {len(self.all_points)} points.")
        logger.info(f"Results written to: {out_prefix}.traj, .xyz, .csv")

        # --- vibration analysis (1D-only)  ---
        if len(self.scan_blocks) == 1:
            traj_final = Trajectory(f"{out_prefix}.traj")
            max_peak_index, peak_indices = detect_peaks_from_df(df=df, prominence=0.01, distance=None, outdir="maxima")

            # vibrational analysis for peak(s)
            if self.peak_vib.lower() == "highest":
                traj_indices = max_peak_index
            elif self.peak_vib.lower() == "all":
                traj_indices = peak_indices
            else:
                logger.info("Vibration analysis for peaks skipped.") 
                traj_indices = []
            
            if len(traj_indices) > 0:
                run_vibrations_for_traj(traj=traj_final, traj_indices=traj_indices, builder=self.builder,
                                        outdir="maxima", out_prefix="scan_final", suffix=suffix)
            else:
                logger.info("No peaks found for vibration analysis. Skipped.")
        else:
            logger.info("Peak detection has not been implemented for multi-dimensional scans. Skipped.")


class IndependentScanRunner(BaseScanRunner):
    """
    Each scan point is modeled independently from the same initial structure (Not gaussian-like).
    Suitable for parallel execution. No dependency between scan points.
    """
    def __init__(self, toml_path: str):
        """
        Initialize the IndependentScanRunner from a TOML input file.

        Parameters:
            toml_path (str): Path to the TOML input file.
        """
        super().__init__(toml_path)

    def generate_initial_structure(self, i: int, point: Tuple[float]) -> Atoms:
        atoms_i = self.atoms.copy()
        atoms_i.info.update(self.atoms.info)
        atoms_i.calc = self.calc_factory()  # Create a new calculator instance for each point

        # Apply scan geometry directly
        for (stype, atomid, _), val in zip(self.scan_blocks, point):
            match stype:
                case "distance":
                    atoms_i.set_distance(*atomid, val)
                case "angle":
                    atoms_i.set_angle(*atomid, val)
                case "dihedral":
                    atoms_i.set_dihedral(*atomid, val)

        atoms_i.set_constraint(self.generate_constraints(point))
        return atoms_i


class SequentialScanRunner(BaseScanRunner):
    """
    Each scan point is modeled from the previous point's optimized structure (gaussian-like).
    Suitable for smooth scans along a path (e.g., 1D or 2D scans in serpentine order).
    """

    def __init__(self, toml_path: str):
        super().__init__(toml_path)
        self.prev_atoms = self.atoms  # Start from the input structure

    def generate_initial_structure(self, i: int, point: Tuple[float]) -> Atoms:
        atoms_i = self.prev_atoms.copy()
        atoms_i.info.update(self.prev_atoms.info)
        atoms_i.calc = self.calc_factory()  # Create a new calculator instance for each point

        for (stype, atomid, _), val in zip(self.scan_blocks, point):
            match stype:
                case "distance":
                    atoms_i.set_distance(*atomid, val)
                case "angle":
                    atoms_i.set_angle(*atomid, val)
                case "dihedral":
                    atoms_i.set_dihedral(*atomid, val)

        atoms_i.set_constraint(self.generate_constraints(point))
        self.prev_atoms = atoms_i  # Update for next step
        return atoms_i


# --- main function ---
def run_scan(toml_path: str):
    """
    Run a multi-dimensional scan based on the provided TOML configuration.

    Parameters:
        toml_path (str): Path to the TOML file containing scan configuration.
    """
    with open(toml_path, "rb") as f:
        data = tomllib.load(f)

    scan_data = data.get("scan", {})
    scan_sequential = scan_data.get("sequential", True)

    distance_blocks = scan_data.get("distance", [])
    angle_blocks = scan_data.get("angle", [])
    dihedral_blocks = scan_data.get("dihedral", [])

    ndim = len(distance_blocks) + len(angle_blocks) + len(dihedral_blocks)

    if ndim > 3:
        raise ValueError("N-dimensional scan (N > 3) is not supported.")

    if scan_sequential:
        if ndim == 3:
            runner = IndependentScanRunner(toml_path)
            logger.info("3D scan is not supported in sequential mode."
                        f"Switching to runner = {runner.__class__.__name__}.")
        else:
            runner = SequentialScanRunner(toml_path)
            logger.info(f"Sequential relaxed scan selected: runner = {runner.__class__.__name__}.")
    else:
        runner = IndependentScanRunner(toml_path)
        logger.info(f"Independent relaxed scan selected: runner = {runner.__class__.__name__}.")

    runner.run(out_prefix="scan_final", log_prefix="scan", log_dir="opt_logs")


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python scan_runner.py input.toml")
        sys.exit(1)

    toml_path = sys.argv[1]

    # log setup
    logger.add("scan_runner.log", level="INFO", rotation="1 MB", backtrace=True, diagnose=True)
    logger.info(f"Running scan with input: {toml_path}")
    runner = IndependentScanRunner(toml_path)  # or SequentialScanRunner
    runner.run()
    logger.info("Scan calculation terminated normally.")
