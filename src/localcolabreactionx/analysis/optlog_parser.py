from __future__ import annotations

# --- general imports ---
import re

# --- third-party modules ---
from loguru import logger
from ase.units import Hartree


def _parse_opt_last_line(line: str) -> tuple:
    """
    Accepts both style:
            Step     Time          Energy          fmax
    FIRE:    0 19:17:30    -1611.208033        0.674131

                    Step[ FC]     Time          Energy          fmax
    BFGSLineSearch:    8[ 16] 12:14:46   -21265.240101       0.0487
      
    Returns: (step, energy_eV, fmax)    
    """
    TIME_RE = re.compile(r'^\d{2}:\d{2}:\d{2}$')
    tokens = line.split()

    # time index 
    for i, token in enumerate(tokens):
        if TIME_RE.match(token):
            time_idx = i
            break

    # energy: time_idx + 1, fmax: time_idx + 2
    energy = float(tokens[time_idx + 1])
    fmax = float(tokens[time_idx + 2])

    # step.
    left_text = " ".join(tokens[:time_idx])

    m = re.search(r'(\d+)\s*\[', left_text)  # ex: '100[102]', '0[  0]'
    if m:
        step = int(m.group(1))
    else:
        # without []
        m = re.search(r'(\d+)\s*$', left_text)
        if not m:
            ints = re.findall(r'\d+', left_text)
            step = int(ints[-1]) if ints else 0
        else:
            step = int(m.group(1))

    return step, energy, fmax


def check_convergence_from_log(logfile_path: str, header: str, fmax_thresh: float, maxstep: int) -> tuple:
    """
    Check convergence status from an ASE optimizer log file.

    Parameters:
        logfile_path (str): Path to the log file.
        header (str): Header message.

    Returns:
        fmax and True if converged, False otherwise.
    """
    try:
        with open(logfile_path, "r") as f:
            lines = f.readlines()

        # Find the last non-empty line
        last_line = ""
        for line in reversed(lines):
            if line.strip():
                last_line = line.strip()
                break
    
        step, energy, fmax = _parse_opt_last_line(last_line)
        Eh = energy / Hartree  # Convert eV to Hartree

        # convergence check
        if step <= maxstep and fmax <= fmax_thresh:
            converged = "Yes"
        else:
            converged = "No"

        # log text
        message = [f"{header}",
                f"Eh = {Eh:.6f} Hartree.",
                f"Converged: {converged} (fmax = {fmax:.4f})"]

        # write log
        if converged == "Yes":
            logger.info(" ".join(message))
        else:
            logger.warning(" ".join(message))

        return fmax, converged

    except Exception as e:
        logger.warning(f"Failed to parse {logfile_path}: {e}")


