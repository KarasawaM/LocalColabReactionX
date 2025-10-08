from __future__ import annotations

from functools import partial
from ase.optimize import BFGS, FIRE, FIRE2, LBFGS, BFGSLineSearch, LBFGSLineSearch, GPMin, MDMin


def get_optimizer(method: str):
    """
    Return a callable constructor for the requested optimizer.
    For most methods this is the class itself; for ABC-FIRE it is a partial
    of FIRE2 with use_abc=True, so callers can keep doing Opt(atoms, ...).
    """
    key = method.strip().upper().replace("-", "").replace("_", "")  # remove - and _.
    match key:
        case "FIRE":
            return FIRE
        case "FIRE2":
            return FIRE2
        case "ABCFIRE" | "FIRE2ABC" | "ABC":  # aliases for convenience.
            return partial(FIRE2, use_abc=True)  # ABC-FIRE is available by use_abc=True.
        case "BFGS":
            return BFGS
        case "LBFGS":
            return LBFGS
        case "BFGSLINESEARCH" | "QUASINEWTON":
            return BFGSLineSearch
        case "LBFGSLINESEARCH":
            return LBFGSLineSearch
        case "GPMIN":
            return GPMin
        case "MDMIN":
            return MDMin
        case _:
            raise ValueError(f"Unsupported optimization method: {method}")

