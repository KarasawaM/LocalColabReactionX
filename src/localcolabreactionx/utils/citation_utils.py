# citation_utils.py
from collections import OrderedDict
from typing import Dict, Iterable
from loguru import logger

# === FULL Citation DB: key -> plain text ==========================
CITATIONS_FULL: Dict[str, str] = {
    # --- general ---
    "ColabReaction": "Karasawa, M.; Leow, C. S.; Yajima, H.; Arai, S.; Nishizaki, H.; Terada, T.; Sato, H. ColabReaction: Accelerating Transition State Searches with Machine Learning Potentials on Google Colaboratory. ChemRxiv 2025. DOI: 10.26434/chemrxiv-2025-zvkqk.",
    "ASE": "Hjorth Larsen, A.; Jørgen Mortensen, J.; Blomqvist, J.; Castelli, I. E.; Christensen, R.; Dułak, M.; Friis, J.; Groves, M. N.; Hammer, B.; Hargus, C.; Hermes, E. D.; Jennings, P. C.; Bjerre Jensen, P.; Kermode, J.; Kitchin, J. R.; Leonhard Kolsbjerg, E.; Kubal, J.; Kaasbjerg, K.; Lysgaard, S.; Bergmann Maronsson, J.; Maxson, T.; Olsen, T.; Pastewka, L.; Peterson, A.; Rostgaard, C.; Schiøtz, J.; Schütt, O.; Strange, M.; Thygesen, K. S.; Vegge, T.; Vilhelmsen, L.; Walter, M.; Zeng, Z.; Jacobsen, K. W. The atomic simulation environment—a Python library for working with atoms. Journal of Physics: Condensed Matter 2017, 29 (27), 273002. DOI: 10.1088/1361-648X/aa680e.",

    # --- calculator / models ---
    "UMA": "Wood, B. M.; Dzamba, M.; Fu, X.; Gao, M.; Shuaibi, M.; Barroso-Luque, L.; Abdelmaqsoud, K.; Gharakhanyan, V.; Kitchin, J. R.; Levine, D. S.; Michel, K.; Sriram, A.; Cohen, T.; Das, A.; Rizvi, A.; Sahoo, S. J.; Ulissi, Z. W.; Zitnick, C. L. UMA: A Family of Universal Models for Atoms. arXiv preprint 2025, arXiv:2506.23971 [cs.LG]. DOI: 10.48550/arXiv.2506.23971.",
    "OMol25": "Levine, D. S.; Shuaibi, M.; Spotte-Smith, E. W. C.; Taylor, M. G.; Hasyim, M. R.; Michel, K.; Batatia, I.; Csányi, G.; Dzamba, M.; Eastman, P.; Frey, N. C.; Fu, X.; Gharakhanyan, V.; Krishnapriyan, A. S.; Rackers, J. A.; Raja, S.; Rizvi, A.; Rosen, A. S.; Ulissi, Z.; Vargas, S.; Zitnick, C. L.; Blau, S. M.; Wood, B. M. The Open Molecules 2025 (OMol25) Dataset, Evaluations, and Models. arXiv preprint 2025, arXiv:2505.08762 [physics.chem-ph]. DOI: 10.48550/arXiv.2505.08762.",
    "Fairchem": "fairchem; https://github.com/facebookresearch/fairchem",
    "TBLite": "TBLite; https://github.com/tblite/tblite",
    "GFN2-xTB": "Bannwarth, C.; Ehlert, S.; Grimme, S. GFN2-xTB—An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions. Journal of Chemical Theory and Computation 2019, 15 (3), 1652-1671. DOI: 10.1021/acs.jctc.8b01176.",

    # --- algorithms / methods ---
    "Minima Hopping": "Goedecker, S. Minima hopping: An efficient search method for the global minimum of the potential energy surface of complex molecular systems. The Journal of Chemical Physics 2004, 120 (21), 9911-9917. DOI: 10.1063/1.1724816.",
    "DMF": "Koda, S.; Saito, S. Locating Transition States by Variational Reaction Path Optimization with an Energy-Derivative-Free Objective Function. Journal of Chemical Theory and Computation 2024, 20 (7), 2798-2811. DOI: 10.1021/acs.jctc.3c01246.",
    "FBENM": "Koda, S.; Saito, S. Flat-Bottom Elastic Network Model for Generating Improved Plausible Reaction Paths. Journal of Chemical Theory and Computation 2024, 20 (16), 7176-7187. DOI: 10.1021/acs.jctc.4c00792.",
    "CFBENM": "Koda, S.; Saito, S. Correlated Flat-Bottom Elastic Network Model for Improved Bond Rearrangement in Reaction Paths. Journal of Chemical Theory and Computation 2025, 21 (7), 3513-3522. DOI: 10.1021/acs.jctc.4c01549.",
    "DMF/UMA": "Nakano, M.; Karasawa, M.; Ohmura, T.; Terada, T.; Sato, H. High-Speed Terpene Pathway Analysis via Double-ended Transition State Search Method Combined with ML Potentials Uncovers the Intricate Rearrangement Reaction in Spiroalbatene Biosynthesis. ChemRxiv 2025. DOI: 10.26434/chemrxiv-2025-md8k6-v2.",
    "FIRE": "Bitzek, E.; Koskinen, P.; Gähler, F.; Moseler, M.; Gumbsch, P. Structural Relaxation Made Simple. Physical Review Letters 2006, 97 (17), 170201. DOI: 10.1103/PhysRevLett.97.170201.",
}



# === META Citation DB: key -> short text ==========================
CITATIONS_META: Dict[str, str] = {
    # --- general ---
    "ColabReaction": "ChemRxiv 2025, DOI: 10.26434/chemrxiv-2025-zvkqk.",
    "ASE": "J. Phys. Condens. Matter 2017, DOI: 10.1088/1361-648X/aa680e.",

    # --- calculator / models ---
    "UMA": "arXiv preprint 2025, DOI: 10.48550/arXiv.2506.23971.",
    "OMol25": "arXiv preprint 2025, DOI: 10.48550/arXiv.2505.08762.",
    "Fairchem": "Fairchem, https://github.com/facebookresearch/fairchem",
    "TBLite": "TBLite, https://github.com/tblite/tblite",
    "GFN2-xTB": "JCTC 2019, DOI: 10.1021/acs.jctc.8b01176.",

    # --- algorithms / methods ---
    "Minima Hopping": "J. Chem. Phys. 2004, DOI: 10.1063/1.1724816.",
    "DMF": "JCTC 2024, DOI: 10.1021/acs.jctc.3c01246.",
    "FBENM": "JCTC 2024, DOI: 10.1021/acs.jctc.4c00792.",
    "CFBENM": "JCTC 2025, DOI: 10.1021/acs.jctc.4c01549.",
    "DMF/UMA": "ChemRxiv 2025, DOI: 10.26434/chemrxiv-2025-md8k6-v2.",
    "FIRE": "PRL 2006, DOI: 10.1103/PhysRevLett.97.170201.",
}

# === rules: TOML for keys =====================
# calculator.type → citation keys
CALCULATOR_TO_CITES = {
    "uma": ["UMA", "OMol25", "Fairchem"],
    "uma-s-1p1": ["UMA", "OMol25", "Fairchem"],
    "uma-m-1p1": ["UMA", "OMol25", "Fairchem"],
    "tblite": ["TBLite", "GFN2-xTB"],
    "xtb": ["TBLite", "GFN2-xTB"],
    "emt": [],  # None
}

# calculation.type → citation key
CALCULATION_TO_CITES = {
    "dmf": ["DMF", "FBENM", "CFBENM"],
    "minhop": ["Minima Hopping"],
    "opt": [],                   
    "scan": [],                 
}

METHOD_TO_CITES = {
    "FIRE": ["FIRE"],
    "BFGS": [],
    "LBFGS": [],
    "BFGSLINESEARCH": [],
}

# === utils =====================================================
def _append_unique(seq: OrderedDict[str, None], keys: Iterable[str]) -> None:
    """Add keys"""
    for k in keys:
        if k and k not in seq:
            seq[k] = None

def collect_citation_keys_from_toml(data: dict) -> list[str]:
    keys = OrderedDict()  # insertion-order set

    # general
    _append_unique(keys, ["ASE", "ColabReaction"])

    calcblk = (data.get("calculation") or {})
    calctype = str(calcblk.get("type", "")).lower()

    calc = (data.get("calculator") or {})
    calc_type = str(calc.get("type", "")).lower()

    # calculators
    _append_unique(keys, CALCULATOR_TO_CITES.get(calc_type, []))

    # calculation.type
    _append_unique(keys, CALCULATION_TO_CITES.get(calctype, []))

    # --- DMF/UMA  ---
    if calctype == "dmf" and calc_type.startswith("uma"):
        _append_unique(keys, ["DMF/UMA"])

    # methods
    method = str(calcblk.get("method", "")).upper()
    _append_unique(keys, METHOD_TO_CITES.get(method, []))

    return list(keys.keys())


def log_citations(keys: Iterable[str], verbose: bool = False) -> None:
    """
    Log citations depending on CLI --verbose flag.
    - verbose=True  -> CITATIONS_FULL
    - verbose=False -> CITATIONS_META
    """
    db = CITATIONS_FULL if verbose else CITATIONS_META
    keys = [k for k in keys if k in db]
    if not keys:
        return

    fmt = "complete" if verbose else "condensed"
    lines = [f"Recommended citations ({fmt} format):"]

    # Longest "[KEY]" width
    maxlen = max(len(f"[{k}]") for k in keys)

    for k in keys:
        key_str = f"[{k}]"
        # pad to maxlen + 2 spaces
        lines.append(f"{key_str:<{maxlen + 2}}{db[k]}")

    logger.info("\n".join(lines))


def collect_and_log_citations(data: dict, verbose: bool = False) -> None:
    keys = collect_citation_keys_from_toml(data)
    log_citations(keys, verbose=verbose)