from __future__ import annotations
import numpy as np
import pandas as pd
from ase.atoms import Atoms
from ase.units import Hartree, Bohr, kcal, mol
from ase.io import read, write, Trajectory
from loguru import logger


def concat_irc_traj(irc_f: str, irc_r: str, prefix: str = "IRC"):
    """
    irc_f: forward direction .traj file
    irc_r: reverse direction .traj file

    Outputs:
      - IRC_FtoR.traj : Foward  to Reverse
      - IRC_RtoF.traj : Reverse to Forward
    """

    f_frames = read(irc_f, index=":")
    if not isinstance(f_frames, list):
        f_frames = [f_frames]
        
    r_frames = read(irc_r, index=":")
    if not isinstance(r_frames, list):
        r_frames = [r_frames]

    if len(f_frames) < 2:
        raise ValueError("IRC-F must have at least 2 frames")

    # reverse IRC-F frames, drop last (first TS)
    f_rev = list(reversed(f_frames))[:-1]

    # --- FtoR (Forward to Reverse) ---
    combined_FtoR = f_rev + r_frames
    with Trajectory("IRC_merged/IRC_FtoR.traj", "w") as traj:
        for frame in combined_FtoR:
            traj.write(frame)

    traj = read("IRC_merged/IRC_FtoR.traj", ":")
    write("IRC_merged/IRC_FtoR.xyz", traj, plain=True)  # write XYZ 

    # --- RtoF (Reverse to Forward) ---
    combined_RtoF = list(reversed(combined_FtoR))
    with Trajectory("IRC_merged/IRC_RtoF.traj", "w") as traj:
        for frame in combined_RtoF:
            traj.write(frame)

    traj = read("IRC_merged/IRC_RtoF.traj", ":")
    write("IRC_merged/IRC_RtoF.xyz", traj, plain=True)    # write XYZ

    logger.info(f"IRC {len(combined_FtoR)} points written to: IRC_merged/IRC_FtoR.traj (Fwd -> Rev), IRC_merged/IRC_RtoF.traj (Rev -> Fwd)")
    logger.info("Directions are arbitrary. Please choose correct reaction path (FtoR or RtoF).")


def _calc_ircstep_distance(a: Atoms, b: Atoms) -> tuple[float, float]:
    """
    a, b: Atoms
    - Δs_mw = sqrt( Σ_i m_i * ||r_i(b) - r_i(a)||^2 )   [units: sqrt(amu) * Å]
    - Δs_plain = sqrt( Σ_i ||r_i(b) - r_i(a)||^2 )       [units: Å]

    Returns
    -------
    (mw_step, plain_step): tuple[float, float]
        mw_step : (sqrt(amu) * Å)
        plain_step: (Å)
    """
    pa = a.get_positions()  # Å
    pb = b.get_positions()  # Å
    if pa.shape != pb.shape:
        raise ValueError("Frame shapes differ (atom count or ordering mismatch).")

    ma = a.get_masses()     # amu
    mb = b.get_masses()     # amu
    if not np.allclose(ma, mb):
        raise ValueError("Atomic masses differ between frames.")

    dr = pb - pa  # Å, shape (N, 3)

    # mass-weighted: Σ_i m_i * (dr_i · dr_i)
    mw_dist = np.einsum("i,ij,ij->", ma, dr, dr)      # (amu) * (Å^2)
    # plain: Σ_i (dr_i · dr_i)
    plain_dist = np.einsum("ij,ij->", dr, dr)         # (Å^2)

    return float(np.sqrt(mw_dist)), float(np.sqrt(plain_dist))


def make_df_oneside_irc(traj):
    traj = read(traj, ":")  # list[Atoms]

    if len(traj) < 2:
        raise ValueError("At least 2 frames required in the trajectory.")

    data_rows = []
    net_rxcoords_mw = 0.0  # [amu^1/2 Angstrom]
    net_rxcoords = 0.0     # [Angstrom]

    for i in range(len(traj)):
        atoms = traj[i]
        e_ev = atoms.get_potential_energy()  # eV

        e_kcalmol = e_ev / (kcal / mol)  # kcal/mol
        Eh = e_ev / Hartree

        if i == 0:
            data_rows.append((net_rxcoords_mw, net_rxcoords, Eh, e_kcalmol))
        else:
            ds_mw, ds_plain = _calc_ircstep_distance(traj[i - 1], traj[i])   # √amu * Å
            net_rxcoords_mw += ds_mw
            net_rxcoords += ds_plain
            data_rows.append((net_rxcoords_mw, net_rxcoords, Eh, e_kcalmol))

    # Print table to df
    df = pd.DataFrame(data_rows, columns=["IRC [amu^1/2 Angstrom]", "IRC [Angstrom]", "Energy [hartree]", "Energy [kcal/mol]"])

    return df


def merge_df_irc(df_minus, df_plus):
    """
    Merge forward and reverse IRC dataframes.
    IRC minus -> plus

    df_minus: intrinsic reaction coords will be minus.
    df_plus: intrinsic reaction coords will be plus.
    """
    # df_minus: Reaction Coordinate will be negative
    df_m = df_minus.copy()
    df_m["IRC [amu^1/2 Angstrom]"] *= -1.0
    df_m["IRC [Angstrom]"] *= -1.0

    # remove 0.0 point from df_m
    df_m = df_m[df_m["IRC [amu^1/2 Angstrom]"] < 0.0]

    # concat df_m and df_plus
    df_all = pd.concat([df_m, df_plus], ignore_index=True)

    # sort by IRC
    df_all = df_all.sort_values(by="IRC [amu^1/2 Angstrom]")
    df_all.reset_index(drop=True, inplace=True)

    # Add IRC in Bohr
    df_all["IRC [amu^1/2 Bohr]"] = df_all["IRC [amu^1/2 Angstrom]"] / Bohr
    df_all["IRC [Bohr]"] = df_all["IRC [Angstrom]"] / Bohr

    # Reorder columns
    df_all = df_all[["IRC [amu^1/2 Angstrom]", "IRC [amu^1/2 Bohr]", "IRC [Angstrom]", "IRC [Bohr]", "Energy [hartree]", "Energy [kcal/mol]"]]

    # TS energy reference. IRC [amu^1/2 Angstrom] = 0.0
    df_all["DeltaE [kcal/mol]"] = df_all["Energy [kcal/mol]"] - df_all.loc[df_all["IRC [amu^1/2 Angstrom]"] == 0.0, "Energy [kcal/mol]"].values[0]

    return df_all


