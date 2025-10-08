from __future__ import annotations

# --- third-party imports ---
import pandas as pd
from ase.units import Hartree, kcal, mol

# ---- LocalColabReactionX imports ----
from ..utils.plot_utlis import make_2Dplot_from_df


def traj_to_csv(traj, out_csv_path="history.csv", return_df=False):
    indices = []
    energies = []
    try:
        for i, atoms in enumerate(traj):
            try:
                energy = atoms.get_potential_energy()
            except Exception:
                energy = None

            indices.append(i+1)
            energies.append((energy))

    except Exception as e:
        print(f"Warning: An error occurred while reading {traj}: {e}")

    # Save energies to CSV
    if out_csv_path is not None:
        df = pd.DataFrame({
            "index": indices,
            "energy [eV]": energies,
            "energy [hartree]": [e / Hartree for e in energies],
            "energy [kcal/mol]": [e / (kcal / mol) for e in energies],
            "Delta E [kcal/mol]": [(e - energies[0]) / (kcal / mol) for e in energies]
        })
        df.to_csv(out_csv_path, index=False)

        # Plot energies [kcal/mol]. read from df
        make_2Dplot_from_df(df=df,
                            df_xlabel="index",
                            df_ylabel="Delta E [kcal/mol]",
                            plot_xlabel="Index",
                            plot_ylabel=r"$\Delta E$ [kcal/mol]",
                            outpath=out_csv_path.replace('.csv', '.pdf'))

    if return_df:
        return df
    else:
        return None
