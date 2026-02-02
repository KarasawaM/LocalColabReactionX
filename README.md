# LCRX: LocalColabReaction eXtended (beta)

A local environment version of **[ColabReaction](https://github.com/BILAB/ColabReaction)** with extended functionalities.

**LocalColabReaction eXtended (LCRX)**:

* Provides a **command-line interface**, making it easy to run on clusters and HPC systems.  
* Supports a wide range of calculation types, including:  
  * **Structure optimization**  
  * **Relaxed scans**
  * **Single- and double-ended transition state searches**  
  * **Molecular dynamics (MD) simulations**

⚠️ **Note:** NVE molecular dynamics is implemented but still under evaluation. Further improvements are planned for future releases.

---

## 🛠️ Installation

Dependencies such as `cyipopt` and `tblite` are most easily installed using **conda**.

We recommend creating a dedicated environment:

```sh
lcrx_env_name="lcrx"  # conda env name
conda create -n ${lcrx_env_name} -y -c conda-forge \
  python=3.12 cyipopt tblite tblite-python && \
conda activate ${lcrx_env_name}
```

Then install **LCRX**:

```sh
git clone https://github.com/KarasawaM/LocalColabReactionX.git
cd LocalColabReactionX
pip install .
```

### macOS users

- `LCRX` installs successfully on macOS, but note that the available `fairchem-core` build does **not** support the Universal Model for Atoms (UMA).
- The extended tight-binding method (xTB) is available through `TBLite`.
- Install `jax` by `conda install jax -c conda-forge`.

---

## 🔑 Using UMA Models

💡 **If you already have access to the UMA model, you can skip this section.**

To use the **UMA** machine learning interatomic potentials, you will need a Hugging Face access token.

For detailed instructions, see:
* [User Guide (English)](https://github.com/BILAB/ColabReaction/blob/main/docs/User_Guide_EN.pdf)
* [ユーザーガイド (Japanese)](https://github.com/BILAB/ColabReaction/blob/main/docs/User_Guide_JP.pdf)

Log in from the terminal with:
```sh
hf auth login
```

You will be prompted to enter your token:
```sh
Enter your token (input will not be visible):
```

After entering the token, confirm your login with:
```sh
hf auth whoami
```

Your Hugging Face account name should be displayed if the login was successful.

---

## 🖥️ Running Calculations

```sh
lcrx run --toml <input>.toml
```

All calculation settings are specified in `<input>.toml` (e.g. `opt.toml`, `dmf.toml`).

Currently supported calculation types include:
* **Geometry optimization**
* **Relaxed scan**
* **Direct MaxFlux (DMF)** — double-ended TS search
* **Nudged Elastic Band (NEB)** — double-ended TS search
* **Transition state optimization + IRC**
* **Minima-Hopping** — global optimization via MD + local relaxation

All workflows are executed through the Atomic Simulation Environment (ASE).

---

## ⚙️ Input Files in TOML Format

You can view example TOML input files with:

```sh
lcrx template opt
```

Replace `opt` with `scan`, `dmf`, `neb`, `ts`, or `minhop` to preview other examples.

To save an example to a file, run:
```sh
lcrx template opt -o
```

This will create `template_opt.toml`.

### 🗂️ TOML Structure

* **[structure]** — input file, charge, and multiplicity
* **[calculator]** — calculator or potential (UMA, TBLite, Gaussian, ORCA)
* **[calculation]** — calculation type (opt, scan, dmf, minhop, …) and parameters

### 📗 Example: Geometry Optimization (OPT)

```toml
[structure]
input = "input.xyz"         # or .pdb, .gjf, .sdf
charge = 0                  # total charge
mult = 1                    # multiplicity (2S+1)

[calculator]
type = "UMA"                # Universal Models for Atoms
model_name = "uma-s-1p1"    # or "uma-m-1p1"
task = "omol"               # omol task for molecules.
device = "gpu"              # or "cpu"

[calculation]
type = "opt"        # Geometry optimization
fmax = 0.02         # maximum force for convergence criterion (eV angstrom^-1)
maxstep = 1000      # maximum number of optimization steps
method = "FIRE"     # or "BFGS", "LBFGS", "BFGSLINESEARCH"
```

### 📗 Example: Direct MaxFlux (DMF)

```toml
[structure]
reactant = "react.xyz"      # or .pdb, .gjf, .sdf
product  = "prod.xyz"
charge = 0                  # total charge
mult = 1                    # multiplicity (2S+1)

[calculator]
type = "UMA"                # Universal Models for Atoms
model_name = "uma-s-1p1"    # or "uma-m-1p1"
task = "omol"               # omol task for molecules.
device = "gpu"              # or "cpu"

[calculation]
type = "dmf"                # Direct MaxFlux calculation
dmf_nmove = 20              # nmove for path optimization. default: 5.
update_teval = false        # or true. default: false
dmf_convergence = "tight"   # or "middle", "loose". default: "tight"
```

See the `examples` directory for additional templates and details:
- [Structure optimization](examples/opt/README.md)
- [Relaxed scan](examples/scan/README.md)
- [Direct MaxFlux (DMF)](examples/dmf/README.md)
- [Nudged Elastic Band (NEB)](examples/neb/README.md)
- [Transition state optimization + IRC](examples/ts_irc/README.md)
- [Minima-Hopping](examples/minhop/README.md)

### 🧩 Constraints

Constraints can be applied to both optimization tasks and molecular dynamics (MD) simulations.

You can preview a template by running:
```sh
lcrx template constraints
```

See the [constraints section](examples/constraints/README.md) for details.

### 🧩 Calculator Backends

`lcrx` supports the following backends:
* **Fairchem UMA model**
* **TBLite (GFN-xTB)**
* **Gaussian**
* **ORCA**

To preview a specific calculator template:
```sh
lcrx template uma
```

See the [calculator section](examples/calculator/README.md) for details.

---

## 🗒️ Running with Shell Scripts

See the [Shell Scripts section](examples/sh/README.md) for details.

---

## License

`LCRX` itself is released under the MIT License.

⚠️ **License Notice**  
- This project **does NOT include or distribute any third-party programs or model weights that require separate licenses**, such as Gaussian, ORCA, or UMA models.  
- Users must obtain these components independently and ensure compliance with the respective license agreements.
