
## 🖥️ Direct MaxFlux (DMF)

`LCRX` supports **double-ended transition state searches** based on the **Direct MaxFlux (DMF)** method.  
See the [DMF GitHub repository](https://github.com/shin1koda/dmf) for details.

You can preview a template by running:
```sh
lcrx template dmf
```

### 📗 Example: DMF Settings

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

# --- initial path ---
fbenm = "cfbenm"            # or "fbenm". "cfbenm" for correlated = true. default: "cfbenm"
fbenm_nmove = 10            # nmove for initial path estimation. default: 10

# --- path optimization ---
dmf_nmove = 20              # nmove for path optimization. default: 5.
update_teval = false        # or true. default: false
dmf_convergence = "tight"   # or "middle", "loose". default: "tight"
peak_vibration = "highest"  # subsequent vibrational analysis for: "highest" energy peak, "all" peaks, or "none". default: "highest"
```

### 💡 Using CPU Calculators for Parallel Processing
When running with CPU calculators, set `parallel = true` in the `[calculation]` block to distribute the `dmf_nmoves` across available CPU cores.

```toml
[structure]
reactant = "react.xyz"      # or .pdb, .gjf, .sdf
product  = "prod.xyz"
charge = 0                  # total charge
mult = 1                    # multiplicity (2S+1)

[calculator]
type = "gaussian"       # Use Gaussian calculator
method = "B3LYP"        # Functional. e.g. "B3LYP", "M062X", "wB97XD"
basis = "6-31G(d)"      # Basis set. e.g. "6-31G(d)", "def2-SVP"
scf = ["xqc", "maxconventionalcycle=80"]  # scf(args1, args2). Comment out if not required.
extra = "nosymm"        # Other inputs. e.g. "nosymm SCRF=... EmpiricalDispersion=..." Comment out if not required.
totalcores = 30         # total CPU cores to use.
totalmem = 60           # in GB. total CPU memory to use.

[calculation]
type = "dmf"                # Direct MaxFlux calculation
dmf_nmove = 5               # nmove for path optimization. default: 5.
update_teval = ture         # or true. default: false
dmf_convergence = "tight"   # or "middle", "loose". default: "tight"
parallel = true             # Enable parallel execution across nmoves. CPU only, not available for GPU. default: false
peak_vibration = "highest"  # subsequent vibrational analysis for: "highest" energy peak, "all" peaks, or "none". default: "highest"
```

The key arguments for parallel calculations are `totalcores`, `totalmem`, and `dmf_nmove`.
In the example above, 30 cores and 60 GB of memory are allocated, and `dmf_nmove = 5` is distributed in parallel.
Thus, each evaluation point receives **6 CPU cores and 12 GB of memory**.

⚠️ **Note:** Ensure that `totalcores` and `totalmem` do not exceed the physical resources available on your machine.


