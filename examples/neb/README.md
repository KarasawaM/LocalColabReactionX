
## 🖥️ Nudged Elastic Band (NEB)

`LCRX` supports **double-ended transition state searches** based on the **Nudged Elastic Band (NEB)** method.  
See the [Nudged elastic band](https://ase-lib.org/ase/neb.html) section in the ASE documentation for details.

You can preview a template by running:
```sh
lcrx template neb
```

### 📗 Example: NEB Settings

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
type = "neb"                # Direct MaxFlux calculation
neb_method = "aseneb"       # or "improvedtangent",  "eb", "spline", "string".
n_images = 5                # total images including endpoints.
climb = true                # climbing image on/off. default: true
idpp = true                 # use IDPP interpolation. default: true
method = "FIRE"             # optimizer ("FIRE", "BFGS", "LBFGS", "BFGSLineSearch")
fmax = 0.05                 # maximum force for convergence criterion (eV angstrom^-1)
maxstep = 1000              # maximum number of optimization steps
parallel = false            # Enable parallel execution across images. CPU only, not available for GPU. default: false
peak_vibration = "highest"  # subsequent vibrational analysis for: "highest" energy peak, "all" peaks, or "none". default: "highest"
k = 0.1                     # (optional) NEB spring constant (eV angstrom^-2). default: 0.1
```

### 💡 Using CPU Calculators for Parallel Processing
When running with CPU calculators, set `parallel = true` in the `[calculation]` block to distribute the `n_images` across available CPU cores.
For more details on parallel execution, see the [Direct MaxFlux (DMF)](../dmf/README.md) section.


