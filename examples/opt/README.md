
## 🖥️ Structure Optimization

`LCRX` supports **structure optimization** based on the ASE `Optimizer` class.

You can preview a template by running:
```sh
lcrx template opt
```

### 📗 Example: Optimization Settings

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

#### Optimization methods
See the [Local optimization](https://ase-lib.org/ase/optimize.html#local-optimization) section in the ASE documentation for details.


### 🧩 Constraints

Constraints can be applied by adding a `[constraints]` block to the TOML file.
See the [constraints section](../constraints/README.md) for details.

