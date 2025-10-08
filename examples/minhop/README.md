
## 🖥️ Minima-Hopping Optimization

`lcrx` supports Minima-Hopping optimization based on the ASE `MinimaHopping` class.

You can preview a template by running:

```sh
lcrx template minhop
```

### 📗 Example: Minima-Hopping
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
type = "minhop"             # Minima-Hopping calculation
Ediff0 = 6                  # Energy acceptance threshold (kcal/mol). New minimum is accepted if E(new) < E(previous) + Ediff0.
T0 = 500.0                  # in K, initial MD temperature. default: 1000 K
totalstep = 50              # total steps (cycles) for MD and subsequent optimization.
minima_threshold = 0.5      # in angstrom. Threshold for considering two minima identical, based on the maximum atomic displacement. default: 0.5 A.
timestep = 1.0              # in fs, timestep for MD simulations. default: 1.0 fs
mdmin = 4                   # criteria to stop MD simulation (no. of minima). Found minima is optimized. default: 2
method = "FIRE"             # local optimizer to use. or "BFGS", "LBFGS", "BFGSLINESEARCH"
fmax = 0.02                 # maximum force for convergence criterion (eV angstrom^-1)
```


### ⚙️ Input Parameters

Parameters are specified in the TOML file using the same variable names as in ASE’s MinimaHopping.
See the [Minima Hopping](https://ase-lib.org/ase/optimize.html#minima-hopping) section in the ASE documentation for details.

### 🧩 Constraints

Constraints can be applied by adding a `[constraints]` block to the TOML file.
See the [constraints section](../constraints/README.md) for details.


### ⚠️ Notes and Recommendations

The efficiency of Minima-Hopping for conformational searches strongly depends on the system and the chosen parameters.
**The termination of a calculation does not necessarily indicate that the global minimum has been located.**

We recommend carefully inspecting the results:

* If the conformational search appears insufficient, try running multiple replicas or increasing `totalstep`.
* To extend a completed calculation, increase `totalstep` in `input.toml` and re-run in the same directory.

For example, if a calculation finished with `totalstep = 10`, then updating it to `totalstep = 20` and re-running will perform an **additional 10 steps** starting from the previous results.

