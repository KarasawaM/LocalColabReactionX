
## 🖥️ Relaxed Scan

`LCRX` supports **relaxed scan calculations** using the ASE methods `Atoms.set_distance()`, `Atoms.set_angle()`, `Atoms.set_dihedral()`, and the `FixInternals` class.

You can preview a template by running:
```sh
lcrx template scan
```


### 📗 Example: Scan Settings

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
type = "scan"               # Scan calculation
fmax = 0.02                 # maximum force for convergence criterion (eV/Å)
maxstep = 1000              # maximum number of optimization steps
method = "FIRE"             # or "BFGS", "LBFGS", "BFGSLINESEARCH"

[scan]  # Start the scan section            
sequential = true  # or false

# 1D scan: set one [[scan.XXX]], 2D scan: set two [[scan.XXX]]
[[scan.distance]]
atomid = "1, 2"         # or [1, 2]. 1-based indices for target atoms.
value = "1.5, 1.0, 6"   # or [1.5, 1.0, 6]. Start, End, Nsteps (in angstrom).

[[scan.angle]]
atomid = "3, 4, 5"
value = "120, 180, 5"   # Start, End, Nsteps (in degrees).

[[scan.dihedral]]
atomid = "6, 7, 8, 9"
value = "90, 180, 10"   # Start, End, Nsteps (in degrees).
```

The argument **`sequential`** controls how the initial structure is prepared for each scan point:

* `true`: the optimized structure from step *N–1* is used as the starting point for step *N*.
* `false`: each step is initialized independently from the original input coordinates.

In most cases, setting **`sequential = true` provides better convergence.**


### 1️⃣ Running One-Dimensional Relaxed Scans

Specify a single `[[scan.XXX]]` block.
For 1D scans, automatic **peak detection** of the energy profile and **vibrational analysis at the maxima** are performed.


### 2️⃣ Running Two-Dimensional Relaxed Scans

Specify two separate `[[scan.XXX]]` blocks.
A 2D scan grid is constructed from the two sets of scan values, and each grid point is evaluated in a **serpentine order**.

### 🔍 Checking Scan Values (optional)

You can check the actual scan values in Python as follows:

```sh
# value = "1.5, 1.0, 6"
python3
>>> import numpy as np
>>> np.linspace(1.5, 1.0, 6)
array([1.5, 1.4, 1.3, 1.2, 1.1, 1. ])
```

In this distance scan example, the bond length is scanned from 1.5 Å down to 1.0 Å in steps of 0.1 Å.


### 🧩 Constraints

Constraints can be applied by adding a `[constraints]` block to the TOML file.
See the [constraints section](../constraints/README.md) for details.

