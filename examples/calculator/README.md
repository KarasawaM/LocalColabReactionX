
## 🧩 Calculator

`LCRX` supports the following backends:
* **Fairchem UMA model**
* **TBLite (GFN-xTB)**
* **Gaussian**
* **ORCA**

You can preview a template by running:
```sh
lcrx template uma
```
Replace `uma` with `tblite`, `gaussian`, or `orca` to preview other examples.


### 📗 Example: UMA Calculator
See the [Fairchem GitHub repository](https://github.com/facebookresearch/fairchem) for details.
```toml
[calculator]
type = "UMA"                # Universal Models for Atoms
model_name = "uma-s-1p1"    # or "uma-m-1p1"
task = "omol"               # omol task for molecules.
device = "gpu"              # or "cpu"
```

### 📗 Example: TBLite Calculator
See the [TBLite documentation](https://tblite.readthedocs.io/en/latest/users/ase.html) for details.
```toml
[calculator]
type = "tblite"                   # Use the TBLite calculator
method = "GFN2-xTB"               # Semiempirical Tight-Binding method
# solvation = ["alpb", "water"]   # Optional: ["solvation method", "solvent name"]
totalcores = 16                   # total CPU cores to use.
totalmem = 60                     # in GB. total CPU memory to use.
```

### 📗 Example: Gaussian Calculator
To use Gaussian as the calculator, **ensure that the `g16` command is available in your PATH.**
See the [Gaussian calculator](https://ase-lib.org/ase/calculators/gaussian.html) in the ASE documentation for details.
```
[calculator]
type = "gaussian"       # Use Gaussian calculator
method = "B3LYP"        # Functional. e.g. "B3LYP", "M062X", "wB97XD".
basis = "6-31G(d)"      # Basis set. e.g. "6-31G(d)", "def2-SVP"
scf = ["xqc", "maxconventionalcycle=80"]  # scf(args1, args2). Comment out if not required.
extra = "nosymm"        # Other inputs. e.g. "nosymm SCRF=... EmpiricalDispersion=..." Comment out if not required.
totalcores = 16         # total CPU cores to use.
totalmem = 60           # in GB. total CPU memory to use.
```

The above settings will create the following `g16.com` input file:
```
%mem=60GB
%chk=g16.chk
%nprocshared=16
#p B3LYP/6-31G(d) ! ASE formatted method and basis
scf(xqc,maxconventionalcycle=80)
nosymm
force

Gaussian input prepared by ASE

(charge multiplicity)
(coordinates)

```
Charge, multiplicity, and coordinates are read from the `[structure]` block.


### 📗 Example: ORCA Calculator
To use ORCA as the calculator, **ensure that `orca` is installed.**
See the [ORCA calculator](https://ase-lib.org/ase/calculators/orca.html#module-ase.calculators.orca) in the ASE documentation for details.

```toml
[calculator]
type = "orca"                               # Use ORCA calculator
orcapath = "/full/path/to/orca"             # Full path to the ORCA executable.
orcasimpleinput = "B3LYP 6-31G(d) ENGRAD"   # ORCA simple input string WITHOUT !. Set gradient type (EnGRAD or NumGrad).
totalcores = 16                             # total CPU cores to use.
totalmem = 60                               # in GB. total CPU memory to use.

orcablocks = """
# Add ORCA blocks below if needed, such as scf, geom etc...
# PAL and Maxcores are controlled based on totalcores and totalmem.
"""
```

The above settings will generate the following `orca.inp` input file:
```
! B3LYP 6-31G(d) ENGRAD 
# Add ORCA blocks below if needed, such as scf, geom etc...
# PAL and Maxcores are controlled by totalmem and totalcores.
%pal nprocs 16 end
%maxcore 3700

*xyz (charge multiplicity)
(coordinates)
*
```
Charge, multiplicity, and coordinates are read from the `[structure]` block.

