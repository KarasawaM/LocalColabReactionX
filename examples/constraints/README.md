
## 🧩 Constraints

`LCRX` supports **constraitns** using the ASE `FixAtoms` and `FixInternal` class. 
Constraints can be applied to both **optimization tasks and molecular dynamics (MD) simulations**.

You can preview a template by running:

```sh
lcrx template constraints
```

### 📗 Example: Constraints Settings
**Positional constraint (FixAtoms in ASE)**
```toml
[constraints]  # Start the constraints section
[[constraints.freeze]]
atomid = "1, 2, 3"      # or [1, 2, 3]. 1-based indices for target atoms.
```

**Distance constraint (FixInternals in ASE)**
```toml
[constraints]  # Start the constraints section
[[constraints.distance]]
atomid = "4, 5"         # or [4, 5]. 1-based indices for target atoms.
value = 1.5             # in angstrom
```

**Angle constraint (FixInternals in ASE)**
```toml
[constraints]  # Start the constraints section
[[constraints.angle]]
atomid = "6, 7, 8"      # or [6, 7, 8]. 1-based indices for target atoms.
value = 120             # in degrees
```

**Dihedral constraint (FixInternals in ASE)**
```toml
[constraints]  # Start the constraints section
[[constraints.dihedral]]
atomid = "9, 10, 11, 12"  # or [9, 10, 11, 12]. 1-based indices for target atoms.
value = 180               # in degrees
```

**Hookean restraint (Hookean in ASE)**
See the [Hookean class](https://ase-lib.org/ase/constraints.html#the-hookean-class) in the ASE documentation for details.  
```toml
[constraints]  # Start the constraints section
[[constraints.hookean]]
atomid = "13, 14"       # or [13, 14]. 1-based indices for target atoms.
k = 10                  # Hooke's law (spring) constant (kcal mol^-1 angstrom^-2)
r0 = 1.5                # in angstrom. When distance > r0, apply harmonic restraint.
```


### 🔗 Multiple Constraints
You can combine multiple constraints in the same input file.
For example, to apply distance constraints to two different bonds, include two separate `[[constraints.distance]]` blocks:

```toml
[constraints]  # Start the constraints section

[[constraints.distance]]
atomid = "4, 5"         # or [4, 5]. 1-based indices for target atoms.
value = 1.5             # in angstrom

[[constraints.distance]]
atomid = "6, 7"         # or [6, 7]. 1-based indices for target atoms.
value = 2.0             # in angstrom
```
