import textwrap


# ----------------------
# structure blocks
# ----------------------
def structure_single():
    return textwrap.dedent('''\
    [structure]
    input = "input.xyz"         # or .pdb, .gjf, .sdf
    charge = 0                  # total charge
    mult = 1                    # multiplicity (2S+1)
    ''')


def structure_double():
    return textwrap.dedent('''\
    [structure]
    reactant = "react.xyz"      # or .pdb, .gjf, .sdf
    product  = "prod.xyz"
    charge = 0                  # total charge
    mult = 1                    # multiplicity (2S+1)
    ''')


# ----------------------
# calculator blocks
# ----------------------
def calculator_uma():
    return textwrap.dedent('''\
    [calculator]
    type = "UMA"                # Universal Models for Atoms
    model_name = "uma-s-1p1"    # or "uma-m-1p1"
    task = "omol"               # omol task for molecules.
    device = "gpu"              # or "cpu"
    ''')


def calculator_uma_verbose():
    return textwrap.dedent('''\
    [calculator]
    type = "UMA"                    # Universal Models for Atoms
    model_name = "uma-s-1p1"        # or "uma-m-1p1"
    task = "omol"                   # omol task for molecules.
    device = "gpu"                  # or "cpu"
    inference_preset = "auto"       # "turbo", "opt", "dmf", etc. "auto" for automatic selection from calculation type.
    # model_path = "/path/to/uma-s-1p1.pt"  # optional, auto-download if omitted.
    ''')


def calculator_tblite():
    return textwrap.dedent('''\
    # --- Example TBLite Calculator (.toml) ---
    # To use TBLite as the calculator, replace the UMA settings with the following.
    # For details, see: https://tblite.readthedocs.io/en/latest/users/ase.html

    [calculator]
    type = "tblite"                   # Use the TBLite calculator
    method = "GFN2-xTB"               # Semiempirical Tight-Binding method
    # solvation = ["alpb", "water"]   # Optional: ["solvation method", "solvent name"]
    totalcores = 16                   # total CPU cores to use.
    totalmem = 60                     # in GB. total CPU memory to use.
    ''')


def calculator_gaussian():
    return textwrap.dedent('''\
    # --- Example Gaussian Calculator (.toml) ---
    # To use Gaussian as the calculator, please ensure that "g16" command is available in your PATH.
    # For details, see: https://ase-lib.org/ase/calculators/gaussian.html

    [calculator]
    type = "gaussian"       # Use Gaussian calculator
    method = "B3LYP"        # Functional. e.g. "B3LYP", "M062X", "wB97XD"
    basis = "6-31G(d)"      # Basis set. e.g. "6-31G(d)", "def2-SVP"
    scf = ["xqc", "maxconventionalcycle=80"]  # scf(args1, args2). Comment out if not required.
    extra = "nosymm"        # Other inputs. e.g. "nosymm SCRF=... EmpiricalDispersion=..." Comment out if not required.
    totalcores = 16         # total CPU cores to use.
    totalmem = 60           # in GB. total CPU memory to use.
    ''')


def calculator_orca():
    return textwrap.dedent('''\
    # --- Example ORCA Calculator (.toml) ---
    # For details, see: https://ase-lib.org/ase/calculators/orca.html#module-ase.calculators.orca

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
    ''')


# ----------------------
# calculation blocks
# ----------------------
def calculation_opt():
    return textwrap.dedent('''\
    [calculation]
    type = "opt"        # Geometry optimization
    fmax = 0.02         # maximum force for convergence criterion (eV angstrom^-1)
    maxstep = 1000      # maximum number of optimization steps
    method = "FIRE"     # or "BFGS", "LBFGS", "BFGSLINESEARCH"
    ''')


def calculation_scan():
    return textwrap.dedent('''\
    [calculation]
    type = "scan"       # Scan calculation
    fmax = 0.02         # maximum force for convergence criterion (eV angstrom^-1)
    maxstep = 1000      # maximum number of optimization steps
    method = "FIRE"     # or "BFGS", "LBFGS", "BFGSLINESEARCH"

    [scan]  # Start the scan section
    # 1D-scan: Set one [[scan.XXX]], 2D-scan: Set two [[scan.XXX]]
    [[scan.distance]]
    atomid = "1, 2"         # or [1, 2]. 1-based indices for target atoms.
    value = "1.5, 1.0, 6"   # or [1.5, 1.0, 6]. Start, End, Nsteps. in angstrom.

    [[scan.angle]]
    atomid = "3, 4, 5"      # or [3, 4, 5]. 1-based indices for target atoms.
    value = "120, 180, 5"   # or [120, 180, 5]. Start, End, Nsteps. in degrees.

    [[scan.dihedral]]
    atomid = "6, 7, 8, 9"   # or [6, 7, 8, 9]. 1-based indices for target atoms.
    value = "90, 180, 10"   # or [90, 180, 10]. start, stop, Nsteps. in degrees.
    ''')


def calculation_scan_verbose():
    return textwrap.dedent('''\
    [calculation]
    type = "scan"       # Scan calculation
    fmax = 0.02         # maximum force for convergence criterion (eV angstrom^-1)
    maxstep = 1000      # maximum number of optimization steps
    method = "FIRE"     # or "BFGS", "LBFGS", "BFGSLINESEARCH"
    peak_vibration = "highest"  # subsequent vibrational analysis for: "highest" energy peak, "all" peaks, or "none". default: "highest"

    [scan]  # Start the scan section
    # true: use previous optimized structure as next initial structure (Gaussian-like)
    # false: all initial structures are modeled from the input coordinates (independent)
    sequential = true  # or false

    # 1D-scan: Set one [[scan.XXX]], 2D-scan: Set two [[scan.XXX]]
    [[scan.distance]]
    atomid = "1, 2"         # or [1, 2]. 1-based indices for target atoms.
    value = "1.5, 1.0, 6"   # or [1.5, 1.0, 6]. Start, End, Nsteps. in angstrom.

    [[scan.angle]]
    atomid = "3, 4, 5"      # or [3, 4, 5]. 1-based indices for target atoms.
    value = "120, 180, 5"   # or [120, 180, 5]. Start, End, Nsteps. in degrees.

    [[scan.dihedral]]
    atomid = "6, 7, 8, 9"   # or [6, 7, 8, 9]. 1-based indices for target atoms.
    value = "90, 180, 10"   # or [90, 180, 10]. start, stop, Nsteps. in degrees.
    ''')


def calculation_dmf():
    return textwrap.dedent('''\
    [calculation]
    type = "dmf"                # Direct MaxFlux calculation
    dmf_nmove = 20              # nmove for path optimization. default: 5.
    update_teval = false        # or true. default: false
    dmf_convergence = "tight"   # or "middle", "loose". default: "tight"
    ''')


def calculation_dmf_verbose():
    return textwrap.dedent('''\
    [calculation]
    type = "dmf"                # Direct MaxFlux calculation

    # --- initial path ---
    fbenm = "cfbenm"            # or "fbenm". "cfbenm" for correlated = true. default: "cfbenm"
    fbenm_nmove = 10            # nmove for initial path estimation. default: 10

    # --- path optimization ---
    dmf_nmove = 20              # nmove for path optimization. default: 5.
    update_teval = false        # or true. default: false
    dmf_convergence = "tight"   # or "middle", "loose". default: "tight"
    parallel = false            # Enable parallel execution across nmoves. CPU only, not available for GPU. default: false
    peak_vibration = "highest"  # subsequent vibrational analysis for: "highest" energy peak, "all" peaks, or "none". default: "highest"
    ''')


def calculation_neb():
    return textwrap.dedent('''\
    [calculation]
    type = "neb"                # NEB calculation
    n_images = 5                # total images including endpoints.
    climb = true                # climbing image on/off. default: true
    idpp = true                 # use IDPP interpolation. default: true
    method = "FIRE"             # optimizer ("FIRE", "BFGS", "LBFGS", "BFGSLineSearch")
    fmax = 0.02                 # maximum force for convergence criterion (eV angstrom^-1)
    maxstep = 1000              # maximum number of optimization steps
    ''')


def calculation_neb_verbose():
    return textwrap.dedent('''\
    [calculation]
    type = "neb"                # NEB calculation
    neb_method = "aseneb"       # or "improvedtangent",  "eb", "spline", "string". https://ase-lib.org/ase/neb.html for details.
    n_images = 5                # total images including endpoints.
    climb = true                # climbing image on/off. default: true
    idpp = true                 # use IDPP interpolation. default: true
    method = "FIRE"             # optimizer ("FIRE", "BFGS", "LBFGS", "BFGSLineSearch")
    fmax = 0.01                 # maximum force for convergence criterion (eV angstrom^-1)
    maxstep = 1000              # maximum number of optimization steps
    parallel = false            # Enable parallel execution across images. CPU only, not available for GPU. default: false
    peak_vibration = "highest"  # subsequent vibrational analysis for: "highest" energy peak, "all" peaks, or "none". default: "highest"
    k = 0.1                     # (optional) NEB spring constant (eV angstrom^-2). default: 0.1
    ''')


def calculation_ts():
    return textwrap.dedent('''\
    [calculation]
    type = "ts"                 # Transition State Search using Sella.
    fmax = 0.001                # maximum force for convergence criterion (eV angstrom^-1)
    maxstep = 1000              # maximum number of optimization steps
    internal = true             # or true. Use internal coordinates. default: true

    [[calculation.irc]]         # Subsequent IRC calculation. Remove this block to skip.
    irc_fmax = 0.05             # maximum force for convergence criterion (eV angstrom^-1). default: 0.05
    irc_dx = 0.1                # in angstrom sqrt(amu). step size for IRC. default: 0.1
    irc_steps = 100             # maximum number of IRC steps. default: 100
    endpoint_opt = "LBFGS"      # or "FIRE", "LBFGS", "BFGSLINESEARCH". Optimize the IRC endpoints.
                                # fmax and maxstep settings are taken from the TS opt. false: skip endpoint optimization.
    ''')


def calculation_ts_verbose():
    return textwrap.dedent('''\
    [calculation]
    type = "ts"                 # Transition State Search using Sella.
    fmax = 0.001                # maximum force for convergence criterion (eV angstrom^-1)
    maxstep = 1000              # maximum number of optimization steps
    internal = true             # or true. Use internal coordinates. default: true
    initial_hess = "vib"        # initial hessian. "vib": ASE Vibrations. "sella": Sella's default estimate. "autograd": UMA only, but not available yet.
                                # or "/path/to/hessian.csv | .npy" hessian_2d() format.
    update_hess = "vib"         # hessian update method. "vib": ASE Vibrations. "sella": Sella's default update. "autograd": UMA only, but not available yet.

    [[calculation.irc]]         # Subsequent IRC calculation. Remove this block to skip.
    irc_fmax = 0.05             # maximum force for convergence criterion (eV angstrom^-1). default: 0.05
    irc_dx = 0.1                # in angstrom sqrt(amu). step size for IRC. default: 0.1
    irc_steps = 100             # maximum number of IRC steps. default: 100
    endpoint_opt = "LBFGS"      # or "FIRE", "LBFGS", "BFGSLINESEARCH". Optimize the IRC endpoints.
                                # fmax and maxstep settings are taken from the TS opt. false: skip endpoint optimization.
    ''')


def calculation_minhop():
    return textwrap.dedent('''\
    # For details, see: https://ase-lib.org/ase/optimize.html#minima-hopping
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
    ''')    


def calculation_nvemd():
    return textwrap.dedent('''\
    # For details, see: https://ase-lib.org/ase/md.html
    [calculation]
    type = "nvemd"      # NVE Molecular Dynamics. Velocity Verlet algorithm
    timestep = 1.0      # in fs, timestep for MD simulations. default: 1.0 fs
    nsteps = 1000       # number of steps for MD simulation. e.g. 1000 * 1 fs = 1 ps 
    gen_temp = 300      # in K, temperature for Maxwell-Boltzmann distribution. default: 300 K
    nstlog = 10         # update log every nstlog steps (log interval).
    ''')


def constraints():
    return textwrap.dedent('''\
    # --- Example Constraints (.toml) ---

    [constraints]  # Start the constraints section

    # apply a positional constraint (FixAtoms in ASE)
    [[constraints.freeze]]
    atomid = "1, 2, 3"      # or [1, 2, 3]. 1-based indices for target atoms.

    # apply a distance constraint (FixInternals in ASE)
    [[constraints.distance]]
    atomid = "4, 5"         # or [4, 5]. 1-based indices for target atoms.
    value = 1.5             # in angstrom

    # apply an angle constraint (FixInternals in ASE)
    [[constraints.angle]]
    atomid = "6, 7, 8"      # or [6, 7, 8]. 1-based indices for target atoms.
    value = 120             # in degrees

    # apply a dihedral constraint (FixInternals in ASE)
    [[constraints.dihedral]]
    atomid = "9, 10, 11, 12"  # or [9, 10, 11, 12]. 1-based indices for target atoms.
    value = 180               # in degrees
    ''')


def constraints_verbose():
    return textwrap.dedent('''\
    # --- Example Constraints (.toml) ---
    # constraints can be applied to the Optimization, Scan, and Minima-Hopping calculations

    [constraints]  # Start the constraints section

    # apply a positional constraint (FixAtoms in ASE)
    [[constraints.freeze]]
    atomid = "1, 2, 3"      # or [1, 2, 3]. 1-based indices for target atoms.

    # apply a distance constraint (FixInternals in ASE)
    [[constraints.distance]]
    atomid = "4, 5"         # or [4, 5]. 1-based indices for target atoms.
    value = 1.5             # in angstrom

    # apply an angle constraint (FixInternals in ASE)
    [[constraints.angle]]
    atomid = "6, 7, 8"      # or [6, 7, 8]. 1-based indices for target atoms.
    value = 120             # in degrees

    # apply a dihedral constraint (FixInternals in ASE)
    [[constraints.dihedral]]
    atomid = "9, 10, 11, 12"  # or [9, 10, 11, 12]. 1-based indices for target atoms.
    value = 180               # in degrees

    # apply a Hookean restraint (Hookean in ASE)
    [[constraints.hookean]]
    atomid = "13, 14"       # or [13, 14]. 1-based indices for target atoms.
    k = 10                  # Hooke's law (spring) constant (kcal mol^-1 angstrom^-2)
    r0 = 1.5                # in angstrom. When distance > r0, apply harmonic restraint.
    ''')

# ----------------------
# main function
# ----------------------
def example_list_for_cli():
    return [
        "opt",
        "scan",
        "dmf",
        "neb",
        "ts",
        "minhop",
        "nvemd",
        "constraints",
        "uma",
        "tblite",
        "gaussian",
        "orca"
    ]


def make_examples(name: str, verbose=False, to_file=False) -> str | None:
    if name == "opt":
        header = "# --- Example OPT input (.toml) ---"
        content = "\n".join([
            header,
            structure_single(),
            calculator_uma_verbose() if verbose else calculator_uma(),
            calculation_opt()
        ])
    
    elif name == "scan":
        header = "# --- Example Scan input (.toml) ---"
        content = "\n".join([
            header,
            structure_single(),
            calculator_uma_verbose() if verbose else calculator_uma(),
            calculation_scan_verbose() if verbose else calculation_scan()
        ])

    elif name == "dmf":
        header = "# --- Example DMF input (.toml) ---"
        content = "\n".join([
            header,
            structure_double(),
            calculator_uma_verbose() if verbose else calculator_uma(),
            calculation_dmf_verbose() if verbose else calculation_dmf()
        ])

    elif name == "neb":
        header = "# --- Example NEB input (.toml) ---"
        content = "\n".join([
            header,
            structure_double(),
            calculator_uma_verbose() if verbose else calculator_uma(),
            calculation_neb_verbose() if verbose else calculation_neb()
        ])

    elif name == "minhop":
        header = "# --- Example Minima-Hopping input (.toml) ---"
        content = "\n".join([
            header,
            structure_single(),
            calculator_uma_verbose() if verbose else calculator_uma(),
            calculation_minhop()
        ])

    elif name == "nvemd":
        header = "# --- Example NVE MD input (.toml) ---"
        content = "\n".join([
            header,
            structure_single(),
            calculator_uma_verbose() if verbose else calculator_uma(),
            calculation_nvemd() if verbose else calculation_nvemd()
        ])

    elif name == "ts":
        header = "# --- Example TS input (.toml) ---"
        content = "\n".join([
            header,
            structure_single(),
            calculator_uma_verbose() if verbose else calculator_uma(),
            calculation_ts_verbose() if verbose else calculation_ts()
        ])

    elif name == "constraints":
        content = constraints_verbose() if verbose else constraints()

    elif name == "uma":
        content = calculator_uma_verbose() if verbose else calculator_uma()

    elif name == "tblite":
        content = calculator_tblite()
    
    elif name == "gaussian":
        content = calculator_gaussian()

    elif name == "orca":
        content = calculator_orca()

    else:
        raise ValueError(f"Unknown example name: {name}. Select from {example_list_for_cli()}")
    
    # --- write to file or return string ---
    if to_file:
        filename = f"template_{name}.toml"
        lines = content.strip().splitlines()
        new_content = "\n".join(lines[1:]) + "\n\n"
        from pathlib import Path
        Path(filename).write_text(new_content)
        print(f"Example TOML written to: {filename}")
        return None
    else:
        return content 
