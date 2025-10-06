from __future__ import annotations

# ---- general imports ----
import tomllib
from pathlib import Path

# ---- third-party imports ----
from loguru import logger
from ase import Atoms
from ase.io import read


# ---- LocalColabReactionX imports ----
from .constraintbuilder import ConstraintBuilder
from ..utils.logging_utils import setup_logger, log_launch_info, log_atoms_info, log_constraints, log_toml_data


class BaseASEBuilder:
    def __init__(self, toml_path: str | Path):
        self.toml_path = Path(toml_path)
        self.data = self._load_toml()

        try:
            self.struct_data = self.data["structure"]
            self.calc_data = self.data["calculator"]
        except KeyError as e:
            raise ValueError(f"Missing '{e.args[0]}' section in TOML data.") from None

        self.scan_data = self.data["scan"] if "scan" in self.data else {}
        self._calc_factory = self._build_calc_factory()  # initialize calculator factory

    def _load_toml(self) -> dict:
        with open(self.toml_path, "rb") as f:
            data = tomllib.load(f)
            log_toml_data(data, label="TOML loaded as follows")
        
        return data

    def _build_calc_factory(self):
        calc_type = self.calc_data["type"].lower()
        # Get total memory and cores for parallelization
        totalmem = self.calc_data.get("totalmem", None)  # in GB. total CPU memory.
        totalcores = self.calc_data.get("totalcores", None)  # total CPU cores.

        # UMA
        if calc_type in {"uma", "uma-s-1p1", "uma-m-1p1"}:
            from fairchem.core import pretrained_mlip, FAIRChemCalculator
            from fairchem.core.units.mlip_unit import load_predict_unit
            from ..utils.inferencepreset_utils import load_inference_preset

            model_name = self.calc_data.get("model_name", None)
            device_key = self.calc_data.get("device", "gpu").lower()
            device = "cuda" if device_key in {"gpu", "cuda"} else "cpu"
            task = self.calc_data.get("task", "omol")
            model_path = self.calc_data.get("model_path", None)
            parallel = self.data.get("calculation", {}).get("parallel", False)
            calculation_type = self.data.get("calculation", {}).get("type", None)  
            inference_preset = self.calc_data.get("inference_preset", "auto")

            # inference settings
            inference_settings = load_inference_preset(inference_preset, calculation_type=calculation_type, calc_data=self.calc_data)

            # parallel gpu run not supported yet.
            if parallel and device == "gpu":
                logger.warning("Parallel execution with GPU is not supported. Switching to serial execution.")
                parallel = False

            # serial run: load predictor once and reuse this.
            if parallel is False:
                if model_path:
                    self._predictor = load_predict_unit(model_path, device=device, inference_settings=inference_settings)
                elif model_name is not None:
                    self._predictor = pretrained_mlip.get_predict_unit(model_name, device=device, inference_settings=inference_settings)
                else:
                    raise ValueError("Either 'model_name' or 'model_path' must be specified in the calculator section.")

            def factory(**overrides):
                if parallel and device == "cpu":
                    # set CPU settings for torch
                    try:
                        import torch
                        torch.set_num_threads(int(self.calc_data.get("torch_threads", 1)))
                        torch.set_num_interop_threads(int(self.calc_data.get("torch_interop_threads", 1)))
                    except Exception:
                        pass

                    # load predictor for each Calculator instance.
                    if model_path:
                        predictor = load_predict_unit(model_path, device=device, inference_settings=inference_settings)
                    else:
                        predictor = pretrained_mlip.get_predict_unit(model_name, device=device, inference_settings=inference_settings)
                else:
                    predictor = self._predictor  # serial run

                return FAIRChemCalculator(predictor, task_name=task)            
            return factory

        # TBLite
        elif calc_type in {"tblite", "xtb"}:
            from tblite.ase import TBLite
            method = self.calc_data.get("method", "GFN2-xTB")
            charge = self.struct_data["charge"]
            multiplicity = self.struct_data["mult"]
            solvation = self.calc_data.get("solvation", None)

            def factory(**overrides):
                # overrides currently not used.
                return TBLite(method=method, charge=charge, 
                              multiplicity=multiplicity, 
                              solvation=solvation, verbosity=0
                              )
            return factory

        # EMT. Only for test.
        elif calc_type == "emt":
            from ase.calculators.emt import EMT
            
            def factory(**overrides):
                # overrides currently not used.
                return EMT()
            return factory

        # Gaussian
        elif calc_type == "gaussian":
            from ase.calculators.gaussian import Gaussian
            charge = int(self.struct_data.get("charge", 0))
            mult = int(self.struct_data.get("mult", 1))
            method = self.calc_data.get("method", "PM7")
            basis = self.calc_data.get("basis", None)
            scf = self.calc_data.get("scf", None)
            extra = self.calc_data.get("extra", None)
            output_type = self.calc_data.get("output_type", "p")

            def factory(**overrides):
                outpath = overrides.get("outpath", "g16/g16")  # output prefix of g16 results.
                mem = overrides.get("mem", f"{totalmem}GB")  # in GB. set as %mem in .gjf file
                nprocshared = overrides.get("nprocshared", str(totalcores))  # set as %nprocshared in .gjf file
                chk = outpath.split("/")[-1] + ".chk" # g16.chk as default
  
                kwargs = dict(
                    charge=charge, mult=mult, chk=chk,
                    method=method, basis=basis,
                    scf=scf, extra=extra,
                    mem=mem, nprocshared=nprocshared,
                    output_type=output_type,
                    label=outpath,
                )

                kwargs = {k: v for k, v in kwargs.items() if v is not None}

                return Gaussian(**kwargs)
            return factory

        # ORCA
        elif calc_type == "orca":
            from ase.calculators.orca import ORCA
            from ase.calculators.orca import OrcaProfile
            charge = int(self.struct_data.get("charge", 0))
            mult = int(self.struct_data.get("mult", 1))
            orcapath = self.calc_data.get("orcapath", None)
            orcasimpleinput = self.calc_data.get("orcasimpleinput", None)
            orcablocks_in = self.calc_data.get("orcablocks", "")

            # Maxcore: in MB / nprocs
            def _make_orca_block(orcablocks, totalcores, totalmem):
                orcablocks += f"%pal nprocs {totalcores} end\n"
                maxcore = (totalmem * 10 // totalcores) * 100  # 60 GB, 16procs -> maxcore 3700 MB.
                orcablocks += f"%maxcore {maxcore}\n"
                return orcablocks

            orcablocks = _make_orca_block(orcablocks_in, totalcores, totalmem) if (totalcores and totalmem) else orcablocks_in

            def factory(**overrides):
                orcablock = overrides.get("orcablocks", orcablocks)  # orcablocks: additional orca input blocks.
                outpath = overrides.get("outpath", "orca")  # orcaout: output directory for orca results.
                profile = OrcaProfile(command=orcapath)
                return ORCA(
                    profile=profile,
                    charge=charge, mult=mult,
                    orcasimpleinput=orcasimpleinput, orcablocks=orcablock,
                    directory=outpath
                )
            return factory

        else:
            raise ValueError(f"Unsupported calculator type: {calc_type}")

    # API: create a calculator factory with the loaded parameters
    def make_calculator(self, **overrides):
        return self._calc_factory(**overrides)

    # Set the calculator to the Atoms object
    def _set_calculator(self, atoms: Atoms, **overrides) -> None:
        atoms.calc = self.make_calculator(**overrides)

    def _apply_constraints(self, atoms: Atoms):
        if "constraints" in self.data:
            constraint_builder = ConstraintBuilder(self.data["constraints"])
            constraint_builder.apply(atoms)


class SingleASEBuilder(BaseASEBuilder):
    def build(self) -> Atoms:
        atoms = read(self.struct_data["input"])
        atoms.info["charge"] = self.struct_data["charge"]
        atoms.info["spin"] = self.struct_data["mult"]
        self._apply_constraints(atoms)
        self._set_calculator(atoms)

        # logging
        logger.info(f"{self.__class__.__name__}: Single structure loaded.")
        log_atoms_info(atoms, label="Input")
        log_constraints(atoms, label="Input")

        return atoms


class DoubleASEBuilder(BaseASEBuilder):
    def build(self) -> tuple[Atoms, Atoms]:
        react = read(self.struct_data["reactant"])
        prod = read(self.struct_data["product"])
        for atoms in (react, prod):
            atoms.info["charge"] = self.struct_data["charge"]
            atoms.info["spin"] = self.struct_data["mult"]
            self._apply_constraints(atoms)
            self._set_calculator(atoms)

        # logging
        logger.info(f"{self.__class__.__name__}: Two structures loaded.")
        log_atoms_info(react, label="Reactant")
        log_constraints(react, label="Reactant")
        log_atoms_info(prod, label="Product")
        log_constraints(prod, label="Product")

        # checking reactant and product
        self._check_reactant_product(react, prod)

        return react, prod

    def _check_reactant_product(self, react: Atoms, prod: Atoms):
        # atom count check
        if len(react) != len(prod):
            logger.error(f"Atom count mismatch detected.")
            raise ValueError("Reactant and product must have the same number of atoms.")

        # atom ordering check
        if react.get_chemical_symbols() != prod.get_chemical_symbols():
            logger.error("Atom ordering mismatch detected.")
            raise ValueError("Reactant and product must have the same atom ordering.")
        
        logger.info("Atom ordering are consistent. Proceeding...")


def select_builder(toml_path: str) -> SingleASEBuilder | DoubleASEBuilder:
    with open(toml_path, "rb") as f:
        data = tomllib.load(f)
    calc_type = data["calculation"]["type"].upper()

    if calc_type in {"NEB", "DMF"}:
        builder = DoubleASEBuilder(toml_path)
    else:
        builder = SingleASEBuilder(toml_path)

    return builder


if __name__ == "__main__":
    # --- Debugging / Example usage ---
    # run with a TOML file path:  python3 toml2ase.py input.toml
    import sys

    toml_path = sys.argv[1]

    setup_logger(verbose=True)  # init log
    log_launch_info()  # log environment info
    builder = select_builder(toml_path)

    if isinstance(builder, SingleASEBuilder):
        atoms = builder.build()
    else:
        react, prod = builder.build()

