from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any

import numpy as np
import torch
from torch.autograd import grad

from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.vibrations.data import VibrationsData

from fairchem.core import FAIRChemCalculator
from fairchem.core.datasets.atomic_data import AtomicData, atomicdata_list_to_batch


@dataclass
class AutogradHessianResult:
    energy: float
    forces: np.ndarray
    hessian: np.ndarray


class FAIRChemAutogradHessian:
    def __init__(
        self,
        predictor: Any,
        task_name: str,
        *,
        r_data_keys: tuple[str, ...] = ("spin", "charge"),
        molecule_cell_size: float | None = None,
        radius: float = 6.0,
        max_neigh: int | None = None,
    ) -> None:
        self.predictor = predictor
        self.task_name = task_name
        self.r_data_keys = tuple(r_data_keys)
        self.molecule_cell_size = molecule_cell_size
        self.radius = radius
        self.max_neigh = max_neigh

    def get_hessian(
        self,
        atoms: Atoms,
        *,
        vmap: bool = False,
        return_energy_forces: bool = False,
    ) -> np.ndarray | AutogradHessianResult:
        batch, pos = self._build_batch_with_grad(atoms)

        with self._temporarily_energy_only_mode():
            with torch.enable_grad():
                pred = self.predictor.predict(batch)

                if "energy" not in pred:
                    raise ValueError(
                        "predictor.predict(batch) did not return 'energy'. "
                        "Autograd Hessian requires an energy-predicting model."
                    )

                energy_tensor = pred["energy"].sum()

                # IMPORTANT: retain_graph=True because Hessian needs second derivatives
                forces_tensor = -grad(
                    energy_tensor,
                    pos,
                    create_graph=True,
                    retain_graph=True,
                )[0]

                hessian_np = self._assemble_hessian(
                    forces=forces_tensor,
                    positions=pos,
                    vmap=vmap,
                )

                energy_value = float(energy_tensor.detach().cpu().item())
                forces_np = forces_tensor.detach().cpu().numpy()

                del pred
                del energy_tensor
                del forces_tensor

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        if return_energy_forces:
            return AutogradHessianResult(
                energy=energy_value,
                forces=forces_np,
                hessian=hessian_np,
            )
        return hessian_np

    def get_vibrations_data(
        self,
        atoms: Atoms,
        *,
        vmap: bool = False,
    ) -> VibrationsData:
        H = self.get_hessian(atoms, vmap=vmap, return_energy_forces=False)
        return VibrationsData.from_2d(atoms, H)

    def attach_calculator(self, atoms: Atoms) -> Calculator:
        calc = FAIRChemCalculator(self.predictor, task_name=self.task_name)
        atoms.calc = calc
        return calc

    def _build_batch_with_grad(self, atoms: Atoms) -> tuple[Any, torch.Tensor]:
        target_dtype = getattr(
            self.predictor.inference_settings,
            "base_precision_dtype",
            torch.float32,
        )

        external_graph_gen = getattr(
            self.predictor.inference_settings,
            "external_graph_gen",
            False,
        )

        kwargs: dict[str, Any] = dict(
            task_name=self.task_name,
            r_data_keys=list(self.r_data_keys),
            target_dtype=target_dtype,
        )

        if self.molecule_cell_size is not None:
            kwargs["molecule_cell_size"] = self.molecule_cell_size

        if external_graph_gen:
            kwargs["r_edges"] = True
            kwargs["radius"] = self.radius
            kwargs["max_neigh"] = self.max_neigh if self.max_neigh is not None else 300

        data = AtomicData.from_ase(atoms, **kwargs)
        batch = atomicdata_list_to_batch([data])

        # fresh leaf tensor for positions
        pos = batch.pos.detach().clone().requires_grad_(True)
        batch.pos = pos
        return batch, pos

    @contextmanager
    def _temporarily_energy_only_mode(self):
        """
        Temporarily restrict predictor/model to energy-only tasks so that
        predictor.predict(batch) does NOT internally consume the autograd graph
        by computing existing force/stress tasks.
        """
        model_module = self.predictor.model.module

        old_tasks = model_module.tasks.copy()
        old_dataset_to_tasks = {
            k: list(v) for k, v in model_module.dataset_to_tasks.items()
        }

        # Save regress_config states if present
        saved_regress = []
        for head in model_module.output_heads.values():
            actual_head = head.head if hasattr(head, "head") else head
            if hasattr(actual_head, "regress_config"):
                rc = actual_head.regress_config
                saved_regress.append(
                    (
                        actual_head,
                        getattr(rc, "forces", None),
                        getattr(rc, "stress", None),
                        getattr(rc, "hessian", None),
                        getattr(rc, "hessian_vmap", None),
                    )
                )

        try:
            # Keep only energy tasks
            energy_tasks = {
                name: task
                for name, task in old_tasks.items()
                if getattr(task, "property", None) == "energy"
            }
            if not energy_tasks:
                raise RuntimeError("No energy task found; cannot compute autograd Hessian.")

            model_module._tasks = energy_tasks

            # rebuild dataset_to_tasks
            new_dataset_to_tasks = {}
            for task in energy_tasks.values():
                for ds in task.datasets:
                    new_dataset_to_tasks.setdefault(ds, []).append(task)
            model_module._dataset_to_tasks = new_dataset_to_tasks

            # Disable derivative regression in heads
            for actual_head, _, _, _, _ in saved_regress:
                rc = actual_head.regress_config
                if hasattr(rc, "forces"):
                    rc.forces = False
                if hasattr(rc, "stress"):
                    rc.stress = False
                if hasattr(rc, "hessian"):
                    rc.hessian = False

            yield

        finally:
            model_module._tasks = old_tasks
            model_module._dataset_to_tasks = old_dataset_to_tasks

            for actual_head, old_forces, old_stress, old_hessian, old_hessian_vmap in saved_regress:
                rc = actual_head.regress_config
                if hasattr(rc, "forces") and old_forces is not None:
                    rc.forces = old_forces
                if hasattr(rc, "stress") and old_stress is not None:
                    rc.stress = old_stress
                if hasattr(rc, "hessian") and old_hessian is not None:
                    rc.hessian = old_hessian
                if hasattr(rc, "hessian_vmap") and old_hessian_vmap is not None:
                    rc.hessian_vmap = old_hessian_vmap

    def _assemble_hessian(
        self,
        *,
        forces: torch.Tensor,
        positions: torch.Tensor,
        vmap: bool,
    ) -> np.ndarray:
        forces_flat = forces.reshape(-1)

        if vmap:
            eye = torch.eye(
                forces_flat.numel(),
                device=forces_flat.device,
                dtype=forces_flat.dtype,
            )

            def _one_component(vec: torch.Tensor) -> torch.Tensor:
                g = grad(
                    -forces_flat,
                    positions,
                    grad_outputs=vec,
                    retain_graph=True,
                    create_graph=False,
                )[0]
                return g.reshape(-1)

            rows = torch.vmap(_one_component)(eye)
            hessian = rows.detach().cpu().numpy()
        else:
            n = forces_flat.numel()
            hessian = np.empty((n, n), dtype=np.float64)
            for i in range(n):
                g = grad(
                    -forces_flat[i],
                    positions,
                    retain_graph=(i < n - 1),
                    create_graph=False,
                )[0]
                hessian[:, i] = g.reshape(-1).detach().cpu().numpy()

        # enforce symmetry slightly
        hessian = 0.5 * (hessian + hessian.T)
        return hessian
    
