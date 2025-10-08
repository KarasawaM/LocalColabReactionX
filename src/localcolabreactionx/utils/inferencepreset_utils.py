from fairchem.core.units.mlip_unit.api.inference import InferenceSettings, inference_settings_default, inference_settings_turbo
from loguru import logger

"""
https://github.com/facebookresearch/fairchem/blob/8fe865a1/docs/core/common_tasks/ase_calculator.md
https://github.com/facebookresearch/fairchem/blob/main/src/fairchem/core/units/mlip_unit/api/inference.py#L92

default:
    return InferenceSettings(
        tf32=False,
        activation_checkpointing=True,
        merge_mole=False,
        compile=False,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )

turbo: compile takes more time, but speeds up the inference.
    return InferenceSettings(
        tf32=True,
        activation_checkpointing=True,
        merge_mole=True,
        compile=True,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )
"""


def inference_settings_dmf():
    return InferenceSettings(
        tf32=False,
        activation_checkpointing=False,
        merge_mole=True,
        compile=False,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )


def inference_settings_opt_scan():
    return InferenceSettings(
        tf32=True,
        activation_checkpointing=False,
        merge_mole=True,
        compile=False,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )


def inference_settings_ts():
    return InferenceSettings(
        tf32=False,
        activation_checkpointing=False,
        merge_mole=True,
        compile=False,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )


def inference_settings_minhop():
    return InferenceSettings(
        tf32=True,
        activation_checkpointing=False,
        merge_mole=True,
        compile=True,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )


def inference_settings_md():
    return InferenceSettings(
        tf32=True,
        activation_checkpointing=False,
        merge_mole=True,
        compile=True,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )


def inference_settings_md_large():
    return InferenceSettings(
        tf32=True,
        activation_checkpointing=True,
        merge_mole=True,
        compile=True,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )


def inference_settings_custom():
    return InferenceSettings(
        tf32=False,
        activation_checkpointing=False,
        merge_mole=False,
        compile=False,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    )


# --- load inference settings ---
def load_inference_preset(preset: str, calculation_type: str, calc_data: str) -> InferenceSettings:
    preset = preset.lower()
    preset_dict = {
        "default": inference_settings_default(),
        "turbo": inference_settings_turbo(),
        "dmf": inference_settings_dmf(),
        "neb": inference_settings_dmf(),
        "opt": inference_settings_opt_scan(),
        "scan": inference_settings_opt_scan(),
        "minhop": inference_settings_minhop(),
        "md": inference_settings_md(),
        "md_large": inference_settings_md_large(),
        "ts": inference_settings_ts(),
    }

    # auto detection
    if preset == "auto":
        # md case
        if calculation_type.lower() in ["md", "nvemd", "nvtmd", "nptmd"]:
            preset = "md"
        else:
            preset = calculation_type.lower()

    # custom case
    if preset == "custom":
        logger.info("Using custom inference settings from [calculation] block.")
        return InferenceSettings(
            tf32=bool(calc_data.get("tf32", False)),
            activation_checkpointing=bool(calc_data.get("activation_checkpointing", False)),
            merge_mole=bool(calc_data.get("merge_mole", False)),
            compile=bool(calc_data.get("compile", False)),
            wigner_cuda=bool(calc_data.get("wigner_cuda", False)),
            external_graph_gen=bool(calc_data.get("external_graph_gen", False)),
            internal_graph_gen_version=int(calc_data.get("internal_graph_gen_version", 2)),
        )

    # search preset
    if preset not in preset_dict:
        raise ValueError(f"Unknown preset: {preset}. Select from {list(preset_dict.keys())}")

    logger.info(f"Using inference preset: {preset}. {preset_dict[preset]}")
    return preset_dict.get(preset, None)

