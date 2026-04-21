import os
from loguru import logger 


def set_environment_variable(toml_data) -> None:
    calculator_data = toml_data.get("calculator", {})
    calculation_data = toml_data.get("calculation", {})

    # get totalcores and totalmem
    totalcores = calculator_data.get("totalcores", None)
    totalmem = calculator_data.get("totalmem", None)

    # --- TBLite --- 
    # For details: https://tblite.readthedocs.io/en/latest/tutorial/parallel.html
    if calculator_data.get("type", "").lower() == "tblite":
        logger.info("TBLite selected: set environment variables.")
        os.environ["OPENBLAS_NUM_THREADS"] = "1"
        os.environ["MKL_NUM_THREADS"] = "1"
        os.environ["BLIS_NUM_THREADS"] = "1"
        os.environ["NUMEXPR_NUM_THREADS"] = "1"
        os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"

        if totalcores is not None:
            os.environ["OMP_NUM_THREADS"] = str(totalcores)
        
        if totalcores is not None and totalmem is not None:
            mem_per_core = int(totalmem) // int(totalcores)
            if mem_per_core >= 8:
                stack = "3G"
            elif mem_per_core >= 4:
                stack = "2G"
            elif mem_per_core >= 2:
                stack = "1G"
            else:                   
                stack = "512M"
            
            os.environ["OMP_STACKSIZE"] = stack

    # logger.debug(f"Environment variables set: {os.environ}")

    # --- jax ---
    if calculation_data.get("type", "").lower() == "ts":
        calc_device = calculation_data.get("device", "cpu").lower()

        if calc_device in {"gpu", "cuda", "GPU"}:
            os.environ["JAX_PLATFORMS"] = "cuda"
            os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
        else:
            os.environ["JAX_PLATFORMS"] = "cpu"

        import jax

        backend = jax.default_backend()
        devices = jax.devices()

        logger.info(f"TS calculation requested device = {calc_device}")
        logger.info(f"JAX_PLATFORMS = {os.environ.get('JAX_PLATFORMS')}")
        logger.info(f"JAX default backend = {backend}")
        logger.info(f"JAX visible devices = {devices}")

        if calc_device in {"gpu", "cuda"}:
            if backend == "gpu":
                logger.info("JAX GPU backend is active.")
            else:
                logger.warning("GPU was requested, but JAX is not using GPU backend.")
        else:
            if backend == "cpu":
                logger.info("JAX CPU backend is active.")
            else:
                logger.warning("CPU was requested, but JAX is not using CPU backend.")
