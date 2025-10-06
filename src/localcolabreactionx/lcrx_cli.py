# ---- general imports -----
from pathlib import Path
import tomllib
import time
import warnings

# --- third-party imports ---
import typer
from loguru import logger

# ---- LocalColabReaction imports ----
from localcolabreactionx.utils.environ_utils import set_environment_variable
from localcolabreactionx.utils.logging_utils import setup_logger, log_launch_info
from localcolabreactionx.utils.toml_examples import make_examples, example_list_for_cli

# to hide future warnings
warnings.simplefilter('ignore', FutureWarning)  

# application instance
app = typer.Typer(help="LocalColabReactionX: ", no_args_is_help=True)


# --- apps: run calculations ---
@app.command(no_args_is_help=True)
def run(
    toml_path: Path = typer.Option(..., "--toml", "-t", "-i", help="Path to input TOML file."),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable DEBUG logging."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress INFO logs.")
):
    """
    Run a calculation from TOML input.
    Example: localcolabreaction run --toml input.toml
    """
    t0 = time.perf_counter()  # start timer

    setup_logger(verbose, quiet)
    log_launch_info()
    logger.info(f"Reading input: {toml_path}")

    try:
        with open(toml_path, "rb") as f:
            data = tomllib.load(f)
    except Exception as e:
        logger.error(f"Failed to parse TOML: {e}")
        raise typer.Exit(code=1)

    # set environment variables
    set_environment_variable(data)

    calc_type = data.get("calculation", {}).get("type", "").lower()
    logger.info(f"Calculation type detected: {calc_type}")

    if calc_type == "opt":
        from localcolabreactionx.runners.opt_runner import run_opt
        logger.info("Running Structure Optimization...")
        run_opt(toml_path)

    elif calc_type == "dmf":
        from localcolabreactionx.runners.dmf_runner import run_dmf
        logger.info("Running Direct MaxFlux calculation...")
        run_dmf(toml_path)

    elif calc_type == "scan":
        from localcolabreactionx.runners.scan_runner import run_scan
        logger.info("Running Relaxed Scanning...")
        run_scan(toml_path)

    elif calc_type == "minhop":
        from localcolabreactionx.runners.minhop_runner import run_minhop
        logger.info("Running Minima-Hopping simulation...")
        run_minhop(toml_path)
    elif calc_type == "nvemd":
        from localcolabreactionx.runners.nvemd_runner import run_nvemd
        logger.info("Running NVE Molecular Dynamics simulation...")
        run_nvemd(toml_path)
    elif calc_type == "ts":
        from localcolabreactionx.runners.ts_runner import run_ts
        logger.info("Running Transition State Search...")
        run_ts(toml_path)
    elif calc_type == "neb":
        from localcolabreactionx.runners.neb_runner import run_neb
        logger.info("Running Nudged Elastic Band calculation...")
        run_neb(toml_path)
    else:
        logger.error(f"Unknown calculation type: {calc_type}")
        raise typer.Exit(code=1)

    t1 = time.perf_counter()  # end timer

    elapsed = t1 - t0
    logger.info(f"Calculation finished in {elapsed:.2f} seconds.")

    # citations
    from localcolabreactionx.utils.citation_utils import collect_and_log_citations
    collect_and_log_citations(data, verbose=verbose)

    logger.info("LocalColabReaction: Normal termination.")


# --- apps: output example .toml ---
@app.command(no_args_is_help=True)
def template(
    name: str = typer.Argument(..., help=f"Name to show template TOML. {example_list_for_cli()}"),
    to_file: bool = typer.Option(False, "--out", "-o", help="Write to file instead of stdout."),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose output.")
):
    """
    Show or save template TOML inputs. Example: localcolabreaction template opt
    """
    name = name.lower()
    try:
        content = make_examples(name, verbose=verbose, to_file=to_file)
        typer.echo(content)
    except KeyError:
        logger.error(f"Unknown template name: {name}")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()


