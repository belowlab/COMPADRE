#!/usr/bin/env python3.14
# Runner for the test to check the memoy and runtime of the perl load
# data function and the python server load_data function
import json  # We are going to use the configuration file
from random import sample
import subprocess
import sys
from pathlib import Path

import pytest

# WE need to get the project root so that we can load the correct functions
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

# Load in the configuration object
CONFIG_PATH = project_root / "benchmark_config.json"
if not CONFIG_PATH.exists():
    raise FileNotFoundError(
        f"Unable to find the configuration json file in the project directory: {project_root}. Please make sure there is a file called benchmark_config.json."
    )

with open(CONFIG_PATH, "r") as config_fh:
    CONFIG = json.load(config_fh)


def run_benchmarks():
    """Main runner script that will run both the perl benchmark and the python benchmark"""
    # Make sure the output directory exist
    synthetic_data_dir = Path(CONFIG.get("output_dir"))

    if not synthetic_data_dir.exists():
        synthetic_data_dir.mkdir(exist_ok=True)

    for config_opts in CONFIG.get("params"):
        sample_size = config_opts["samples"]
        # We will in time need to add the pairs argument here
        related_fraction = config_opts["related_fraction"]

        # Now we need to run the script to generate the synthetic data
        gen_cmd = [
            sys.executable,
            "scratch/generate_synthetic_data.py",
            "--samples",
            f"{N}",
            "--pairs",
            f"{pairs}",
            "--related-fraction",
            "0.95",
            "--output-dir",
            f"{synthetic_data_dir}",
        ]

        subprocess.run(gen_cmd, stdout=subprocess.STDOUT)

        # Get the path for the synthetic data
        synthetic_genome = (
            synthetic_data_dir / f"synthetic_input_{sample_size}_samples.genome"
        )
        synthetic_segments = (
            synthetic_data_dir / f"synthetic_segments_{sample_size}_samples.txt"
        )

        # We can now run the perl command
        perl_cmd = [
            "perl",
            "-Ilib/perl_modules",
            "t/performance/profile_perl_memory.pl",
            f"{synthetic_genome}",
        ]

        perl_result = subprocess.run(perl_cmd, capture_output=True, text=True)


if __name__ == "__main__":
    run_benchmarks()
