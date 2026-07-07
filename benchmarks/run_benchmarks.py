#!/usr/bin/env python3.14
# Runner for the test to check the memory and runtime of the perl load
# data function and the python server load_data function
import os
import json
import subprocess
import sys
from pathlib import Path

# We need to get the project root so that we can load the correct functions
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

# Load in the configuration object from the benchmarks folder
CONFIG_PATH = project_root / "benchmarks/benchmark_config.json"
if not CONFIG_PATH.exists():
    raise FileNotFoundError(
        f"Unable to find the configuration json file: {CONFIG_PATH}."
    )

with open(CONFIG_PATH, "r") as config_fh:
    CONFIG = json.load(config_fh)


def run_benchmarks():
    """Main runner script that will run both the perl benchmark and the python benchmark"""
    # Make sure the output directory exists
    synthetic_data_dir = project_root / CONFIG.get("output_dir", "benchmark_data")

    if not synthetic_data_dir.exists():
        synthetic_data_dir.mkdir(parents=True, exist_ok=True)

    for config_opts in CONFIG.get("params"):
        label = config_opts["name"]
        sample_size = config_opts["samples"]
        related_fraction = config_opts["related_fraction"]
        pair_count = config_opts["pairs"]

        print(f"Generating synthetic data for {sample_size} samples")

        # Now we need to run the script to generate the synthetic data
        gen_cmd = [
            sys.executable,
            str(project_root / "benchmarks/generate_synthetic_data.py"),
            "--samples",
            f"{sample_size}",
            "--pairs",
            f"{pair_count}",
            "--related-fraction",
            f"{related_fraction}",
            "--output-dir",
            f"{synthetic_data_dir}",
        ]

        subprocess.run(gen_cmd, check=True)

        # Get the paths for the generated synthetic files
        synthetic_genome = (
            synthetic_data_dir / f"synthetic_input_{sample_size}_samples.genome"
        )
        synthetic_segments = (
            synthetic_data_dir / f"synthetic_segments_{sample_size}_samples.txt"
        )

        # Run the pytest command using pytest-monitor (active by default if installed)
        db_path = CONFIG.get("db", "db/benchmark_metrics.db")
        python_server_label = CONFIG.get("python_server_label", "python_server")
        db_parent = project_root / db_path
        db_parent.parent.mkdir(parents=True, exist_ok=True)

        env = os.environ.copy()
        env["BENCHMARK_CONFIG_PATH"] = str(CONFIG_PATH)

        pytest_cmd = [
            sys.executable,
            "-m",
            "pytest",
            str(project_root / "benchmarks/test_python_benchmarks.py"),
            "-k",
            f"[{label}]",
            f"--db={db_parent}",
            f"--force-component={python_server_label}"
        ]
        
        print(f"Running Python load & lookup benchmarks for {label}...")
        subprocess.run(pytest_cmd, env=env, check=True)

        # Run the perl profiling command
        perl_cmd = [
            "perl",
            "-I",
            str(project_root / "lib/perl_modules"),
            str(project_root / "benchmarks/profile_perl_memory.pl"),
            str(synthetic_genome),
        ]

        print(f"Running Perl load memory benchmark for {label}...")
        subprocess.run(perl_cmd, check=True)


if __name__ == "__main__":
    run_benchmarks()
