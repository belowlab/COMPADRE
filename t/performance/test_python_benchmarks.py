from pathlib import Path
from subprocess import CompletedProcess
import sys
import json
import pytest

# Add the project root to sys.path to import compadre.py
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from compadre import load_segment_information

# Load benchmark configuration
CONFIG_PATH = project_root / "benchmark_config.json"
if not CONFIG_PATH.exists():
    raise FileNotFoundError(f"Unable to find the configuration file: {CONFIG_PATH}")

with open(CONFIG_PATH, "r") as f:
    BENCH_CONFIG = json.load(f)

RUN_OPTS = {s["name"]: s for s in BENCH_CONFIG["params"]}

RUN_LABELS = list(RUN_OPTS.keys())


@pytest.fixture(params=RUN_LABELS)
def segment_file(request):
    scale = request.param
    run_opts = RUN_OPTS[scale]
    n_samples = run_opts["samples"]
    output_dir = Path(BENCH_CONFIG.get("output_dir", "benchmark_data"))

    # Check for uncompressed file path
    file_path = output_dir / f"synthetic_segments_{n_samples}_samples.txt"

    # Fallback to gzipped file if present
    if not file_path.exists() and not file_path.with_suffix(".txt..gz").exists():
        raise FileNotFoundError(
            f"Failed to find the file {file_path} or {file_path}.gz. This like means the setup script was not run. Make sure to run the generate_synthetic_data.py script before running this script"
        )
    elif file_path.with_suffix(".txt.gz").exists():
        file_path = file_path.with_suffix(".txt.gz")

    return str(file_path)


def test_benchmark_segment_loading(benchmark, segment_file):
    """
    Benchmark the loading and parsing of the entire segment file.
    This measures the file parsing throughput and dictionary construction time.
    """
    segment_dict, ibd2_status = benchmark(
        load_segment_information, segment_file, min_cm_options=2.5
    )
    assert isinstance(segment_dict, dict)


def test_benchmark_segment_lookup(benchmark, segment_file):
    """
    Benchmark the lookup speed of a single pair key.
    This simulates standard query resolution time on the socket server.
    """
    # Load the segment data once outside the benchmark block
    segment_dict, _ = load_segment_information(segment_file, min_cm_options=2.5)

    # Grab an existing key, or fallback to default
    test_key = list(segment_dict.keys())[0] if segment_dict else "id1:id2"

    # Benchmark the dictionary lookup operation
    result = benchmark(segment_dict.get, test_key)
    assert result is not None
