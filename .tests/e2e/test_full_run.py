import pytest
import shutil
import subprocess as sp
import sys
import tempfile
from pathlib import Path


@pytest.fixture
def setup():
    temp_dir = Path(tempfile.mkdtemp())

    reads_fp = Path(".tests/data/reads/").resolve()
    dummy_wolr_fp = temp_dir / "wolr"
    dummy_wolr_fp.mkdir()
    dummy_phy_fp = temp_dir / "dummy_phy"
    dummy_phy_fp.touch()

    project_dir = temp_dir / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"sbx_shotgun_unifrac: {{wolr_fp: '{str(dummy_wolr_fp)}'}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_shotgun_unifrac: {{tree_fp: '{str(dummy_phy_fp)}'}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup
    output_fp = project_dir / "sunbeam_output"
    log_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    # Run the test job
    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_template",
            "--directory",
            temp_dir,
            "-n",
        ]
    )

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    assert True  # Dryrun successful
