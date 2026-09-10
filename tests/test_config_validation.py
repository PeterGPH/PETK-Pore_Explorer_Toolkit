"""
Behavioural tests for the `distance_metric` / `chamfer_depth` / `profile`
config plumbing in `sem.config.validate_config` and `create_example_config`,
plus the `create_config` CLI subcommand in `sem.cli`. Imports only
`sem.config` (and `sem.cli`, whose argument-parsing path must also stay
importable without dolfinx) -- no dolfinx dependency anywhere in this file.
"""

import json

import pytest

import sem.cli
from sem.config import create_example_config, validate_config


def _example_config(tmp_path, pore_type):
    """Load an example config for `pore_type` via the real template writer,
    so these tests exercise the same defaults `create_config` produces."""
    out = tmp_path / f"{pore_type}_example.json"
    create_example_config(pore_type, str(out))
    with open(out) as fh:
        return json.load(fh)


def _write_csv(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write(",".join(str(v) for v in row) + "\n")


# ---------------------------------------------------------------------------
# distance_metric
# ---------------------------------------------------------------------------
def test_distance_metric_defaults_to_euclidean_when_absent(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    del cfg["pore_geometry"]["distance_metric"]
    assert validate_config(cfg, require_analyte=False)
    assert cfg["pore_geometry"]["distance_metric"] == "euclidean"


def test_distance_metric_is_lower_cased_and_written_back(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"]["distance_metric"] = "EUCLIDEAN"
    assert validate_config(cfg, require_analyte=False)
    assert cfg["pore_geometry"]["distance_metric"] == "euclidean"


def test_invalid_distance_metric_value_fails(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"]["distance_metric"] = "radial"
    assert validate_config(cfg, require_analyte=False) is False


# ---------------------------------------------------------------------------
# corner_radius / chamfer_depth
# ---------------------------------------------------------------------------
def test_negative_corner_radius_fails(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"]["corner_radius"] = -5.0
    assert validate_config(cfg, require_analyte=False) is False


def test_chamfer_depth_non_positive_fails(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"]["corner_radius"] = 10.0

    cfg["pore_geometry"]["chamfer_depth"] = 0.0
    assert validate_config(cfg, require_analyte=False) is False

    cfg["pore_geometry"]["chamfer_depth"] = -1.0
    assert validate_config(cfg, require_analyte=False) is False


def test_chamfer_depth_positive_passes_and_warns_for_non_cylindrical(tmp_path, caplog):
    cfg = _example_config(tmp_path, "double_cone")
    cfg["pore_geometry"]["chamfer_depth"] = 5.0
    with caplog.at_level("WARNING"):
        assert validate_config(cfg, require_analyte=False) is True
    assert any("chamfer_depth" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# membrane_thickness
# ---------------------------------------------------------------------------
def test_missing_membrane_thickness_for_cylindrical_fails(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    del cfg["pore_geometry"]["membrane_thickness"]
    assert validate_config(cfg, require_analyte=False) is False


def test_cylindrical_membrane_thickness_zero_is_invalid(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"]["membrane_thickness"] = 0.0
    assert validate_config(cfg, require_analyte=False) is False


def test_bin_file_membrane_thickness_zero_is_valid(tmp_path):
    cfg = _example_config(tmp_path, "bin_file")
    bin_path = tmp_path / "pore_structure.bin"
    bin_path.write_bytes(b"\x00" * 16)
    cfg["pore_geometry"]["bin_file_path"] = str(bin_path)
    cfg["pore_geometry"]["membrane_thickness"] = 0.0
    assert validate_config(cfg, require_analyte=False) is True


# ---------------------------------------------------------------------------
# profile pore type
# ---------------------------------------------------------------------------
def test_profile_without_profile_path_fails(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"] = {
        "pore_type": "profile",
        "membrane_thickness": 200.0,
    }
    assert validate_config(cfg, require_analyte=False) is False


def test_profile_with_valid_csv_passes(tmp_path):
    csv_path = tmp_path / "profile.csv"
    _write_csv(csv_path, [
        ("z", "r"),
        (-100.0, 50.0),
        (0.0, 25.0),
        (100.0, 50.0),
    ])
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"] = {
        "pore_type": "profile",
        "membrane_thickness": 200.0,
        "profile_path": str(csv_path),
    }
    assert validate_config(cfg, require_analyte=False) is True


def test_profile_with_nonpositive_radius_fails(tmp_path):
    csv_path = tmp_path / "profile.csv"
    _write_csv(csv_path, [
        (-100.0, 50.0),
        (0.0, 0.0),
        (100.0, 50.0),
    ])
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"] = {
        "pore_type": "profile",
        "membrane_thickness": 200.0,
        "profile_path": str(csv_path),
    }
    assert validate_config(cfg, require_analyte=False) is False


def test_profile_with_non_monotonic_z_fails(tmp_path):
    csv_path = tmp_path / "profile.csv"
    _write_csv(csv_path, [
        (-100.0, 50.0),
        (-100.0, 40.0),
        (100.0, 50.0),
    ])
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"] = {
        "pore_type": "profile",
        "membrane_thickness": 200.0,
        "profile_path": str(csv_path),
    }
    assert validate_config(cfg, require_analyte=False) is False


def test_profile_with_missing_file_fails(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    cfg["pore_geometry"] = {
        "pore_type": "profile",
        "membrane_thickness": 200.0,
        "profile_path": str(tmp_path / "does_not_exist.csv"),
    }
    assert validate_config(cfg, require_analyte=False) is False


# ---------------------------------------------------------------------------
# create_example_config
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("pore_type", ["cylindrical", "double_cone", "conical", "profile"])
def test_create_example_config_contains_distance_metric(tmp_path, pore_type):
    cfg = _example_config(tmp_path, pore_type)
    assert cfg["pore_geometry"]["distance_metric"] == "euclidean"


def test_create_example_config_cylindrical_has_chamfer_depth_null(tmp_path):
    cfg = _example_config(tmp_path, "cylindrical")
    assert cfg["pore_geometry"]["chamfer_depth"] is None


def test_create_example_config_profile_has_profile_path(tmp_path):
    cfg = _example_config(tmp_path, "profile")
    assert "profile_path" in cfg["pore_geometry"]
    assert cfg["pore_geometry"]["pore_type"] == "profile"


# ---------------------------------------------------------------------------
# create_config CLI subcommand (argument parsing only -- no dolfinx needed)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("pore_type", ["conical", "profile"])
def test_create_config_cli_accepts_new_pore_types(tmp_path, pore_type):
    out_path = tmp_path / f"{pore_type}_config.json"
    sem.cli.main(["create_config", pore_type, "-o", str(out_path)])
    assert out_path.exists()
    with open(out_path) as fh:
        cfg = json.load(fh)
    assert cfg["pore_geometry"]["pore_type"] == pore_type
    assert cfg["pore_geometry"]["distance_metric"] == "euclidean"
