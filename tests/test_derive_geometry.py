"""
Behavioural (round-trip) tests for `sem.scripts.derive_geometry`: fitting a
parametric pore (or a radius-vs-z `profile` table) back out of a synthetic
all-atom-style conductivity/distance map built from a known `PoreProfile`.

Numpy/scipy only -- no dolfinx dependency anywhere in this file. Grids are
kept to 81**3 (box +-40 A, 1 A spacing) so the whole file runs quickly.
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

import sem.scripts.derive_geometry as dg
from sem.config import create_example_config, validate_config
from sem.geometry_profiles import PoreProfile
from sem.grid_io import condfrac, write_binary_file

GRID_HALF = 40.0
SPACING = 1.0


def _make_bin(tmp_path, profile, *, name="pore.bin", dx=0.0, dy=0.0):
    """Build a synthetic .bin distance map from a PoreProfile on a box
    +-40 A, 1 A grid (81**3), optionally shifting the pore off-axis by
    (dx, dy)."""
    n = int(round(2 * GRID_HALF / SPACING)) + 1
    coords = np.linspace(-GRID_HALF, GRID_HALF, n)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing="ij")
    R = np.sqrt((X - dx) ** 2 + (Y - dy) ** 2)
    d = profile.distance(R, Z)
    path = tmp_path / name
    write_binary_file(
        d.transpose(2, 1, 0).astype(np.float32),
        origin=(-GRID_HALF, -GRID_HALF, -GRID_HALF),
        resolution=SPACING,
        filename=str(path),
    )
    return str(path)


def _make_conductivity_bin(tmp_path, profile, *, name="pore_sigma.bin", bulk=10.5):
    n = int(round(2 * GRID_HALF / SPACING)) + 1
    coords = np.linspace(-GRID_HALF, GRID_HALF, n)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing="ij")
    R = np.sqrt(X ** 2 + Y ** 2)
    d = profile.distance(R, Z)
    sigma = bulk * condfrac(d)
    path = tmp_path / name
    write_binary_file(
        sigma.transpose(2, 1, 0).astype(np.float32),
        origin=(-GRID_HALF, -GRID_HALF, -GRID_HALF),
        resolution=SPACING,
        filename=str(path),
    )
    return str(path)


# ---------------------------------------------------------------------------
# Static guard
# ---------------------------------------------------------------------------
def test_module_source_has_no_dolfinx_or_utils_import():
    source = Path(dg.__file__).read_text()
    assert not re.search(r"^\s*(import|from)\s+.*\bdolfinx\b", source, re.MULTILINE)
    assert not re.search(r"^\s*from\s+\.\.\s+import\s+utils\b", source, re.MULTILINE)
    assert not re.search(r"^\s*from\s+sem\.utils\b", source, re.MULTILINE)
    assert not re.search(r"^\s*from\s+\.\.utils\b", source, re.MULTILINE)


# ---------------------------------------------------------------------------
# Round-trip: plain cylinder
# ---------------------------------------------------------------------------
def test_derive_cylindrical_round_trip(tmp_path):
    profile = PoreProfile.cylindrical(pore_radius=15.0, half_thickness=20.0)
    bin_path = _make_bin(tmp_path, profile, name="cyl.bin")

    result = dg.derive(bin_path, "cylindrical")

    assert result.pore_geometry["pore_type"] == "cylindrical"
    assert result.pore_geometry["pore_radius"] == pytest.approx(15.0, abs=0.05)
    assert result.pore_geometry["membrane_thickness"] == pytest.approx(40.0, abs=0.1)
    assert abs(result.derivation["predicted_deviation_pct"]) < 0.2


# ---------------------------------------------------------------------------
# Round-trip: double cone
# ---------------------------------------------------------------------------
def test_derive_double_cone_round_trip(tmp_path):
    profile = PoreProfile.double_cone(inner_radius=12.0, outer_radius=24.0, half_thickness=20.0)
    bin_path = _make_bin(tmp_path, profile, name="dcone.bin")

    result = dg.derive(bin_path, "double_cone")

    assert result.pore_geometry["pore_type"] == "double_cone"
    assert result.pore_geometry["pore_radius"] == pytest.approx(12.0, abs=0.2)
    assert result.pore_geometry["outer_radius"] == pytest.approx(24.0, abs=0.2)


# ---------------------------------------------------------------------------
# Round-trip: conical
# ---------------------------------------------------------------------------
def test_derive_conical_round_trip(tmp_path):
    profile = PoreProfile.conical(bottom_radius=10.0, top_radius=18.0, half_thickness=20.0)
    bin_path = _make_bin(tmp_path, profile, name="cone.bin")

    result = dg.derive(bin_path, "conical")

    assert result.pore_geometry["pore_type"] == "conical"
    assert result.pore_geometry["bottom_radius"] == pytest.approx(10.0, abs=0.2)
    assert result.pore_geometry["top_radius"] == pytest.approx(18.0, abs=0.2)


# ---------------------------------------------------------------------------
# Round-trip: chamfered cylinder (--fit-corner)
# ---------------------------------------------------------------------------
def test_derive_cylindrical_with_corner_round_trip(tmp_path):
    profile = PoreProfile.cylindrical(
        pore_radius=15.0, half_thickness=20.0, corner_radius=5.0, chamfer_depth=5.0
    )
    bin_path = _make_bin(tmp_path, profile, name="chamfer.bin")

    result = dg.derive(bin_path, "cylindrical", fit_corner=True)

    assert result.pore_geometry["corner_radius"] == pytest.approx(5.0, abs=1.0)
    assert result.pore_geometry["chamfer_depth"] == pytest.approx(5.0, abs=1.0)


# ---------------------------------------------------------------------------
# `profile` output reproduces A_bin(z) on membrane slices
# ---------------------------------------------------------------------------
def test_derive_profile_reproduces_area_bin(tmp_path):
    # A conical (single-frustum) wall: sloped, so the profile table
    # exercises the slope-aware R_map inversion, but -- unlike double_cone
    # -- it has no apex kink, where a discontinuous true slope would
    # defeat the smoothed local-slope estimate used to build the table.
    profile = PoreProfile.conical(bottom_radius=10.0, top_radius=18.0, half_thickness=20.0)
    bin_path = _make_bin(tmp_path, profile, name="profile_src.bin")

    result = dg.derive(bin_path, "profile")
    z_tab, r_tab = result.profile_table
    half_L = result.derivation["membrane_thickness_used"] / 2.0
    assert len(z_tab) == len(r_tab)
    assert z_tab[0] == pytest.approx(-half_L, abs=1e-6)
    assert z_tab[-1] == pytest.approx(half_L, abs=1e-6)

    rebuilt = PoreProfile.from_table(z_tab, r_tab, half_thickness=half_L)
    slices = result.slices
    z_center = result.derivation["z_center"]
    # Restrict to slices genuinely inside the fitted membrane extent: the
    # f_far < 0.5 threshold also catches a couple of Angstroms of ramp
    # overhang just past the true face (RAMP_MIN/RAMP_MAX wide), which
    # neither a table nor a parametric fit is meant to reproduce.
    membrane = np.abs(slices.z - z_center) <= half_L

    x = np.linspace(-GRID_HALF, GRID_HALF, int(round(2 * GRID_HALF / SPACING)) + 1)
    y = x.copy()
    A_param = dg.parametric_area_profile(
        rebuilt, (x, y, slices.z - z_center), slices.spacing
    )
    rel_err = np.abs(A_param[membrane] - slices.area_bin[membrane]) / slices.area_bin[membrane]
    assert np.max(rel_err) < 0.005


# ---------------------------------------------------------------------------
# Conductivity-units bin
# ---------------------------------------------------------------------------
def test_derive_conductivity_units_matches_distance_units(tmp_path):
    profile = PoreProfile.cylindrical(pore_radius=15.0, half_thickness=20.0)
    bin_path = _make_conductivity_bin(tmp_path, profile, name="cyl_sigma.bin", bulk=10.5)

    result = dg.derive(bin_path, "cylindrical", units="conductivity")

    assert result.pore_geometry["pore_radius"] == pytest.approx(15.0, abs=0.05)


# ---------------------------------------------------------------------------
# Off-axis pore
# ---------------------------------------------------------------------------
def test_derive_off_axis_shift_recovers_offset(tmp_path):
    profile = PoreProfile.cylindrical(pore_radius=15.0, half_thickness=20.0)
    bin_path = _make_bin(tmp_path, profile, name="offaxis.bin", dx=2.0, dy=0.0)

    result = dg.derive(bin_path, "cylindrical")

    assert result.pore_geometry["pore_radius"] == pytest.approx(15.0, abs=0.05)
    x0, y0 = result.derivation["axis_offset"]
    assert x0 == pytest.approx(2.0, abs=0.3)
    assert y0 == pytest.approx(0.0, abs=0.3)


# ---------------------------------------------------------------------------
# CLI: main() writes the four files; JSON block merges into a valid config
# ---------------------------------------------------------------------------
def test_main_writes_files_and_config_merge_validates(tmp_path):
    profile = PoreProfile.cylindrical(pore_radius=15.0, half_thickness=20.0)
    bin_path = _make_bin(tmp_path, profile, name="cli.bin")
    prefix = str(tmp_path / "derived")

    rc = dg.main([bin_path, "--pore-type", "cylindrical", "--output-prefix", prefix])
    assert rc == 0

    geom_path = Path(prefix + "_pore_geometry.json")
    profile_csv = Path(prefix + "_profile.csv")
    profile_pore_csv = Path(prefix + "_profile_pore.csv")
    derivation_path = Path(prefix + "_derivation.json")
    assert geom_path.is_file()
    assert profile_csv.is_file()
    assert profile_pore_csv.is_file()
    assert derivation_path.is_file()

    with open(geom_path) as fh:
        pore_geometry = json.load(fh)
    with open(derivation_path) as fh:
        derivation = json.load(fh)
    assert derivation["predicted_deviation_pct"] is not None

    template_path = tmp_path / "template.json"
    create_example_config("cylindrical", str(template_path))
    with open(template_path) as fh:
        cfg = json.load(fh)
    cfg["pore_geometry"] = pore_geometry
    assert validate_config(cfg, require_analyte=False)


# ---------------------------------------------------------------------------
# CLI help works through both entry points
# ---------------------------------------------------------------------------
def test_main_help_exits_cleanly():
    with pytest.raises(SystemExit) as excinfo:
        dg.main(["--help"])
    assert excinfo.value.code == 0
