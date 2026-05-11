"""
Structure preparation utilities to add hydrogens and generate PQR files
without relying on external pdb2pqr binaries.
"""

import logging
import shutil
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple
import subprocess

from .van_der_waals import VanDerWaalsRadii

logger = logging.getLogger(__name__)


@dataclass
class PreparedStructure:
    """
    Container for intermediate files generated during structure preparation.
    """

    pdb_with_h: Path
    pqr_file: Path
    workdir: Path

    def cleanup(self) -> None:
        """
        Remove the temporary working directory created during preparation.
        """
        if self.workdir.exists():
            shutil.rmtree(self.workdir, ignore_errors=True)


def _normalise_radius_overrides(overrides: Optional[Iterable[str]]) -> Dict[str, float]:
    """
    Convert an iterable of key=value strings into a dictionary of overrides.
    """
    if overrides is None:
        return {}
    result: Dict[str, float] = {}
    for item in overrides:
        key, value = item.split("=")
        result[key.strip().upper()] = float(value)
    return result


def _split_radius_overrides(
    overrides: Dict[str, float],
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """
    Split overrides into element-based and atom-name-based dictionaries.
    Keys with the prefix 'ATOM:' (case-insensitive) or longer than two characters
    are treated as atom-name overrides.
    """
    element_overrides: Dict[str, float] = {}
    atom_overrides: Dict[str, float] = {}
    for raw_key, value in overrides.items():
        key = raw_key.strip().upper()
        if key.startswith("ATOM:"):
            atom_overrides[key[5:].strip()] = value
        elif len(key) <= 2:
            element_overrides[key] = value
        else:
            atom_overrides[key] = value
    return element_overrides, atom_overrides


def add_hydrogens_with_pdbfixer(
    pdb_in: Path,
    pdb_out: Path,
    *,
    ph: float = 7.0,
    keep_waters: bool = True,
    add_missing_heavy: bool = True,
    add_missing_loops: bool = True,
    add_missing_terminals: bool = True,
    remove_heterogens: bool = False,
) -> None:
    """
    Add missing atoms and hydrogens using PDBFixer and write the result to pdb_out.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        raise ImportError(
            "pdbfixer and openmm are required for structure preparation."
        ) from exc

    logger.info("Loading PDB with PDBFixer: %s", pdb_in)
    fixer = PDBFixer(filename=str(pdb_in))

    if remove_heterogens:
        fixer.removeHeterogens(keepWater=keep_waters)

    if add_missing_loops or add_missing_terminals:
        fixer.findMissingResidues()
    fixer.findMissingAtoms()
    if add_missing_heavy:
        fixer.addMissingAtoms()

    if not keep_waters:
        fixer.removeWater()

    logger.info("Adding hydrogens at pH %.2f", ph)
    fixer.addMissingHydrogens(pH=ph)

    logger.info("Writing hydrogenated PDB: %s", pdb_out)
    with open(pdb_out, "w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)


def _resolve_atom_radius(
    atom_name: str,
    element: str,
    element_radii: Dict[str, float],
    atom_radii: Optional[Dict[str, float]],
    default_radius: float,
    record: str = "",
) -> float:
    """
    Determine the radius for an atom using atom-name overrides first, then element data.
    """
    atom_key = atom_name.strip().upper()
    record_key = (record or "").strip().upper()
    element_key_temp = (element or "").strip().upper()

    if atom_key and atom_radii and atom_key in atom_radii:
        if not (
            element_key_temp
            and element_key_temp == atom_key
            and element_key_temp in VanDerWaalsRadii.VDW_RADII
            and record_key == "HETATM"
        ):
            return float(atom_radii[atom_key])

    if record_key == "ATOM" and atom_key:
        canonical_map = {
            "C": element_radii.get("C", VanDerWaalsRadii.get_radius("C", default_radius)),
            "H": element_radii.get("H", VanDerWaalsRadii.get_radius("H", default_radius)),
            "N": element_radii.get("N", VanDerWaalsRadii.get_radius("N", default_radius)),
            "O": element_radii.get("O", VanDerWaalsRadii.get_radius("O", default_radius)),
            "S": element_radii.get("S", VanDerWaalsRadii.get_radius("S", default_radius)),
        }
        prefix = atom_key[0]
        if prefix in canonical_map:
            return float(canonical_map[prefix])

    element_key = element_key_temp
    if element_key and element_key in element_radii:
        return float(element_radii[element_key])

    if atom_key and atom_radii and atom_key in atom_radii:
        return float(atom_radii[atom_key])

    return float(VanDerWaalsRadii.get_radius(element_key or atom_key, default_radius))


def pdb_to_pqr(
    pdb_in: Path,
    pqr_out: Path,
    *,
    radius_table: Optional[Dict[str, float]] = None,
    atom_radius_table: Optional[Dict[str, float]] = None,
    default_radius: float = 1.5,
    default_charge: float = 0.0,
) -> None:
    """
    Convert a PDB file into a PQR file using the supplied radius table.
    """
    if radius_table is None:
        radius_table = VanDerWaalsRadii.VDW_RADII
    if atom_radius_table is None:
        atom_radius_table = VanDerWaalsRadii.ATOM_NAME_RADII

    logger.info("Converting PDB to PQR: %s -> %s", pdb_in, pqr_out)
    with open(pdb_in, "r") as fin, open(pqr_out, "w") as fout:
        for line in fin:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                fout.write(line)
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                fout.write(line)
                continue

            serial = line[6:11]
            name = line[12:16]
            alt_loc = line[16:17]
            res_name = line[17:20]
            chain_id = line[21:22]
            res_seq = line[22:26]
            i_code = line[26:27]

            element = VanDerWaalsRadii.guess_element_from_pdb_line(line)
            radius = _resolve_atom_radius(
                name.strip(),
                element,
                radius_table,
                atom_radius_table,
                default_radius,
                record,
            )
            charge = float(default_charge)

            pqr_line = (
                f"{record:<6}{serial:>5} {name:<4}{alt_loc}{res_name:>3} {chain_id}"
                f"{res_seq:>4}{i_code}   {x:>8.3f}{y:>8.3f}{z:>8.3f} "
                f"{charge:>7.4f} {radius:>7.4f}\n"
            )
            fout.write(pqr_line)


def pdb_add_radii_in_bfactor(
    pdb_in: Path,
    pdb_out: Path,
    *,
    radius_table: Optional[Dict[str, float]] = None,
    atom_radius_table: Optional[Dict[str, float]] = None,
    default_radius: float = 1.5,
    zero_occupancy: bool = True,
) -> None:
    """
    Write van der Waals radii into the B-factor field of a PDB file.
    """
    if radius_table is None:
        radius_table = VanDerWaalsRadii.VDW_RADII
    if atom_radius_table is None:
        atom_radius_table = VanDerWaalsRadii.ATOM_NAME_RADII

    logger.info("Embedding radii in B-factor column: %s -> %s", pdb_in, pdb_out)
    with open(pdb_in, "r") as fin, open(pdb_out, "w") as fout:
        for line in fin:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                fout.write(line)
                continue

            element = VanDerWaalsRadii.guess_element_from_pdb_line(line)
            radius = _resolve_atom_radius(
                line[12:16].strip(),
                element,
                radius_table,
                atom_radius_table,
                default_radius,
                record,
            )

            row = list(line.rstrip("\n"))
            if len(row) < 66:
                row += [" "] * (66 - len(row))

            occupancy = f"{0.00:6.2f}" if zero_occupancy else line[54:60]
            for idx, char in enumerate(occupancy):
                row[54 + idx] = char

            bfactor = f"{radius:6.2f}"
            for idx, char in enumerate(bfactor):
                row[60 + idx] = char

            fout.write("".join(row) + "\n")


def prepare_structure(
    pdb_in: Path,
    *,
    ph: float = 7.0,
    default_radius: float = 1.5,
    default_charge: float = 0.0,
    radius_overrides: Optional[Iterable[str]] = None,
    keep_waters: bool = True,
    add_missing_heavy: bool = True,
    add_missing_loops: bool = True,
    add_missing_terminals: bool = True,
    remove_heterogens: bool = False,
    temp_root: Optional[Path] = None,
    use_external_pdb2pqr: bool = False,
    pdb2pqr_force_field: str = "PARSE",
    pdb2pqr_extra_flags: Optional[Iterable[str]] = None,
) -> PreparedStructure:
    """
    High-level helper that performs hydrogenation and PQR generation.
    """
    workdir = Path(
        tempfile.mkdtemp(prefix="sem_prep_", dir=str(temp_root) if temp_root else None)
    )

    pdb_in = Path(pdb_in)
    pdb_with_h = workdir / f"{pdb_in.stem}_withH.pdb"
    pqr_out = workdir / f"{pdb_in.stem}.pqr"

    radius_table = dict(VanDerWaalsRadii.VDW_RADII)
    atom_radius_table = dict(VanDerWaalsRadii.ATOM_NAME_RADII)
    overrides = _normalise_radius_overrides(radius_overrides)
    element_overrides, atom_overrides = _split_radius_overrides(overrides)
    radius_table.update(element_overrides)
    atom_radius_table.update(atom_overrides)

    if use_external_pdb2pqr:
        logger.info(
            "Running external pdb2pqr (force field: %s) for hydrogenation.",
            pdb2pqr_force_field,
        )
        pdb2pqr_cmd = [
            "pdb2pqr",
            "--ff",
            pdb2pqr_force_field,
        ]
        extra_flags = list(pdb2pqr_extra_flags) if pdb2pqr_extra_flags else [
            "--nodebump",
            "--noopt",
        ]
        pdb2pqr_cmd.extend(extra_flags)
        generated_pqr = workdir / f"{pdb_in.stem}_pdb2pqr_raw.pqr"
        pdb2pqr_cmd.extend([str(pdb_in), str(generated_pqr)])

        try:
            subprocess.run(
                pdb2pqr_cmd,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            logger.error("pdb2pqr failed: %s", exc.stderr)
            raise RuntimeError("pdb2pqr execution failed") from exc

        _rewrite_pqr_with_custom_radii(
            generated_pqr,
            pqr_out,
            element_radius_table=radius_table,
            atom_radius_table=atom_radius_table,
            default_radius=default_radius,
        )

        try:
            import MDAnalysis as mda

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                universe = mda.Universe(str(pqr_out))
                universe.atoms.write(str(pdb_with_h))
        except Exception as exc:
            logger.warning(
                "Could not generate PDB from PQR using MDAnalysis: %s; "
                "using PQR as placeholder for pdb_with_h.",
                exc,
            )
            pdb_with_h.write_text(pqr_out.read_text())

        if generated_pqr.exists():
            generated_pqr.unlink(missing_ok=True)

    else:
        add_hydrogens_with_pdbfixer(
            pdb_in,
            pdb_with_h,
            ph=ph,
            keep_waters=keep_waters,
            add_missing_heavy=add_missing_heavy,
            add_missing_loops=add_missing_loops,
            add_missing_terminals=add_missing_terminals,
            remove_heterogens=remove_heterogens,
        )

        pdb_to_pqr(
            pdb_with_h,
            pqr_out,
            radius_table=radius_table,
            atom_radius_table=atom_radius_table,
            default_radius=default_radius,
            default_charge=default_charge,
        )

    return PreparedStructure(
        pdb_with_h=pdb_with_h,
        pqr_file=pqr_out,
        workdir=workdir,
    )


__all__ = [
    "PreparedStructure",
    "add_hydrogens_with_pdbfixer",
    "pdb_to_pqr",
    "pdb_add_radii_in_bfactor",
    "prepare_structure",
]


def _parse_pqr_atom_line(line: str) -> Dict[str, Any]:
    """
    Parse a PQR atom line into its constituent fields using robust tokenisation.
    Raises ValueError if the line cannot be interpreted.
    """
    record = line[:6].strip()
    if record not in ("ATOM", "HETATM"):
        raise ValueError("Not an atom or heteroatom record")

    serial = line[6:11].strip()
    atom_name = line[12:16].strip()
    alt_loc = line[16:17]
    res_name = line[17:20].strip()
    chain_id = line[21:22].strip()
    res_seq = line[22:26].strip()
    i_code = line[26:27]

    numeric_fields = line[30:].split()
    if len(numeric_fields) < 5:
        raise ValueError("Failed to parse numeric PQR fields")

    try:
        x = float(numeric_fields[0])
        y = float(numeric_fields[1])
        z = float(numeric_fields[2])
        charge = float(numeric_fields[3])
        radius = float(numeric_fields[4])
    except (TypeError, ValueError) as exc:
        raise ValueError("Failed to parse numeric PQR fields") from exc

    if not serial:
        raise ValueError("Missing atom serial number")

    if not atom_name:
        # Attempt to recover atom name from serial token if needed (e.g., '54HH11')
        serial_token = line[6:].split()[0]
        digits = "".join(ch for ch in serial_token if ch.isdigit())
        suffix = serial_token[len(digits):]
        if digits and suffix:
            serial = digits
            atom_name = suffix

    if not atom_name:
        raise ValueError("Missing atom name in PQR line")

    return {
        "record": record,
        "serial": serial,
        "atom_name": atom_name,
        "alt_loc": alt_loc,
        "res_name": res_name,
        "chain_id": chain_id,
        "res_seq": res_seq,
        "i_code": i_code,
        "x": x,
        "y": y,
        "z": z,
        "charge": charge,
        "radius": radius,
    }


def _format_serial_field(serial: str) -> str:
    """
    Format the atom serial field ensuring at least one leading space.
    """
    serial_clean = str(serial).strip()
    if not serial_clean:
        return "      "
    try:
        serial_int = int(float(serial_clean))
        field = f"{serial_int:>5}"
    except ValueError:
        field = serial_clean[-5:]
        field = field.rjust(5)
    if not field.startswith(" "):
        field = " " + field
    else:
        field = field.rjust(6)
    return field


def _format_atom_name_field(atom_name: str) -> str:
    """
    Format the atom name according to PDB alignment heuristics.
    """
    name = atom_name.strip()
    if not name:
        return "    "
    if len(name) == 4:
        return name[:4]
    if name[0].isdigit():
        return f"{name:<4}"
    return f"{name:>4}"


def _format_res_seq_field(res_seq: str) -> str:
    """
    Format the residue sequence number into a 4-character field.
    """
    seq = res_seq.strip()
    if len(seq) > 4:
        seq = seq[-4:]
    return f"{seq:>4}" if seq else "   1"


def _format_pqr_atom_line(fields: Dict[str, Any], radius: float) -> str:
    """
    Build a normalised PQR atom line from parsed fields.
    """
    record = fields["record"]
    serial_field = _format_serial_field(fields["serial"])
    name_field = _format_atom_name_field(fields["atom_name"])
    alt_loc = fields.get("alt_loc", " ")
    if not alt_loc or alt_loc == "":
        alt_loc = " "
    res_name_field = f"{fields['res_name'].strip()[:3]:>3}"
    chain_value = (fields.get("chain_id") or "").strip() or "X"
    chain_field = chain_value[:1]
    res_seq_field = _format_res_seq_field(fields.get("res_seq", ""))
    i_code = fields.get("i_code", " ")
    if not i_code or i_code == "":
        i_code = " "

    x = float(fields["x"])
    y = float(fields["y"])
    z = float(fields["z"])
    charge = float(fields["charge"])

    return (
        f"{record:<6}{serial_field} {name_field}{alt_loc}{res_name_field} {chain_field}"
        f"{res_seq_field}{i_code}   {x:>8.3f}{y:>8.3f}{z:>8.3f} {charge:>8.4f} {radius:>7.4f}\n"
    )


def _rewrite_pqr_with_custom_radii(
    pqr_in: Path,
    pqr_out: Path,
    element_radius_table: Dict[str, float],
    atom_radius_table: Dict[str, float],
    default_radius: float,
) -> None:
    """
    Rewrite PQR file replacing radius column with custom values from radius_table.
    """
    with open(pqr_in, "r") as fin, open(pqr_out, "w") as fout:
        for line in fin:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                fout.write(line)
                continue

            try:
                fields = _parse_pqr_atom_line(line)
            except ValueError:
                fout.write(line)
                continue

            element = VanDerWaalsRadii.guess_element_from_pdb_line(line)
            radius = _resolve_atom_radius(
                fields["atom_name"],
                element,
                element_radius_table,
                atom_radius_table,
                default_radius,
                record,
            )

            try:
                updated_line = _format_pqr_atom_line(fields, radius)
            except Exception:
                fout.write(line)
                continue

            fout.write(updated_line)
