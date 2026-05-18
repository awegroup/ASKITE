"""
Adjust PSM wing element rest lengths to match current nodal positions.

By default this updates the mapped stretched PSM YAML in-place. Only the l0
field in wing_elements is rewritten; all other YAML text is preserved.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import numpy as np
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_YAML_PATH = (
    SCRIPT_DIR
    / "results"
    / "struc_geometry_PSM_reduced_stretched_mapped_to_FEM.yaml"
)


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def wing_node_coordinates(geometry: dict[str, Any]) -> dict[int, np.ndarray]:
    return {
        int(row[0]): np.array([float(row[1]), float(row[2]), float(row[3])])
        for row in geometry["wing_particles"]["data"]
    }


def wing_connection_lengths(geometry: dict[str, Any]) -> dict[str, list[float]]:
    coords = wing_node_coordinates(geometry)
    lengths_by_name: dict[str, list[float]] = {}

    for row in geometry["wing_connections"]["data"]:
        name = str(row[0])
        ci = int(row[1])
        cj = int(row[2])
        if ci not in coords or cj not in coords:
            raise ValueError(
                f"wing_connection {name} references missing node(s): {ci}, {cj}"
            )
        length = float(np.linalg.norm(coords[cj] - coords[ci]))
        lengths_by_name.setdefault(name, []).append(length)

    return lengths_by_name


def mean_l0_by_element_name(
    lengths_by_name: dict[str, list[float]],
    *,
    duplicate_tolerance: float,
) -> dict[str, float]:
    l0_by_name: dict[str, float] = {}
    for name, lengths in sorted(lengths_by_name.items()):
        values = np.asarray(lengths, dtype=float)
        spread = float(values.max() - values.min())
        if len(values) > 1 and spread > duplicate_tolerance:
            print(
                "[adjust_wing_element_lenghts_to_nodal_position] warning: "
                f"{name} has {len(values)} connection lengths with spread "
                f"{spread:.6e} m; using mean."
            )
        l0_by_name[name] = float(values.mean())
    return l0_by_name


def split_line_ending(line: str) -> tuple[str, str]:
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n"):
        return line[:-1], "\n"
    return line, ""


def format_l0(value: float) -> str:
    value = float(value)
    if abs(value) < 5e-13:
        value = 0.0
    return f"{value:.6f}"


def replace_l0_in_wing_element_row(
    body: str,
    l0_by_name: dict[str, float],
) -> tuple[str, str | None]:
    row_re = re.compile(r"^(\s*-\s*)\[(.*)\](.*)$")
    match = row_re.match(body)
    if match is None:
        return body, None

    inner = match.group(2)
    values = [value.strip() for value in inner.split(",")]
    if len(values) < 2:
        return body, None

    name = values[0]
    if name not in l0_by_name:
        return body, None

    values[1] = format_l0(l0_by_name[name])
    return f"{match.group(1)}[{', '.join(values)}]{match.group(3)}", name


def update_wing_element_l0_text(
    input_path: Path,
    output_path: Path,
    l0_by_name: dict[str, float],
) -> set[str]:
    top_level_key_re = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\s*:")
    active_wing_elements = False
    updated_names: set[str] = set()
    output_lines: list[str] = []

    for line in input_path.read_text(encoding="utf-8").splitlines(keepends=True):
        body, ending = split_line_ending(line)
        top_level_match = top_level_key_re.match(body)
        if top_level_match is not None:
            active_wing_elements = top_level_match.group(1) == "wing_elements"

        if active_wing_elements:
            body, updated_name = replace_l0_in_wing_element_row(body, l0_by_name)
            if updated_name is not None:
                updated_names.add(updated_name)

        output_lines.append(body + ending)

    missing = sorted(set(l0_by_name) - updated_names)
    if missing:
        raise ValueError(f"Did not update wing_elements rows for: {missing}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(output_lines), encoding="utf-8")
    return updated_names


def adjust_wing_element_lengths(
    yaml_path: Path,
    output_path: Path,
    *,
    duplicate_tolerance: float,
) -> set[str]:
    geometry = load_yaml(yaml_path)
    lengths_by_name = wing_connection_lengths(geometry)
    l0_by_name = mean_l0_by_element_name(
        lengths_by_name,
        duplicate_tolerance=duplicate_tolerance,
    )
    return update_wing_element_l0_text(yaml_path, output_path, l0_by_name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "yaml_path",
        nargs="?",
        type=Path,
        default=DEFAULT_YAML_PATH,
        help="YAML file to update.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output YAML path. Defaults to overwriting yaml_path.",
    )
    parser.add_argument(
        "--duplicate-tolerance",
        type=float,
        default=1e-9,
        help="Warn when repeated wing element names differ by more than this.",
    )
    args = parser.parse_args()

    yaml_path = args.yaml_path.resolve()
    output_path = args.output.resolve() if args.output is not None else yaml_path
    updated_names = adjust_wing_element_lengths(
        yaml_path,
        output_path,
        duplicate_tolerance=args.duplicate_tolerance,
    )

    print(
        "[adjust_wing_element_lenghts_to_nodal_position] updated "
        f"{len(updated_names)} wing element l0 values in {output_path}"
    )


if __name__ == "__main__":
    main()
