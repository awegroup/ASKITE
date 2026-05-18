"""
Map reduced PSM wing nodes onto FEM wing-side bridle attachment endpoints.

The output is a copy of the stretched PSM YAML where only selected
wing_particles coordinates are replaced. PSM forward tip nodes keep their
original x position, but their y and z coordinates are aligned with their
paired rearward tip node.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
GEOMETRY_PATH = SCRIPT_DIR / "results" / "struc_geometry_FEM_full_stretched.yaml"
PSM_GEOMETRY_PATH = SCRIPT_DIR / "results" / "struc_geometry_PSM_reduced_stretched.yaml"
OUTPUT_PATH = (
    SCRIPT_DIR / "results" / "struc_geometry_PSM_reduced_stretched_mapped_to_FEM.yaml"
)

# PSM nodes 1 and 19 are the forward tip nodes. Their x coordinate should stay
# unchanged, while their y and z should match their paired rearward tip nodes.
# PSM_FORWARD_TIP_NODE_IDS = {1, 19}
PSM_FORWARD_TIP_NODE_IDS = {22}


@dataclass
class StationEndpoint:
    station_name: str
    forward_node_id: int
    forward_coord: np.ndarray
    rearward_node_id: int
    rearward_coord: np.ndarray


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def coordinate_map(geometry: dict[str, Any]) -> dict[int, np.ndarray]:
    coords: dict[int, np.ndarray] = {}
    for table_name in ("wing_particles", "bridle_particles"):
        for row in geometry.get(table_name, {}).get("data", []):
            coords[int(row[0])] = np.array(
                [float(row[1]), float(row[2]), float(row[3])],
                dtype=float,
            )
    return coords


def bridle_connection_node_ids(geometry: dict[str, Any]) -> set[int]:
    node_ids: set[int] = set()
    for row in geometry.get("bridle_connections", {}).get("data", []):
        node_ids.update(int(value) for value in row[1:] if int(value) != 0)
    return node_ids


def wing_particle_node_ids(geometry: dict[str, Any]) -> set[int]:
    return {int(row[0]) for row in geometry.get("wing_particles", {}).get("data", [])}


def station_y_mean(node_ids: Iterable[int], coords: dict[int, np.ndarray]) -> float:
    station_coords = [
        coords[int(node_id)] for node_id in node_ids if int(node_id) in coords
    ]
    if not station_coords:
        raise ValueError("Station has no known node coordinates.")
    return float(np.mean([coord[1] for coord in station_coords]))


def fem_station_endpoints(geometry: dict[str, Any]) -> list[StationEndpoint]:
    coords = coordinate_map(geometry)
    connection_nodes = bridle_connection_node_ids(geometry)
    direct_wing_attachments = wing_particle_node_ids(geometry).intersection(
        connection_nodes
    )

    station_entries = []
    assigned_nodes: set[int] = set()
    for station_i, row in enumerate(geometry.get("strut_tubes", {}).get("data", [])):
        station_name = str(row[0])
        station_node_ids = [int(node_id) for node_id in row[6]]
        attachment_ids = [
            node_id for node_id in station_node_ids if node_id in connection_nodes
        ]
        assigned_nodes.update(attachment_ids)
        station_entries.append(
            {
                "station_i": station_i,
                "name": station_name,
                "node_ids": station_node_ids,
                "attachment_ids": attachment_ids,
                "y_mean": station_y_mean(station_node_ids, coords),
            }
        )

    # Include direct bridle attachments to wing_particles, such as the tip-side
    # a5_equiv nodes, by assigning them to the nearest spanwise station.
    for node_id in sorted(direct_wing_attachments):
        if node_id in assigned_nodes:
            continue
        coord = coords[node_id]
        nearest_station = min(
            station_entries,
            key=lambda station: abs(float(coord[1]) - float(station["y_mean"])),
        )
        nearest_station["attachment_ids"].append(node_id)
        assigned_nodes.add(node_id)

    endpoints: list[StationEndpoint] = []
    for station in station_entries:
        attachment_ids = sorted(
            station["attachment_ids"],
            key=lambda node_id: (coords[int(node_id)][0], int(node_id)),
        )
        if not attachment_ids:
            raise ValueError(
                f"No FEM bridle attachment nodes found for station {station['name']}."
            )
        forward_node_id = int(attachment_ids[0])
        rearward_node_id = int(attachment_ids[-1])
        endpoints.append(
            StationEndpoint(
                station_name=str(station["name"]),
                forward_node_id=forward_node_id,
                forward_coord=coords[forward_node_id],
                rearward_node_id=rearward_node_id,
                rearward_coord=coords[rearward_node_id],
            )
        )

    return endpoints


def psm_wing_node_pairs(geometry: dict[str, Any]) -> list[tuple[int, int]]:
    wing_node_ids = sorted(
        int(row[0]) for row in geometry.get("wing_particles", {}).get("data", [])
    )
    if len(wing_node_ids) % 2 != 0:
        raise ValueError("Expected an even number of PSM wing nodes.")
    return [
        (wing_node_ids[i], wing_node_ids[i + 1])
        for i in range(0, len(wing_node_ids), 2)
    ]


def build_psm_to_fem_mapping(
    psm_geometry: dict[str, Any],
    fem_endpoints: list[StationEndpoint],
) -> dict[int, np.ndarray]:
    psm_coords = coordinate_map(psm_geometry)
    psm_pairs = psm_wing_node_pairs(psm_geometry)
    if len(psm_pairs) != len(fem_endpoints):
        raise ValueError(
            f"Cannot map {len(psm_pairs)} PSM wing stations to "
            f"{len(fem_endpoints)} FEM stations."
        )

    mapping: dict[int, np.ndarray] = {}
    for (psm_forward, psm_rearward), fem_endpoint in zip(psm_pairs, fem_endpoints):
        mapping[psm_rearward] = fem_endpoint.rearward_coord
        if psm_forward in PSM_FORWARD_TIP_NODE_IDS:
            mapped_forward = psm_coords[psm_forward].copy()
            mapped_forward[1:] = fem_endpoint.rearward_coord[1:]
            mapping[psm_forward] = mapped_forward
        else:
            mapping[psm_forward] = fem_endpoint.forward_coord
    return mapping


def split_line_ending(line: str) -> tuple[str, str]:
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n"):
        return line[:-1], "\n"
    return line, ""


def format_yaml_float(value: float) -> str:
    value = float(value)
    if abs(value) < 5e-13:
        value = 0.0
    text = f"{value:.8f}".rstrip("0").rstrip(".")
    if text in {"", "-0"}:
        return "0"
    return text


def format_yaml_position(coord: np.ndarray) -> str:
    return ", ".join(format_yaml_float(value) for value in coord)


def write_psm_copy_with_mapped_wing_nodes(
    input_path: Path,
    output_path: Path,
    mapping: dict[int, np.ndarray],
) -> set[int]:
    node_row_re = re.compile(r"^(\s*-\s*)\[\s*([+-]?\d+)\s*,[^\]]*\](.*)$")
    top_level_key_re = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\s*:")

    active_wing_table = False
    updated_nodes: set[int] = set()
    output_lines: list[str] = []

    for line in input_path.read_text(encoding="utf-8").splitlines(keepends=True):
        body, ending = split_line_ending(line)
        top_level_key = top_level_key_re.match(body)
        if top_level_key is not None:
            active_wing_table = top_level_key.group(1) == "wing_particles"

        if active_wing_table:
            node_match = node_row_re.match(body)
            if node_match is not None:
                node_id = int(node_match.group(2))
                if node_id in mapping:
                    body = (
                        f"{node_match.group(1)}"
                        f"[{node_id}, {format_yaml_position(mapping[node_id])}]"
                        f"{node_match.group(3)}"
                    )
                    updated_nodes.add(node_id)

        output_lines.append(body + ending)

    missing_nodes = sorted(set(mapping) - updated_nodes)
    if missing_nodes:
        raise ValueError(f"Did not update mapped PSM wing nodes: {missing_nodes}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(output_lines), encoding="utf-8")
    return updated_nodes


def main() -> None:
    fem_geometry = load_yaml(GEOMETRY_PATH)
    psm_geometry = load_yaml(PSM_GEOMETRY_PATH)
    fem_endpoints = fem_station_endpoints(fem_geometry)
    mapping = build_psm_to_fem_mapping(psm_geometry, fem_endpoints)
    updated_nodes = write_psm_copy_with_mapped_wing_nodes(
        PSM_GEOMETRY_PATH,
        OUTPUT_PATH,
        mapping,
    )

    print(f"[map_psm_wing_nodes_to_fem_wing_nodes] wrote {OUTPUT_PATH}")
    print(
        "[map_psm_wing_nodes_to_fem_wing_nodes] mapped PSM wing nodes: "
        f"{sorted(updated_nodes)}"
    )
    print(
        "[map_psm_wing_nodes_to_fem_wing_nodes] kept x for PSM forward tip nodes: "
        f"{sorted(PSM_FORWARD_TIP_NODE_IDS)}"
    )


if __name__ == "__main__":
    main()
