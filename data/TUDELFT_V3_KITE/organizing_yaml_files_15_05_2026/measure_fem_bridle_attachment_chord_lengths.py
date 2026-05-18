"""
Measure chordwise distances between bridle attachment nodes on the FEM wing.

The script writes a markdown report in results/ and can show a small 3D
inspection plot with node labels when called with --show.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
GEOMETRY_PATH = SCRIPT_DIR / "results" / "struc_geometry_FEM_full_stretched.yaml"
PSM_GEOMETRY_PATH = (
    SCRIPT_DIR
    / "results"
    / "struc_geometry_PSM_reduced_stretched_mapped_to_FEM_stretched.yaml"
)
RESULTS_DIR = SCRIPT_DIR / "results"
OUTPUT_MD = RESULTS_DIR / "fem_bridle_attachment_chord_lengths.md"


@dataclass
class AttachmentNode:
    node_id: int
    coord: np.ndarray
    table: str
    connection_names: list[str]


@dataclass
class StationResult:
    name: str
    sort_index: int
    y_mean: float
    attachments: list[AttachmentNode]

    @property
    def sorted_attachments(self) -> list[AttachmentNode]:
        return sorted(self.attachments, key=lambda node: (node.coord[0], node.node_id))

    @property
    def forward_node(self) -> AttachmentNode | None:
        attachments = self.sorted_attachments
        return attachments[0] if attachments else None

    @property
    def rearward_node(self) -> AttachmentNode | None:
        attachments = self.sorted_attachments
        return attachments[-1] if attachments else None

    @property
    def distance(self) -> float | None:
        forward = self.forward_node
        rearward = self.rearward_node
        if forward is None or rearward is None or forward.node_id == rearward.node_id:
            return None
        return float(np.linalg.norm(rearward.coord - forward.coord))


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def load_wing_nodes(path: Path) -> dict[int, np.ndarray]:
    geometry = load_yaml(path)
    return {
        int(row[0]): np.array([float(row[1]), float(row[2]), float(row[3])])
        for row in geometry.get("wing_particles", {}).get("data", [])
    }


def build_coordinate_map(
    geometry: dict[str, Any],
) -> tuple[dict[int, np.ndarray], dict[int, str]]:
    coords: dict[int, np.ndarray] = {}
    tables: dict[int, str] = {}

    for table_name in ("wing_particles", "bridle_particles"):
        for row in geometry.get(table_name, {}).get("data", []):
            node_id = int(row[0])
            coords[node_id] = np.array(
                [float(row[1]), float(row[2]), float(row[3])],
                dtype=float,
            )
            tables[node_id] = table_name

    return coords, tables


def build_connection_lookup(
    geometry: dict[str, Any],
) -> tuple[set[int], dict[int, list[str]]]:
    connection_nodes: set[int] = set()
    connection_names_by_node: dict[int, list[str]] = {}

    for row in geometry.get("bridle_connections", {}).get("data", []):
        name = str(row[0])
        for value in row[1:]:
            node_id = int(value)
            if node_id == 0:
                continue
            connection_nodes.add(node_id)
            connection_names_by_node.setdefault(node_id, []).append(name)

    return connection_nodes, connection_names_by_node


def direct_wing_particle_attachments(
    geometry: dict[str, Any],
    connection_nodes: set[int],
) -> set[int]:
    wing_node_ids = {
        int(row[0]) for row in geometry.get("wing_particles", {}).get("data", [])
    }
    return wing_node_ids.intersection(connection_nodes)


def build_station_results(geometry: dict[str, Any]) -> list[StationResult]:
    coords, tables = build_coordinate_map(geometry)
    connection_nodes, connection_names_by_node = build_connection_lookup(geometry)
    assigned_nodes: set[int] = set()
    stations: list[StationResult] = []

    for station_i, row in enumerate(geometry.get("strut_tubes", {}).get("data", [])):
        name = str(row[0])
        station_node_ids = [int(node_id) for node_id in row[6]]
        station_coords = [
            coords[node_id] for node_id in station_node_ids if node_id in coords
        ]
        y_mean = float(np.mean([coord[1] for coord in station_coords]))

        attachment_ids = [
            node_id for node_id in station_node_ids if node_id in connection_nodes
        ]
        assigned_nodes.update(attachment_ids)

        stations.append(
            StationResult(
                name=name,
                sort_index=station_i,
                y_mean=y_mean,
                attachments=[
                    AttachmentNode(
                        node_id=node_id,
                        coord=coords[node_id],
                        table=tables[node_id],
                        connection_names=connection_names_by_node.get(node_id, []),
                    )
                    for node_id in attachment_ids
                ],
            )
        )

    # A few bridle connections attach directly to wing_particles instead of nodes
    # listed in strut_tubes.node_indices. Assign them to the closest station in y.
    for node_id in sorted(direct_wing_particle_attachments(geometry, connection_nodes)):
        if node_id in assigned_nodes:
            continue
        coord = coords[node_id]
        station = min(stations, key=lambda item: abs(coord[1] - item.y_mean))
        station.attachments.append(
            AttachmentNode(
                node_id=node_id,
                coord=coord,
                table=tables[node_id],
                connection_names=connection_names_by_node.get(node_id, []),
            )
        )
        assigned_nodes.add(node_id)

    return stations


def format_coord(coord: np.ndarray) -> str:
    return f"[{coord[0]:.6f}, {coord[1]:.6f}, {coord[2]:.6f}]"


def format_node(node: AttachmentNode | None) -> str:
    if node is None:
        return "-"
    return f"{node.node_id} {format_coord(node.coord)}"


def write_markdown_report(
    stations: list[StationResult],
    *,
    geometry_path: Path,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# FEM Bridle Attachment Chord Lengths",
        "",
        f"Input geometry: `{geometry_path.name}`",
        "",
        "Attachment nodes are nodes that appear in `bridle_connections` and also "
        "belong to a `strut_tubes.node_indices` station. Direct `wing_particles` "
        "attachments are assigned to the nearest station by y-coordinate.",
        "",
        "| Station | Attachment nodes sorted by x | Most forward node [x,y,z] | Most rearward node [x,y,z] | Distance [m] |",
        "| --- | --- | --- | --- | ---: |",
    ]

    for station in stations:
        attachments = station.sorted_attachments
        attachment_label = ", ".join(str(node.node_id) for node in attachments)
        distance = station.distance
        distance_text = "-" if distance is None else f"{distance:.6f}"
        lines.append(
            "| "
            f"{station.name} | "
            f"{attachment_label or '-'} | "
            f"{format_node(station.forward_node)} | "
            f"{format_node(station.rearward_node)} | "
            f"{distance_text} |"
        )

    lines.extend(
        [
            "",
            "## Attachment Details",
            "",
            "| Station | Node | Source table | Position [x,y,z] | Bridle connections |",
            "| --- | ---: | --- | --- | --- |",
        ]
    )

    for station in stations:
        for node in station.sorted_attachments:
            connection_label = ", ".join(node.connection_names)
            lines.append(
                "| "
                f"{station.name} | "
                f"{node.node_id} | "
                f"{node.table} | "
                f"{format_coord(node.coord)} | "
                f"{connection_label} |"
            )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[measure_fem_bridle_attachment_chord_lengths] wrote {output_path}")


def plot_attachment_nodes(
    stations: list[StationResult],
    psm_wing_nodes: dict[int, np.ndarray],
) -> None:
    fig = plt.figure(figsize=(9.0, 6.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("FEM bridle attachment nodes and PSM wing nodes")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")

    fem_label_added = False
    for station in stations:
        attachments = station.sorted_attachments
        if not attachments:
            continue
        points = np.array([node.coord for node in attachments])
        ax.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            s=30,
            color="tab:blue",
            label="FEM bridle attachment" if not fem_label_added else None,
        )
        fem_label_added = True
        for node in attachments:
            ax.text(
                node.coord[0],
                node.coord[1],
                node.coord[2],
                f"FEM {node.node_id}",
                fontsize=8,
                color="tab:blue",
            )

        forward = station.forward_node
        rearward = station.rearward_node
        if (
            forward is not None
            and rearward is not None
            and forward.node_id != rearward.node_id
        ):
            segment = np.array([forward.coord, rearward.coord])
            ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], color="black", lw=1.0)

    if psm_wing_nodes:
        psm_ids = sorted(psm_wing_nodes)
        psm_points = np.array([psm_wing_nodes[node_id] for node_id in psm_ids])
        ax.scatter(
            psm_points[:, 0],
            psm_points[:, 1],
            psm_points[:, 2],
            s=22,
            color="tab:orange",
            marker="^",
            label="PSM wing node",
        )
        for node_id in psm_ids:
            coord = psm_wing_nodes[node_id]
            ax.text(
                coord[0],
                coord[1],
                coord[2],
                f"PSM {node_id}",
                fontsize=8,
                color="tab:orange",
            )

    all_coords = [node.coord for station in stations for node in station.attachments]
    all_coords.extend(psm_wing_nodes.values())
    all_points = np.array(all_coords)
    if len(all_points) > 0:
        center = 0.5 * (all_points.min(axis=0) + all_points.max(axis=0))
        half_range = 0.55 * float(
            np.max(all_points.max(axis=0) - all_points.min(axis=0))
        )
        ax.set_xlim(center[0] - half_range, center[0] + half_range)
        ax.set_ylim(center[1] - half_range, center[1] + half_range)
        ax.set_zlim(center[2] - half_range, center[2] + half_range)
        ax.set_box_aspect((1.0, 1.0, 1.0))

    ax.view_init(elev=20, azim=-120)
    ax.legend(loc="upper right")
    fig.tight_layout()
    plt.show()


def main(show: bool = False) -> None:
    geometry = load_yaml(GEOMETRY_PATH)
    stations = build_station_results(geometry)
    write_markdown_report(
        stations,
        geometry_path=GEOMETRY_PATH,
        output_path=OUTPUT_MD,
    )
    if show:
        plot_attachment_nodes(stations, load_wing_nodes(PSM_GEOMETRY_PATH))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--show", action="store_true", help="show the 3D node plot")
    args = parser.parse_args()
    main(show=args.show)
