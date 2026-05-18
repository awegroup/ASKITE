#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_WING_MASS_KG = 11.0
REFERENCE_FULL_BRIDLE_MASS_KG = 2.6


@dataclass
class MassSummary:
    path: Path
    model_kind: str
    wing_mass_kg: float
    wing_mass_source: str
    bridle_line_mass_kg: float
    pulley_mass_kg: float
    pulley_count: int
    bridle_connection_count: int
    bridle_volume_m3: float
    kite_mass_kg: float
    kcu_mass_kg: float
    dyneema_density_kg_m3: float | None
    required_dyneema_density_for_bridle_target_kg_m3: float | None
    wing_mass_scale_to_reference: float | None
    notes: list[str]
    wing_group_masses: dict[str, float]

    @property
    def bridle_mass_kg(self) -> float:
        return self.bridle_line_mass_kg + self.pulley_mass_kg

    @property
    def kite_plus_kcu_mass_kg(self) -> float:
        return self.kite_mass_kg + self.kcu_mass_kg


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError(f"{path} did not parse as a YAML mapping")
    return cfg


def table_rows(cfg: dict[str, Any], section: str) -> list[dict[str, Any]]:
    table = cfg.get(section)
    if not isinstance(table, dict):
        return []
    headers = [str(header) for header in table.get("headers", [])]
    rows = []
    for row in table.get("data", []):
        if not isinstance(row, list):
            continue
        item = dict(zip(headers, row))
        item["_raw"] = row
        rows.append(item)
    return rows


def table_by_name(cfg: dict[str, Any], section: str) -> dict[str, dict[str, Any]]:
    return {str(row["name"]): row for row in table_rows(cfg, section) if "name" in row}


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def material_density(cfg: dict[str, Any], material: Any) -> float:
    material_cfg = cfg.get(str(material), {})
    if isinstance(material_cfg, dict):
        density = as_float(material_cfg.get("density"))
        if density > 0.0:
            return density
    return 970.0


def line_length(line: dict[str, Any]) -> float:
    return as_float(line.get("rest_length", line.get("l0")))


def line_diameter_m(line: dict[str, Any]) -> float:
    diameter = as_float(line.get("diameter", line.get("d")))
    # Be tolerant of old files that store bridle diameters in mm.
    if diameter > 0.05:
        diameter /= 1000.0
    return diameter


def line_density(cfg: dict[str, Any], line: dict[str, Any]) -> float:
    explicit_density = as_float(line.get("density"))
    if explicit_density > 0.0:
        return explicit_density
    return material_density(cfg, line.get("material", "dyneema"))


def line_cross_section_area_m2(line: dict[str, Any]) -> float:
    return math.pi * (line_diameter_m(line) / 2.0) ** 2


def line_mass_kg(cfg: dict[str, Any], line: dict[str, Any]) -> float:
    return line_density(cfg, line) * line_cross_section_area_m2(line) * line_length(line)


def is_pulley_connection(conn: dict[str, Any]) -> bool:
    raw = conn.get("_raw", [])
    return isinstance(raw, list) and len(raw) >= 4 and raw[3] not in (None, "", 0)


def wing_mass(cfg: dict[str, Any]) -> tuple[float, str, dict[str, float]]:
    if "wing_elements" in cfg:
        elements = table_by_name(cfg, "wing_elements")
        group_masses: dict[str, float] = defaultdict(float)
        total = 0.0
        for conn in table_rows(cfg, "wing_connections"):
            name = str(conn.get("name", ""))
            element = elements.get(name)
            if element is None:
                continue
            mass = as_float(element.get("m"))
            total += mass
            group_masses[name.split("_", 1)[0]] += mass
        return total, "sum wing_elements.m over wing_connections", dict(group_masses)

    if "mass_without_bridles" in cfg:
        return (
            as_float(cfg["mass_without_bridles"]),
            "mass_without_bridles",
            {},
        )

    return 0.0, "not found", {}


def bridle_mass(cfg: dict[str, Any]) -> tuple[float, float, int, int, float, list[str]]:
    element_section = "bridle_lines" if "bridle_lines" in cfg else "bridle_elements"
    elements = table_by_name(cfg, element_section)
    line_total = 0.0
    volume_total = 0.0
    pulley_count = 0
    notes: list[str] = []
    connection_keys: list[tuple[Any, ...]] = []

    for conn in table_rows(cfg, "bridle_connections"):
        raw = conn.get("_raw", [])
        if isinstance(raw, list):
            connection_keys.append(tuple(raw))
        line = elements.get(str(conn.get("name", "")))
        if line is None:
            notes.append(f"missing bridle line definition for {conn.get('name')!r}")
            continue
        line_total += line_mass_kg(cfg, line)
        volume_total += line_cross_section_area_m2(line) * line_length(line)
        if is_pulley_connection(conn):
            pulley_count += 1

    duplicate_rows = [
        (row, count) for row, count in Counter(connection_keys).items() if count > 1
    ]
    if duplicate_rows:
        duplicate_text = ", ".join(
            f"{list(row)} x{count}" for row, count in sorted(duplicate_rows)
        )
        notes.append(f"duplicate bridle connection rows: {duplicate_text}")

    pulley_total = pulley_count * as_float(cfg.get("pulley_mass"))
    return (
        line_total,
        pulley_total,
        pulley_count,
        len(connection_keys),
        volume_total,
        notes,
    )


def summarize(path: Path) -> MassSummary:
    cfg = load_yaml(path)
    notes: list[str] = []

    wing_total, wing_source, wing_group_masses = wing_mass(cfg)
    (
        bridle_line_total,
        pulley_total,
        pulley_count,
        bridle_connection_count,
        bridle_volume,
        bridle_notes,
    ) = bridle_mass(cfg)
    notes.extend(bridle_notes)

    dyneema_density = None
    if isinstance(cfg.get("dyneema"), dict):
        dyneema_density = as_float(cfg["dyneema"].get("density"))

    if "wing_elements" in cfg:
        model_kind = "PSM reduced"
    elif "strut_tubes" in cfg or "leading_edge_tubes" in cfg:
        model_kind = "FEM full"
    else:
        model_kind = "unknown"

    if abs(wing_total - REFERENCE_WING_MASS_KG) > 0.1:
        notes.append(
            f"wing mass differs from 11.0 kg by {wing_total - REFERENCE_WING_MASS_KG:+.3f} kg"
        )
    bridle_total = bridle_line_total + pulley_total
    if abs(bridle_total - REFERENCE_FULL_BRIDLE_MASS_KG) > 0.1:
        notes.append(
            "bridle mass differs from full-flight 2.6 kg by "
            f"{bridle_total - REFERENCE_FULL_BRIDLE_MASS_KG:+.3f} kg"
        )

    required_density = None
    if dyneema_density and bridle_line_total > 0.0:
        target_line_mass = REFERENCE_FULL_BRIDLE_MASS_KG - pulley_total
        if target_line_mass > 0.0:
            required_density = dyneema_density * target_line_mass / bridle_line_total

    wing_mass_scale = None
    if "wing_elements" in cfg and wing_total > 0.0:
        wing_mass_scale = REFERENCE_WING_MASS_KG / wing_total

    return MassSummary(
        path=path,
        model_kind=model_kind,
        wing_mass_kg=wing_total,
        wing_mass_source=wing_source,
        bridle_line_mass_kg=bridle_line_total,
        pulley_mass_kg=pulley_total,
        pulley_count=pulley_count,
        bridle_connection_count=bridle_connection_count,
        bridle_volume_m3=bridle_volume,
        kite_mass_kg=wing_total + bridle_total,
        kcu_mass_kg=as_float(cfg.get("kcu_mass")),
        dyneema_density_kg_m3=dyneema_density,
        required_dyneema_density_for_bridle_target_kg_m3=required_density,
        wing_mass_scale_to_reference=wing_mass_scale,
        notes=notes,
        wing_group_masses=wing_group_masses,
    )


def fmt(value: float) -> str:
    return f"{value:.3f}"


def print_markdown_table(summaries: list[MassSummary]) -> None:
    print(
        "| file | kind | wing kg | bridle kg | line kg | pulley kg | pulleys | "
        "kite kg | kite+KCU kg | rho Dyneema | rho for 2.6 kg bridle | wing source |"
    )
    print(
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"
    )
    for item in summaries:
        density = (
            f"{item.dyneema_density_kg_m3:.0f}"
            if item.dyneema_density_kg_m3 is not None
            else ""
        )
        required_density = (
            f"{item.required_dyneema_density_for_bridle_target_kg_m3:.0f}"
            if item.required_dyneema_density_for_bridle_target_kg_m3 is not None
            else ""
        )
        print(
            f"| {item.path.name} | {item.model_kind} | {fmt(item.wing_mass_kg)} | "
            f"{fmt(item.bridle_mass_kg)} | {fmt(item.bridle_line_mass_kg)} | "
            f"{fmt(item.pulley_mass_kg)} | {item.pulley_count} | "
            f"{fmt(item.kite_mass_kg)} | {fmt(item.kite_plus_kcu_mass_kg)} | "
            f"{density} | {required_density} | {item.wing_mass_source} |"
        )


def print_details(summaries: list[MassSummary]) -> None:
    print()
    print("Reference values used:")
    print(f"- wing mass without bridles: {REFERENCE_WING_MASS_KG:.3f} kg")
    print(f"- full bridle system mass: {REFERENCE_FULL_BRIDLE_MASS_KG:.3f} kg")
    print(
        f"- full kite mass, excluding KCU: "
        f"{REFERENCE_WING_MASS_KG + REFERENCE_FULL_BRIDLE_MASS_KG:.3f} kg"
    )
    print()
    print("Details:")
    for item in summaries:
        print(f"- {item.path.name}")
        print(f"  - bridle connections counted: {item.bridle_connection_count}")
        print(f"  - bridle line volume: {item.bridle_volume_m3:.9f} m^3")
        if item.wing_group_masses:
            groups = ", ".join(
                f"{name}={mass:.3f} kg"
                for name, mass in sorted(item.wing_group_masses.items())
            )
            print(f"  - PSM wing mass by connection group: {groups}")
        if item.wing_mass_scale_to_reference is not None:
            print(
                "  - uniform PSM wing-element mass scale for 11.0 kg: "
                f"{item.wing_mass_scale_to_reference:.6f}"
            )
        for note in item.notes:
            print(f"  - note: {note}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate full-kite mass totals carried by the current "
            "struc_geometry YAML files."
        )
    )
    parser.add_argument(
        "yaml_files",
        nargs="*",
        type=Path,
        help="YAML files to inspect. Defaults to struc_geometry_*.yaml in this folder.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = args.yaml_files or sorted(SCRIPT_DIR.glob("struc_geometry_*.yaml"))
    summaries = [summarize(path) for path in paths]
    print_markdown_table(summaries)
    print_details(summaries)


if __name__ == "__main__":
    main()
