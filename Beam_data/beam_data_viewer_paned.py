#!/usr/bin/env python3
"""
Beam data viewer: MCC + ASC + dose-engine CSV overlay UI

Loads measured water-tank .mcc/.asc files and dose-engine .csv files, auto-detects
PDD vs profile curves, applies simple normalization, and overlays PDDs and
profiles in separate tabs. Curves can be toggled on/off from the left table.

Requirements:
    python3 -m pip install numpy matplotlib

Usage:
    python3 beam_data_viewer.py
    python3 beam_data_viewer.py file1.mcc file2.asc dose_cax.csv dose_crossline_50mm.csv
    python3 beam_data_viewer.py --summary *.mcc *.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# GUI imports are intentionally kept after core imports so --summary can fail
# gracefully on systems without a display.


@dataclass
class Curve:
    path: Path
    name: str
    kind: str  # "PDD" or "PROFILE"
    x: np.ndarray
    y: np.ndarray
    metadata: Dict[str, str] = field(default_factory=dict)
    source: str = ""
    x_label: str = "Position / depth (mm)"
    y_label: str = "Signal"

    def label(self) -> str:
        bits = []
        energy = self.metadata.get("ENERGY", "")
        filt = self.metadata.get("FILTER", "")
        field_ip = self.metadata.get("FIELD_INPLANE", "")
        field_cp = self.metadata.get("FIELD_CROSSPLANE", "")
        depth = self.metadata.get("SCAN_DEPTH", "")
        curve_type = self.metadata.get("SCAN_CURVETYPE", self.kind)

        if energy:
            # ENERGY in MCC is usually numeric MV.
            try:
                bits.append(f"{float(energy):g} MV")
            except ValueError:
                bits.append(str(energy))
        if filt:
            bits.append(filt)
        if field_ip and field_cp:
            try:
                bits.append(f"{float(field_ip):g}x{float(field_cp):g} mm")
            except ValueError:
                bits.append(f"{field_ip}x{field_cp} mm")
        if depth and self.kind == "PROFILE":
            try:
                bits.append(f"d={float(depth):g} mm")
            except ValueError:
                bits.append(f"d={depth} mm")
        if curve_type and curve_type.upper() not in {"PDD", "PROFILE"}:
            bits.append(curve_type.replace("_", " ").title())

        prefix = " | ".join(bits)
        return f"{self.name}" + (f" ({prefix})" if prefix else "")


def _parse_float(text: str) -> Optional[float]:
    try:
        return float(text.strip())
    except Exception:
        return None


def _safe_divide(numer: np.ndarray, denom: np.ndarray) -> np.ndarray:
    out = np.full_like(numer, np.nan, dtype=float)
    mask = np.isfinite(numer) & np.isfinite(denom) & (denom != 0)
    out[mask] = numer[mask] / denom[mask]
    return out


def _sort_by_x(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    order = np.argsort(x)
    return x[order], y[order]


def parse_mcc(path: Path) -> List[Curve]:
    """Parse PTW/IBA-style MCC exports with one or more BEGIN_SCAN blocks."""
    curves: List[Curve] = []
    lines = path.read_text(errors="replace").splitlines()

    i = 0
    scan_counter = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("BEGIN_SCAN"):
            i += 1
            continue

        scan_counter += 1
        metadata: Dict[str, str] = {}
        data_rows: List[Tuple[float, float, float]] = []
        i += 1
        in_data = False

        while i < len(lines):
            raw = lines[i]
            stripped = raw.strip()

            if stripped.startswith("END_SCAN"):
                break

            if stripped == "BEGIN_DATA":
                in_data = True
                i += 1
                continue

            if stripped == "END_DATA":
                in_data = False
                i += 1
                continue

            if in_data:
                parts = stripped.split()
                if len(parts) >= 2:
                    vals = [_parse_float(p) for p in parts[:3]]
                    if vals[0] is not None and vals[1] is not None:
                        # MCC commonly stores: position, field detector, ref detector.
                        # If reference is missing, use 1.0.
                        data_rows.append((vals[0], vals[1], vals[2] if vals[2] is not None else 1.0))
            else:
                if "=" in stripped:
                    key, value = stripped.split("=", 1)
                    metadata[key.strip()] = value.strip()
            i += 1

        if data_rows:
            arr = np.array(data_rows, dtype=float)
            x = arr[:, 0]
            field_signal = arr[:, 1]
            ref_signal = arr[:, 2]
            y = _safe_divide(field_signal, ref_signal)
            x, y = _sort_by_x(x, y)

            curve_type = metadata.get("SCAN_CURVETYPE", "").upper()
            if "PDD" in curve_type:
                kind = "PDD"
                x_label = "Depth (mm)"
            else:
                kind = "PROFILE"
                x_label = "Off-axis position (mm)"

            pretty_type = metadata.get("SCAN_CURVETYPE", kind).replace("_", " ").title()
            name = f"{path.stem} | scan {scan_counter}: {pretty_type}"
            curves.append(
                Curve(
                    path=path,
                    name=name,
                    kind=kind,
                    x=x,
                    y=y,
                    metadata=metadata,
                    source="MCC",
                    x_label=x_label,
                    y_label="Field / reference signal",
                )
            )

        i += 1

    return curves



def _parse_energy_from_filename(path: Path) -> str:
    """Best-effort energy extraction from names such as X06, X10 FFF, 6MV."""
    stem = path.stem.upper()
    m = re.search(r"\bX\s*0?(\d{1,2})(?:\s|_|-|$)", stem)
    if not m:
        m = re.search(r"\b0?(\d{1,2})\s*MV\b", stem)
    return m.group(1) if m else ""


def parse_asc(path: Path) -> List[Curve]:
    """Parse simple water-tank .ASC exports containing $STOM/$ENOM blocks.

    This format commonly stores one measurement block per scan:
        %TYPE OPD / %AXIS Z  -> PDD
        %TYPE OPP / %AXIS X/Y -> profile

    Data rows are expected as <x y z dose>, with coordinates in mm and dose
    already usually normalized/smoothed by the water-tank software.
    """
    curves: List[Curve] = []
    lines = path.read_text(errors="replace").splitlines()

    scan_counter = 0
    i = 0
    while i < len(lines):
        if not lines[i].strip().startswith("$STOM"):
            i += 1
            continue

        scan_counter += 1
        i += 1
        metadata: Dict[str, str] = {}
        data_rows: List[Tuple[float, float, float, float]] = []

        while i < len(lines):
            stripped = lines[i].strip()
            if stripped.startswith("$ENOM") or stripped.startswith("$ENOF"):
                break
            if stripped.startswith("$") and not stripped.startswith("$STOM"):
                # Unknown block marker; leave this scan cleanly.
                break

            if stripped.startswith("%"):
                parts = stripped[1:].split(None, 1)
                if parts:
                    key = parts[0].strip().upper()
                    value = parts[1].strip() if len(parts) > 1 else ""
                    metadata[key] = value
            elif stripped.startswith("<") and stripped.endswith(">"):
                vals = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", stripped)
                if len(vals) >= 4:
                    try:
                        data_rows.append(tuple(float(v) for v in vals[:4]))  # type: ignore[arg-type]
                    except ValueError:
                        pass
            i += 1

        if data_rows:
            arr = np.array(data_rows, dtype=float)
            axis = metadata.get("AXIS", "").upper()
            asc_type = metadata.get("TYPE", "").upper()

            if axis == "Z" or asc_type == "OPD":
                kind = "PDD"
                x = arr[:, 2]
                x_label = "Depth (mm)"
                curve_type = "PDD"
            else:
                kind = "PROFILE"
                if axis == "Y":
                    x = arr[:, 1]
                    curve_type = "Inplane profile" if axis == "Y" else "Profile"
                else:
                    x = arr[:, 0]
                    curve_type = "Crossplane profile" if axis == "X" else "Profile"
                x_label = "Off-axis position (mm)"

            y = arr[:, 3]
            x, y = _sort_by_x(x.astype(float), y.astype(float))

            # Map useful ASC metadata into the same keys used by MCC/CSV labels.
            common_metadata: Dict[str, str] = {f"ASC_{k}": v for k, v in metadata.items()}
            common_metadata["SCAN_CURVETYPE"] = curve_type
            common_metadata["ASC_TYPE"] = asc_type
            common_metadata["ASC_AXIS"] = axis
            if metadata.get("DPTH") and kind == "PROFILE":
                common_metadata["SCAN_DEPTH"] = metadata["DPTH"]
            if metadata.get("SSD"):
                common_metadata["SSD"] = metadata["SSD"]
            if metadata.get("FLSZ"):
                fmatch = re.match(r"\s*([0-9.]+)\s*[*xX]\s*([0-9.]+)", metadata["FLSZ"])
                if fmatch:
                    common_metadata["FIELD_INPLANE"] = fmatch.group(1)
                    common_metadata["FIELD_CROSSPLANE"] = fmatch.group(2)
                common_metadata["FIELD_SIZE"] = metadata["FLSZ"]
            energy = _parse_energy_from_filename(path)
            if energy:
                common_metadata["ENERGY"] = energy
            if "FFF" in path.stem.upper():
                common_metadata["FILTER"] = "FFF"

            name = f"{path.stem} | scan {scan_counter}: {curve_type}"
            curves.append(
                Curve(
                    path=path,
                    name=name,
                    kind=kind,
                    x=x,
                    y=y,
                    metadata=common_metadata,
                    source="ASC",
                    x_label=x_label,
                    y_label="Dose (%) / exported value",
                )
            )

        i += 1

    return curves

def _numeric_columns(rows: List[Dict[str, str]], fieldnames: Sequence[str]) -> List[str]:
    numeric: List[str] = []
    for col in fieldnames:
        ok = 0
        for row in rows[:20]:
            if _parse_float(row.get(col, "")) is not None:
                ok += 1
        if ok >= max(1, min(5, len(rows[:20]))):
            numeric.append(col)
    return numeric


def _guess_depth_from_filename(path: Path) -> str:
    # Examples: dose_crossline_50mm.csv, profile_depth100.csv
    m = re.search(r"(?:depth[_ -]?)?(\d+(?:\.\d+)?)\s*mm", path.stem, flags=re.IGNORECASE)
    return m.group(1) if m else ""


def parse_csv(path: Path) -> List[Curve]:
    """Parse simple dose-engine CSV files.

    Auto-detects PDD when a z/depth column exists; otherwise profile when an
    x/y/off-axis column exists. Dose is taken from 'dose' if present, otherwise
    the last numeric non-axis column.
    """
    curves: List[Curve] = []
    with path.open("r", newline="", errors="replace") as f:
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t ")
        reader = csv.DictReader(f, dialect=dialect)
        if reader.fieldnames is None:
            return curves
        rows = list(reader)

    if not rows:
        return curves

    fieldnames = [fn.strip() for fn in reader.fieldnames]
    # Normalize row keys if DictReader kept spaces.
    normalized_rows: List[Dict[str, str]] = []
    for row in rows:
        normalized_rows.append({(k.strip() if k else ""): v for k, v in row.items()})
    rows = normalized_rows

    lower_map = {c.lower(): c for c in fieldnames}

    axis_col = None
    kind = None
    x_label = "Position / depth (mm)"

    for candidate in ["z_mm", "depth_mm", "depth", "z"]:
        if candidate in lower_map:
            axis_col = lower_map[candidate]
            kind = "PDD"
            x_label = "Depth (mm)"
            break

    if axis_col is None:
        for candidate in ["x_mm", "y_mm", "position_mm", "pos_mm", "offaxis_mm", "x", "y"]:
            if candidate in lower_map:
                axis_col = lower_map[candidate]
                kind = "PROFILE"
                x_label = "Off-axis position (mm)"
                break

    if axis_col is None:
        # Fall back to the first numeric column as x.
        numeric = _numeric_columns(rows, fieldnames)
        if len(numeric) < 2:
            return curves
        axis_col = numeric[0]
        kind = "PROFILE"
        x_label = axis_col

    dose_col = None
    for candidate in ["dose", "signal", "value", "reading", "terma"]:
        if candidate in lower_map and lower_map[candidate] != axis_col:
            dose_col = lower_map[candidate]
            break

    if dose_col is None:
        numeric = [c for c in _numeric_columns(rows, fieldnames) if c != axis_col]
        if not numeric:
            return curves
        dose_col = numeric[-1]

    x_vals: List[float] = []
    y_vals: List[float] = []
    for row in rows:
        x = _parse_float(row.get(axis_col, ""))
        y = _parse_float(row.get(dose_col, ""))
        if x is not None and y is not None:
            x_vals.append(x)
            y_vals.append(y)

    if not x_vals:
        return curves

    x = np.array(x_vals, dtype=float)
    y = np.array(y_vals, dtype=float)
    x, y = _sort_by_x(x, y)

    metadata = {
        "CSV_AXIS_COLUMN": axis_col,
        "CSV_VALUE_COLUMN": dose_col,
        "SCAN_CURVETYPE": kind or "PROFILE",
    }
    if kind == "PROFILE":
        depth = _guess_depth_from_filename(path)
        if depth:
            metadata["SCAN_DEPTH"] = depth

    curves.append(
        Curve(
            path=path,
            name=path.stem,
            kind=kind or "PROFILE",
            x=x,
            y=y,
            metadata=metadata,
            source="CSV",
            x_label=x_label,
            y_label=dose_col,
        )
    )
    return curves


def load_curves(paths: Iterable[Path]) -> List[Curve]:
    curves: List[Curve] = []
    for p in paths:
        p = Path(p)
        if not p.exists() or not p.is_file():
            continue
        suffix = p.suffix.lower()
        try:
            if suffix == ".mcc":
                curves.extend(parse_mcc(p))
            elif suffix in {".asc", ".prs"}:
                curves.extend(parse_asc(p))
            elif suffix == ".csv":
                curves.extend(parse_csv(p))
        except Exception as exc:
            print(f"Could not parse {p}: {exc}", file=sys.stderr)
    return curves


def normalize_y(curve: Curve, mode: str) -> np.ndarray:
    """Return normalized y according to UI mode."""
    y = curve.y.astype(float).copy()
    if mode == "Raw":
        return y

    if curve.kind == "PDD":
        # PDD convention: normalize to maximum dose = 100%.
        denom = np.nanmax(y) if np.any(np.isfinite(y)) else np.nan
        return y / denom * 100.0 if np.isfinite(denom) and denom != 0 else y

    if mode == "Max = 100%":
        denom = np.nanmax(y) if np.any(np.isfinite(y)) else np.nan
    elif mode == "CAX = 100%":
        # Interpolate at x=0 when possible; fallback to nearest point to zero.
        x = curve.x
        y_sorted = y
        if len(x) >= 2 and np.nanmin(x) <= 0 <= np.nanmax(x):
            denom = float(np.interp(0.0, x, y_sorted))
        else:
            denom = y_sorted[int(np.nanargmin(np.abs(x)))]
    else:
        denom = np.nanmax(y) if np.any(np.isfinite(y)) else np.nan

    return y / denom * 100.0 if np.isfinite(denom) and denom != 0 else y


def curve_summary(curves: Sequence[Curve]) -> str:
    lines = []
    for c in curves:
        x_min = float(np.nanmin(c.x)) if len(c.x) else math.nan
        x_max = float(np.nanmax(c.x)) if len(c.x) else math.nan
        depth = c.metadata.get("SCAN_DEPTH", "")
        energy = c.metadata.get("ENERGY", "")
        filt = c.metadata.get("FILTER", "")
        field_ip = c.metadata.get("FIELD_INPLANE", "")
        field_cp = c.metadata.get("FIELD_CROSSPLANE", "")
        bits = [c.kind, f"n={len(c.x)}", f"x=[{x_min:g}, {x_max:g}] mm"]
        if energy:
            bits.append(f"E={energy} MV")
        if filt:
            bits.append(f"filter={filt}")
        if field_ip and field_cp:
            bits.append(f"field={field_ip}x{field_cp} mm")
        if depth:
            bits.append(f"depth={depth} mm")
        lines.append(f"- {c.name}: " + ", ".join(bits))
    return "\n".join(lines)


def _first_existing(candidates: Iterable[Path]) -> Optional[Path]:
    """Return the first existing file from a candidate path list."""
    seen = set()
    for candidate in candidates:
        try:
            c = Path(candidate).expanduser()
            key = str(c.resolve(strict=False))
        except Exception:
            c = Path(candidate)
            key = str(c)
        if key in seen:
            continue
        seen.add(key)
        if c.exists() and c.is_file():
            return c
    return None


def default_autoload_paths() -> List[Path]:
    """Default files for Rohit's head_model dose-engine comparison.

    Expected launch location:
        .../head_model/dose_engine/Beam_data

    Reference measured file:
        FINAL_DATA/PROFILE_PDD/X06 OPEN 10X10 PDDXY SMOOTHED.ASC

    Dose-engine files:
        ../build/dose_crossline_200mm.csv
        ../build/dose_crossline_100mm.csv
        ../build/dose_cax.csv
    """
    cwd = Path.cwd()
    script_dir = Path(__file__).resolve().parent

    ref_rel = Path("FINAL_DATA") / "PROFILE_PDD" / "X06 OPEN 10X10 PDDXY SMOOTHED.ASC"
    dose_names = ["dose_crossline_200mm.csv", "dose_crossline_100mm.csv", "dose_cax.csv"]

    paths: List[Path] = []

    ref = _first_existing([
        cwd / ref_rel,
        cwd.parent / ref_rel,
        script_dir / ref_rel,
        script_dir.parent / ref_rel,
        # Useful when testing with files in the same directory as this script.
        cwd / ref_rel.name,
        script_dir / ref_rel.name,
    ])
    if ref is not None:
        paths.append(ref)

    for name in dose_names:
        found = _first_existing([
            cwd.parent / "build" / name,
            cwd / ".." / "build" / name,
            script_dir.parent / "build" / name,
            script_dir / ".." / "build" / name,
            # Useful when testing with files in the same directory as this script.
            cwd / name,
            script_dir / name,
        ])
        if found is not None:
            paths.append(found)

    return paths


def _depth_as_float_mm(value: str) -> Optional[float]:
    """Parse depth strings such as '015', '100', '100 mm'."""
    if not value:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?", str(value))
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


class BeamDataViewer:
    def __init__(self, initial_files: Optional[Sequence[Path]] = None, apply_compare_view: bool = False) -> None:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure

        self.tk = tk
        self.ttk = ttk
        self.filedialog = filedialog
        self.messagebox = messagebox
        self.FigureCanvasTkAgg = FigureCanvasTkAgg
        self.NavigationToolbar2Tk = NavigationToolbar2Tk
        self.Figure = Figure

        self.root = tk.Tk()
        self.root.title("Beam Data Viewer - MCC/ASC + Dose Engine CSV")
        self.root.geometry("1250x820")

        self.curves: List[Curve] = []
        self.curve_visible: List[bool] = []
        self.profile_norm_mode = tk.StringVar(value="CAX = 100%")
        self.pdd_norm_mode = tk.StringVar(value="Max = 100%")
        self.selected_only = tk.BooleanVar(value=False)
        self.show_grid = tk.BooleanVar(value=True)

        self._build_ui()

        if initial_files:
            self.add_files(initial_files)
            if apply_compare_view:
                self.apply_default_compare_visibility()

    def _build_ui(self) -> None:
        tk = self.tk
        ttk = self.ttk

        root = self.root
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        root.rowconfigure(1, weight=0)

        # Main layout: draggable horizontal split.
        # The left sidebar is initially placed at ~25% of the window width
        # by _set_initial_sash_position(), but the user can drag the sash.
        self.main_pane = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        self.main_pane.grid(row=0, column=0, sticky="nsew")

        left = ttk.Frame(self.main_pane, padding=8)
        left.rowconfigure(2, weight=1)
        left.columnconfigure(0, weight=1)

        right = ttk.Frame(self.main_pane, padding=(0, 8, 8, 8))
        right.rowconfigure(0, weight=1)
        right.columnconfigure(0, weight=1)

        # Weight 1:3 gives the split sensible behaviour as the window grows.
        self.main_pane.add(left, weight=1)
        self.main_pane.add(right, weight=3)

        buttons = ttk.Frame(left)
        buttons.grid(row=0, column=0, sticky="ew")
        ttk.Button(buttons, text="Load files", command=self.choose_files).grid(row=0, column=0, padx=(0, 4))
        ttk.Button(buttons, text="Auto-load", command=self.autoload_defaults).grid(row=0, column=1, padx=(0, 4))
        ttk.Button(buttons, text="Clear", command=self.clear).grid(row=0, column=2, padx=(0, 4))
        ttk.Button(buttons, text="Save PNGs", command=self.save_pngs).grid(row=0, column=3)

        opts = ttk.LabelFrame(left, text="Plot options", padding=8)
        opts.grid(row=1, column=0, sticky="ew", pady=(8, 8))
        ttk.Label(opts, text="PDD norm").grid(row=0, column=0, sticky="w")
        ttk.Combobox(
            opts,
            textvariable=self.pdd_norm_mode,
            values=["Max = 100%", "Raw"],
            width=15,
            state="readonly",
        ).grid(row=0, column=1, sticky="ew", padx=(6, 0))

        ttk.Label(opts, text="Profile norm").grid(row=1, column=0, sticky="w", pady=(4, 0))
        ttk.Combobox(
            opts,
            textvariable=self.profile_norm_mode,
            values=["CAX = 100%", "Max = 100%", "Raw"],
            width=15,
            state="readonly",
        ).grid(row=1, column=1, sticky="ew", padx=(6, 0), pady=(4, 0))

        vis_buttons = ttk.Frame(opts)
        vis_buttons.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        ttk.Button(vis_buttons, text="Show all", command=self.show_all_curves).grid(row=0, column=0, padx=(0, 4))
        ttk.Button(vis_buttons, text="Hide all", command=self.hide_all_curves).grid(row=0, column=1, padx=(0, 4))
        ttk.Button(vis_buttons, text="Invert selected", command=self.invert_selected_curves).grid(row=0, column=2, padx=(0, 4))
        ttk.Button(vis_buttons, text="10x10 compare", command=self.apply_default_compare_visibility).grid(row=0, column=3)

        ttk.Checkbutton(opts, text="Grid", variable=self.show_grid, command=self.redraw).grid(
            row=3, column=0, columnspan=2, sticky="w", pady=(6, 0)
        )
        opts.columnconfigure(1, weight=1)

        self.pdd_norm_mode.trace_add("write", lambda *_: self.redraw())
        self.profile_norm_mode.trace_add("write", lambda *_: self.redraw())

        table_frame = ttk.LabelFrame(left, text="Loaded curves", padding=4)
        table_frame.grid(row=2, column=0, sticky="nsew")
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)

        columns = ("show", "kind", "source", "points", "range", "details")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="tree headings", selectmode="extended", height=18)
        self.tree.heading("#0", text="Name")
        self.tree.heading("show", text="Plot")
        self.tree.heading("kind", text="Type")
        self.tree.heading("source", text="Src")
        self.tree.heading("points", text="N")
        self.tree.heading("range", text="Range")
        self.tree.heading("details", text="Details")
        # Keep the sidebar compact.  Wide fixed columns force the left pane
        # to occupy nearly half the screen, so only the descriptive columns stretch.
        self.tree.column("#0", width=210, minwidth=90, stretch=True)
        self.tree.column("show", width=42, minwidth=36, stretch=False, anchor="center")
        self.tree.column("kind", width=62, minwidth=50, stretch=False, anchor="center")
        self.tree.column("source", width=45, minwidth=40, stretch=False, anchor="center")
        self.tree.column("points", width=48, minwidth=42, stretch=False, anchor="e")
        self.tree.column("range", width=82, minwidth=65, stretch=False)
        self.tree.column("details", width=125, minwidth=70, stretch=True)
        self.tree.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scroll.set)
        scroll.grid(row=0, column=1, sticky="ns")
        self.tree.bind("<<TreeviewSelect>>", lambda event: self._show_metadata())
        self.tree.bind("<Button-1>", self._on_tree_click)
        self.tree.bind("<Double-1>", self._on_tree_double_click)
        self.tree.bind("<space>", lambda event: self.invert_selected_curves())

        meta_frame = ttk.LabelFrame(left, text="Selected metadata", padding=4)
        meta_frame.grid(row=3, column=0, sticky="ew", pady=(8, 0))
        self.meta_text = tk.Text(meta_frame, height=9, width=35, wrap="none")
        self.meta_text.pack(fill="both", expand=True)

        self.notebook = ttk.Notebook(right)
        self.notebook.grid(row=0, column=0, sticky="nsew")

        self.pdd_fig = self.Figure(figsize=(7, 5), dpi=100)
        self.pdd_ax = self.pdd_fig.add_subplot(111)
        self.pdd_canvas, self.pdd_toolbar = self._add_plot_tab("PDD overlay", self.pdd_fig)

        self.prof_fig = self.Figure(figsize=(7, 5), dpi=100)
        self.prof_ax = self.prof_fig.add_subplot(111)
        self.prof_canvas, self.prof_toolbar = self._add_plot_tab("Profile overlay", self.prof_fig)

        self.status = tk.StringVar(value="Load .mcc, .asc, and/or .csv files to begin.")
        ttk.Label(root, textvariable=self.status, anchor="w").grid(row=1, column=0, sticky="ew", padx=8, pady=(0, 4))

        # Apply the initial 25% split after Tk has calculated the real window size.
        root.after(100, self._set_initial_sash_position)

    def _set_initial_sash_position(self) -> None:
        """Place the left sidebar at about 25% of the available window width."""
        try:
            self.root.update_idletasks()
            width = self.main_pane.winfo_width() or self.root.winfo_width()
            if width <= 1:
                return
            left_width = max(280, int(width * 0.25))
            # Leave enough room for the plot pane even on smaller screens.
            left_width = min(left_width, max(280, width - 520))
            self.main_pane.sashpos(0, left_width)
        except Exception:
            # Non-fatal: the UI still works; only the initial sash placement failed.
            pass

    def _add_plot_tab(self, title: str, fig):
        frame = self.ttk.Frame(self.notebook)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        canvas = self.FigureCanvasTkAgg(fig, master=frame)
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar_frame = self.ttk.Frame(frame)
        toolbar_frame.grid(row=1, column=0, sticky="ew")
        toolbar = self.NavigationToolbar2Tk(canvas, toolbar_frame)
        toolbar.update()
        self.notebook.add(frame, text=title)
        return canvas, toolbar

    def choose_files(self) -> None:
        filenames = self.filedialog.askopenfilenames(
            title="Choose MCC/ASC/CSV files",
            filetypes=[("Beam data files", "*.mcc *.asc *.ASC *.csv"), ("MCC files", "*.mcc"), ("ASC files", "*.asc *.ASC"), ("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if filenames:
            self.add_files([Path(f) for f in filenames])

    def autoload_defaults(self) -> None:
        """Clear current data, load the default reference + dose-engine files, then apply compare visibility."""
        paths = default_autoload_paths()
        if not paths:
            self.messagebox.showwarning(
                "Auto-load files not found",
                "Could not find the default reference ASC or dose-engine CSV files.\n\n"
                "Expected reference:\n"
                "  FINAL_DATA/PROFILE_PDD/X06 OPEN 10X10 PDDXY SMOOTHED.ASC\n\n"
                "Expected dose-engine CSVs:\n"
                "  ../build/dose_crossline_200mm.csv\n"
                "  ../build/dose_crossline_100mm.csv\n"
                "  ../build/dose_cax.csv",
            )
            return
        self.clear()
        self.add_files(paths)
        self.apply_default_compare_visibility()

    def apply_default_compare_visibility(self) -> None:
        """Show dose-engine PDD/crosslines and matching measured ASC crossline/PDD curves.

        The preset is intended for the current 6 MV 10x10 comparison:
          - Show all CSV curves: dose_cax, dose_crossline_100mm, dose_crossline_200mm.
          - Show measured ASC PDD.
          - Show measured ASC crossplane/crossline profiles only at depths matching the CSV crosslines.
          - Hide measured inplane profiles and non-matching depths.
        """
        if not self.curves:
            return

        # Depths present in dose-engine profile CSVs, e.g. 100 and 200 mm.
        csv_profile_depths = set()
        for c in self.curves:
            if c.source == "CSV" and c.kind == "PROFILE":
                d = _depth_as_float_mm(c.metadata.get("SCAN_DEPTH", ""))
                if d is not None:
                    csv_profile_depths.add(round(d, 3))

        changed = 0
        shown = 0
        for idx, c in enumerate(self.curves):
            show = False
            if c.source == "CSV":
                # Dose-engine curves are the curves being compared, so keep them visible.
                show = True
            elif c.source == "ASC" and c.kind == "PDD":
                # Reference measured PDD for dose_cax.
                show = True
            elif c.source == "ASC" and c.kind == "PROFILE":
                curve_type = c.metadata.get("SCAN_CURVETYPE", "").lower()
                is_crossline = "cross" in curve_type or c.metadata.get("ASC_AXIS", "").upper() == "X"
                d = _depth_as_float_mm(c.metadata.get("SCAN_DEPTH", ""))
                depth_matches = d is not None and round(d, 3) in csv_profile_depths
                show = is_crossline and depth_matches

            if idx < len(self.curve_visible) and self.curve_visible[idx] != show:
                changed += 1
            if show:
                shown += 1
            self._set_curve_visible(idx, show)

        self.redraw()
        self.status.set(f"Applied 10x10 compare view: showing {shown}/{len(self.curves)} curves; changed {changed} toggles.")

    def add_files(self, files: Sequence[Path]) -> None:
        new_curves = load_curves(files)
        if not new_curves:
            self.messagebox.showwarning("No curves found", "No PDD or profile curves could be read from the selected files.")
            return
        self.curves.extend(new_curves)
        self.curve_visible.extend([True] * len(new_curves))
        self._refresh_tree()
        self.redraw()
        n_pdd = sum(c.kind == "PDD" for c in self.curves)
        n_prof = sum(c.kind == "PROFILE" for c in self.curves)
        self.status.set(f"Loaded {len(self.curves)} curves: {n_pdd} PDD, {n_prof} profile.")

    def clear(self) -> None:
        self.curves.clear()
        self.curve_visible.clear()
        self._refresh_tree()
        self.meta_text.delete("1.0", "end")
        self.redraw()
        self.status.set("Cleared all loaded curves.")

    def _refresh_tree(self) -> None:
        self.tree.delete(*self.tree.get_children())
        for idx, c in enumerate(self.curves):
            x_min = float(np.nanmin(c.x)) if len(c.x) else math.nan
            x_max = float(np.nanmax(c.x)) if len(c.x) else math.nan
            details = []
            if c.metadata.get("ENERGY"):
                details.append(f"{c.metadata['ENERGY']} MV")
            if c.metadata.get("FILTER"):
                details.append(c.metadata["FILTER"])
            if c.metadata.get("FIELD_INPLANE") and c.metadata.get("FIELD_CROSSPLANE"):
                details.append(f"{c.metadata['FIELD_INPLANE']}x{c.metadata['FIELD_CROSSPLANE']} mm")
            if c.metadata.get("SCAN_DEPTH"):
                details.append(f"d={c.metadata['SCAN_DEPTH']} mm")
            self.tree.insert(
                "",
                "end",
                iid=str(idx),
                text=c.name,
                values=("☑" if self.curve_visible[idx] else "☐", c.kind, c.source, len(c.x), f"{x_min:g} to {x_max:g}", "; ".join(details)),
            )

    def _visible_curves(self) -> List[Curve]:
        return [c for i, c in enumerate(self.curves) if i < len(self.curve_visible) and self.curve_visible[i]]

    def _set_curve_visible(self, idx: int, visible: bool) -> None:
        if idx < 0 or idx >= len(self.curve_visible):
            return
        self.curve_visible[idx] = visible
        if self.tree.exists(str(idx)):
            values = list(self.tree.item(str(idx), "values"))
            if values:
                values[0] = "☑" if visible else "☐"
                self.tree.item(str(idx), values=values)

    def _toggle_curve_visible(self, idx: int) -> None:
        if idx < 0 or idx >= len(self.curve_visible):
            return
        self._set_curve_visible(idx, not self.curve_visible[idx])

    def _on_tree_click(self, event) -> None:
        # Toggle only when the user clicks the Plot checkbox column.
        region = self.tree.identify("region", event.x, event.y)
        column = self.tree.identify_column(event.x)
        row = self.tree.identify_row(event.y)
        if region == "cell" and column == "#1" and row:
            try:
                self._toggle_curve_visible(int(row))
                self.redraw()
            except Exception:
                pass
            return "break"
        return None

    def _on_tree_double_click(self, event) -> str:
        row = self.tree.identify_row(event.y)
        if row:
            try:
                self._toggle_curve_visible(int(row))
                self.redraw()
            except Exception:
                pass
        return "break"

    def show_all_curves(self) -> None:
        for i in range(len(self.curve_visible)):
            self._set_curve_visible(i, True)
        self.redraw()

    def hide_all_curves(self) -> None:
        for i in range(len(self.curve_visible)):
            self._set_curve_visible(i, False)
        self.redraw()

    def invert_selected_curves(self) -> None:
        selection = self.tree.selection()
        if not selection:
            return
        for iid in selection:
            try:
                self._toggle_curve_visible(int(iid))
            except Exception:
                pass
        self.redraw()

    def _show_metadata(self) -> None:
        self.meta_text.delete("1.0", "end")
        selection = self.tree.selection()
        if not selection:
            return
        try:
            c = self.curves[int(selection[0])]
        except Exception:
            return
        lines = [f"Name: {c.name}", f"File: {c.path}", f"Type: {c.kind}", f"Source: {c.source}", f"Points: {len(c.x)}", ""]
        for key in sorted(c.metadata):
            lines.append(f"{key}: {c.metadata[key]}")
        self.meta_text.insert("1.0", "\n".join(lines))

    def redraw(self) -> None:
        self._plot_pdds()
        self._plot_profiles()

    def _plot_pdds(self) -> None:
        ax = self.pdd_ax
        ax.clear()
        curves = [c for c in self._visible_curves() if c.kind == "PDD"]
        for c in curves:
            y = normalize_y(c, self.pdd_norm_mode.get())
            ax.plot(c.x, y, label=c.label(), linewidth=1.8)
        ax.set_title("PDD overlay")
        ax.set_xlabel("Depth (mm)")
        ax.set_ylabel("Dose (%)" if self.pdd_norm_mode.get() != "Raw" else "Dose / signal")
        if self.show_grid.get():
            ax.grid(True, alpha=0.35)
        if curves:
            ax.legend(fontsize=8, loc="best")
        else:
            ax.text(0.5, 0.5, "No PDD curves loaded", ha="center", va="center", transform=ax.transAxes)
        self.pdd_fig.tight_layout()
        self.pdd_canvas.draw_idle()

    def _plot_profiles(self) -> None:
        ax = self.prof_ax
        ax.clear()
        curves = [c for c in self._visible_curves() if c.kind == "PROFILE"]
        for c in curves:
            y = normalize_y(c, self.profile_norm_mode.get())
            ax.plot(c.x, y, label=c.label(), linewidth=1.8)
        ax.set_title("Profile overlay")
        ax.set_xlabel("Off-axis position (mm)")
        ax.set_ylabel("Dose (%)" if self.profile_norm_mode.get() != "Raw" else "Dose / signal")
        if self.show_grid.get():
            ax.grid(True, alpha=0.35)
        if curves:
            ax.legend(fontsize=8, loc="best")
        else:
            ax.text(0.5, 0.5, "No profile curves loaded", ha="center", va="center", transform=ax.transAxes)
        self.prof_fig.tight_layout()
        self.prof_canvas.draw_idle()

    def save_pngs(self) -> None:
        directory = self.filedialog.askdirectory(title="Choose folder for PNG export")
        if not directory:
            return
        out_dir = Path(directory)
        pdd_path = out_dir / "pdd_overlay.png"
        prof_path = out_dir / "profile_overlay.png"
        self.pdd_fig.savefig(pdd_path, dpi=150)
        self.prof_fig.savefig(prof_path, dpi=150)
        self.status.set(f"Saved {pdd_path.name} and {prof_path.name}.")

    def run(self) -> None:
        self.root.mainloop()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Overlay measured MCC/ASC and dose-engine CSV PDD/profile data.")
    parser.add_argument("files", nargs="*", type=Path, help="Optional MCC/ASC/CSV files to load at startup.")
    parser.add_argument("--summary", action="store_true", help="Print parsed curve summary and exit without opening the UI.")
    parser.add_argument("--no-autoload", action="store_true", help="Do not auto-load the default 10x10 reference/dose-engine files when no files are supplied.")
    args = parser.parse_args(argv)

    startup_files = list(args.files)
    apply_compare_view = False
    if not startup_files and not args.no_autoload:
        startup_files = default_autoload_paths()
        apply_compare_view = bool(startup_files)

    if args.summary:
        curves = load_curves(startup_files)
        print(curve_summary(curves) if curves else "No curves found.")
        return 0

    app = BeamDataViewer(initial_files=startup_files, apply_compare_view=apply_compare_view)
    app.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
