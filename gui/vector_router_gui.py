#!/usr/bin/env python3
"""UK Vector Router - desktop performance lab.

Cross-platform Tkinter front-end for the uk_vector_router engine.
Features:
  * import/preprocess OSM PBF or XML data through the C engine
  * render the preprocessed road network with pan / zoom
  * route between coordinates, measure wall-clock and in-engine time
  * compare bidirectional, unidirectional and Dijkstra searches
  * show a deterministic efficiency score and the route on the map

Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import json
import math
import mmap
import os
import queue
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import threading
import time
import tkinter as tk
from array import array
from dataclasses import dataclass, field
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText
from typing import Callable, Optional
import zlib

APP_TITLE = "UK Vector Router - Performance Lab"
APP_VERSION = "1.0"

VGRAPH_MAGIC = 0x564D5232  # "VMR2"
VGRAPH_HEADER = struct.Struct("=9I4f")

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = ROOT / "data" / "europe" / "uk" / "uk_roadmap.db"
DEFAULT_VGRAPH = ROOT / "data" / "europe" / "uk" / "uk_roadmap.vgraph"

# Dark palette
C_BG = "#0d1117"
C_PANEL = "#151b23"
C_CARD = "#1b2430"
C_BORDER = "#2a3441"
C_FG = "#e6edf3"
C_MUTED = "#8b98a5"
C_ACCENT = "#ffb000"
C_ACCENT2 = "#3fb950"
C_RED = "#f85149"
C_BLUE = "#58a6ff"
C_NET = "#24405c"
C_NET_HI = "#2f5a80"

MODES = [("bidir", "Bidirectional 8D"), ("unidir", "Unidirectional 8D"), ("dijkstra", "Dijkstra")]
PROFILES = ["balanced", "fastest", "direct", "highway", "local"]
CITIES = {
    "London": (51.5074, -0.1278),
    "Edinburgh": (55.9533, -3.1883),
    "Manchester": (53.4808, -2.2426),
    "Cardiff": (51.5072, -3.1836),
    "Belfast": (54.5973, -5.9301),
    "Inverness": (57.4778, -4.2247),
    "Dover": (51.1279, 1.3134),
    "Penzance": (50.1186, -5.5436),
}


# ---------------------------------------------------------------------------
# Engine wrapper
# ---------------------------------------------------------------------------

@dataclass
class RouteResult:
    mode: str
    found: bool
    dist_km: float = 0.0
    direct_km: float = 0.0
    est_min: float = 0.0
    hops: int = 0
    expanded: int = 0
    compute_ms: float = 0.0
    total_ms: float = 0.0
    score: float = 0.0
    coords: list[tuple[float, float]] = field(default_factory=list)
    error: str = ""


def parse_route_output(text: str) -> dict:
    out: dict = {
        "found": False, "start_osm": None, "start_lat": None, "start_lon": None,
        "goal_osm": None, "goal_lat": None, "goal_lon": None,
        "dist_km": 0.0, "direct_km": 0.0, "est_min": 0.0, "est_s": 0.0,
        "hops": 0, "expanded": 0, "compute_ms": 0.0,
    }

    m = re.search(r"Start Node\s*:\s*osm=(-?\d+)\s*\(([-\d.]+),\s*([-\d.]+)\)", text)
    if m:
        out["start_osm"] = int(m.group(1))
        out["start_lat"] = float(m.group(2))
        out["start_lon"] = float(m.group(3))

    m = re.search(r"Goal Node\s*:\s*osm=(-?\d+)\s*\(([-\d.]+),\s*([-\d.]+)\)", text)
    if m:
        out["goal_osm"] = int(m.group(1))
        out["goal_lat"] = float(m.group(2))
        out["goal_lon"] = float(m.group(3))

    if re.search(r"Found\s*:\s*YES", text, re.IGNORECASE):
        out["found"] = True

    m = re.search(r"Road Dist\s*:\s*([-\d.]+)\s*km", text)
    if m:
        out["dist_km"] = float(m.group(1))

    m = re.search(r"direct TM chord:\s*([-\d.]+)\s*km", text)
    if m:
        out["direct_km"] = float(m.group(1))

    m = re.search(r"Est Time\s*:\s*([-\d.]+)\s*min\s*\(([-\d.]+)\s*s\)", text)
    if m:
        out["est_min"] = float(m.group(1))
        out["est_s"] = float(m.group(2))

    m = re.search(r"Performance\s*:\s*(\d+)\s*hops\s*\|\s*(\d+)\s*nodes expanded\s*\|\s*([-\d.]+)\s*ms", text)
    if m:
        out["hops"] = int(m.group(1))
        out["expanded"] = int(m.group(2))
        out["compute_ms"] = float(m.group(3))

    return out


def parse_geojson(path: str | os.PathLike) -> list[tuple[float, float]]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    geom = data["features"][0]["geometry"]
    if geom["type"] == "Point":
        lon, lat = geom["coordinates"]
        return [(float(lon), float(lat))]
    return [(float(lon), float(lat)) for lon, lat in geom["coordinates"]]


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def score_route(res: dict, fallback_direct_km: float) -> float:
    """Deterministic 0-100 score.

    55% ETA vs an 80 km/h straight-line baseline, 30% direct/road efficiency,
    15% in-engine compute speed.  Higher is better.
    """
    if not res.get("found"):
        return 0.0
    direct = res.get("direct_km") or fallback_direct_km
    dist = res.get("dist_km") or 0.0
    est = res.get("est_min") or 0.0
    compute = res.get("compute_ms") or 0.0
    if direct <= 0 or dist <= 0 or est <= 0:
        return 0.0
    baseline_min = direct / 80.0 * 60.0
    time_score = max(0.0, min(1.0, baseline_min / est))
    eff_score = max(0.0, min(1.0, direct / dist))
    speed_score = max(0.0, min(1.0, 300.0 / compute)) if compute > 0 else 0.0
    return round(100.0 * (0.55 * time_score + 0.30 * eff_score + 0.15 * speed_score), 1)


def find_executable(explicit: Optional[str]) -> Optional[Path]:
    names = ["uk_vector_router.exe", "uk_vector_router"]
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))
    candidates.extend(ROOT / n for n in names)
    candidates.append(ROOT / "build" / names[0])
    candidates.append(ROOT / "build" / names[1])
    for c in candidates:
        if c.is_file():
            return c
    for n in names:
        found = shutil.which(n)
        if found:
            return Path(found)
    return None


class Engine:
    def __init__(self, exe: Optional[Path], db: Path, vgraph: Path):
        self.exe = exe
        self.db = db
        self.vgraph = vgraph

    def available(self) -> bool:
        return self.exe is not None and self.exe.is_file()

    def base_args(self) -> list[str]:
        if self.vgraph.is_file():
            return ["--load-bin", str(self.vgraph)]
        return ["--db", str(self.db)]

    def route(self, mode: str, profile: str, lat1: float, lon1: float,
              lat2: float, lon2: float, geojson: Optional[str] = None) -> tuple[RouteResult, str, float]:
        if not self.available():
            raise RuntimeError("uk_vector_router executable not found")
        args = self.base_args() + [
            "--route-coord", f"{lat1:.6f}", f"{lon1:.6f}", f"{lat2:.6f}", f"{lon2:.6f}",
            "--mode", mode, "--profile", profile,
        ]
        if geojson:
            args += ["--geojson", geojson]
        t0 = time.perf_counter()
        proc = subprocess.run([str(self.exe)] + args, capture_output=True, text=True,
                              errors="replace", timeout=900)
        elapsed = (time.perf_counter() - t0) * 1000.0
        text = (proc.stdout or "") + ("\n" + proc.stderr if proc.stderr else "")
        parsed = parse_route_output(text)
        parsed["score"] = score_route(parsed, haversine_km(lat1, lon1, lat2, lon2))
        parsed["total_ms"] = elapsed
        res = RouteResult(
            mode=mode, found=parsed["found"], dist_km=parsed["dist_km"],
            direct_km=parsed["direct_km"], est_min=parsed["est_min"],
            hops=parsed["hops"], expanded=parsed["expanded"],
            compute_ms=parsed["compute_ms"], total_ms=elapsed,
            score=parsed["score"],
        )
        if not parsed["found"]:
            res.error = "no route found"
        return res, text, elapsed

    def stream(self, args: list[str], on_line: Callable[[str], None]) -> int:
        if not self.available():
            raise RuntimeError("uk_vector_router executable not found")
        proc = subprocess.Popen([str(self.exe)] + args, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, errors="replace", bufsize=1)
        assert proc.stdout is not None
        for line in proc.stdout:
            on_line(line.rstrip())
        proc.wait()
        return proc.returncode


# ---------------------------------------------------------------------------
# Map data
# ---------------------------------------------------------------------------

@dataclass
class MapData:
    lon: array
    lat: array
    efrom: array
    eto: array
    node_count: int
    edge_count: int
    min_lon: float
    max_lon: float
    min_lat: float
    max_lat: float
    raw_node_count: int = 0
    raw_edge_count: int = 0
    contracted: int = 0
    peeled: int = 0

    @classmethod
    def load(cls, path: Path) -> "MapData":
        with open(path, "rb") as fh:
            mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
            try:
                hdr = VGRAPH_HEADER.unpack_from(mm, 0)
                magic, version, n, m = hdr[0], hdr[1], hdr[2], hdr[3]
                if magic != VGRAPH_MAGIC:
                    raise ValueError(f"bad .vgraph magic 0x{magic:08x}")
                raw_n, raw_m, contracted, peeled = hdr[4], hdr[5], hdr[6], hdr[7]
                off = VGRAPH_HEADER.size
                off += n * 8  # node_osm_id
                lat = array("d")
                lat.frombytes(mm[off:off + n * 8])
                off += n * 8
                lon = array("d")
                lon.frombytes(mm[off:off + n * 8])
                off += n * 8
                # node_x, node_y, node_comp_id, node_flags, fwd_offsets(n + 1)
                off += n * 4 + n * 4 + n * 4 + n + (n + 1) * 4
                efrom = array("I")
                efrom.frombytes(mm[off:off + m * 4])
                off += m * 4
                eto = array("I")
                eto.frombytes(mm[off:off + m * 4])
            finally:
                mm.close()
        if sys.byteorder != "little":
            for a in (lat, lon, efrom, eto):
                a.byteswap()
        if n == 0:
            return cls(lon, lat, efrom, eto, 0, 0, 0, 0, 0, 0)
        min_lon = min(lon)
        max_lon = max(lon)
        min_lat = min(lat)
        max_lat = max(lat)
        return cls(lon, lat, efrom, eto, n, m, min_lon, max_lon, min_lat, max_lat,
                   raw_n, raw_m, contracted, peeled)


def synthetic_map() -> MapData:
    """Small synthetic network used by the GUI smoke test."""
    lon: array = array("d")
    lat: array = array("d")
    efrom: array = array("I")
    eto: array = array("I")
    nodes: list[tuple[float, float]] = []
    for i in range(40):
        for j in range(30):
            nodes.append((-6.0 + i * 0.2, 50.0 + j * 0.3))
    for i in range(40):
        for j in range(30):
            idx = i * 30 + j
            lon.append(nodes[idx][0])
            lat.append(nodes[idx][1])
            if i + 1 < 40:
                efrom.append(idx)
                eto.append(idx + 30)
            if j + 1 < 30:
                efrom.append(idx)
                eto.append(idx + 1)
    return MapData(lon, lat, efrom, eto, len(lon), len(efrom),
                   min(lon), max(lon), min(lat), max(lat))


# ---------------------------------------------------------------------------
# Map canvas
# ---------------------------------------------------------------------------

class MapView:
    def __init__(self, canvas: tk.Canvas, app: "VectorRouterApp"):
        self.canvas = canvas
        self.app = app
        self.data: Optional[MapData] = None
        self.center_lon = -2.5
        self.center_lat = 54.5
        self.k = 12.0  # pixels per degree latitude
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.step = 1
        self.route: list[tuple[float, float]] = []
        self.start: Optional[tuple[float, float]] = None
        self.goal: Optional[tuple[float, float]] = None
        self.dragging = False
        self.drag_x = 0
        self.drag_y = 0
        self.last_redraw = 0.0
        self.after_id: Optional[str] = None
        self.pick_mode = tk.StringVar(value="off")

        canvas.configure(background=C_BG, highlightthickness=0, bd=0)
        canvas.bind("<ButtonPress-1>", self.on_press)
        canvas.bind("<B1-Motion>", self.on_drag)
        canvas.bind("<ButtonRelease-1>", self.on_release)
        canvas.bind("<Motion>", self.on_motion)
        canvas.bind("<Configure>", lambda _e: self.redraw())
        canvas.bind("<MouseWheel>", lambda e: self.on_wheel(e))
        canvas.bind("<Button-4>", lambda e: self.on_wheel(e, 1))
        canvas.bind("<Button-5>", lambda e: self.on_wheel(e, -1))
        canvas.bind("<Button-2>", self.on_middle)

    # -- geometry ---------------------------------------------------------
    def size(self) -> tuple[int, int]:
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()
        return max(w, 1), max(h, 1)

    def aspect(self) -> float:
        return max(0.2, math.cos(math.radians(self.center_lat)))

    def project(self, lon: float, lat: float) -> tuple[float, float]:
        w, h = self.size()
        x = (lon - self.center_lon) * self.k * self.aspect() + w * 0.5 + self.pan_x
        y = -(lat - self.center_lat) * self.k + h * 0.5 + self.pan_y
        return x, y

    def unproject(self, x: float, y: float) -> tuple[float, float]:
        w, h = self.size()
        lon = (x - w * 0.5 - self.pan_x) / (self.k * self.aspect()) + self.center_lon
        lat = -(y - h * 0.5 - self.pan_y) / self.k + self.center_lat
        return lon, lat

    def fit(self, min_lon: float, max_lon: float, min_lat: float, max_lat: float, pad: float = 0.08):
        w, h = self.size()
        if w <= 1 or h <= 1:
            return
        span_lon = max(max_lon - min_lon, 1e-3)
        span_lat = max(max_lat - min_lat, 1e-3)
        self.center_lon = (min_lon + max_lon) * 0.5
        self.center_lat = (min_lat + max_lat) * 0.5
        self.k = min(w * (1 - pad) / (span_lon * self.aspect()),
                     h * (1 - pad) / span_lat)
        self.pan_x = 0.0
        self.pan_y = 0.0

    def fit_data(self):
        if not self.data:
            return
        self.fit(self.data.min_lon, self.data.max_lon, self.data.min_lat, self.data.max_lat)
        self.redraw()

    def fit_route(self):
        if not self.route:
            return
        lons = [p[0] for p in self.route]
        lats = [p[1] for p in self.route]
        self.fit(min(lons), max(lons), min(lats), max(lats))
        self.redraw()

    # -- rendering --------------------------------------------------------
    def redraw(self):
        if self.after_id is not None:
            try:
                self.canvas.after_cancel(self.after_id)
            except Exception:
                pass
            self.after_id = None
        c = self.canvas
        c.delete("all")
        w, h = self.size()
        if w <= 1 or h <= 1:
            return
        self.draw_background(w, h)
        if self.data:
            self.render_raster(w, h)
        self.draw_route()
        self.draw_markers()
        self.draw_scalebar(w, h)
        self.draw_hud(w, h)

    def draw_background(self, w: int, h: int):
        # Kept for API compatibility; the raster renderer draws the backdrop.
        self.canvas.create_rectangle(0, 0, w, h, fill=C_BG, outline="", width=0, tags="bg")

    @staticmethod
    def line(buf: bytearray, w: int, h: int, x0: float, y0: float,
             x1: float, y1: float, rgb: tuple[int, int, int]):
        xi0, yi0, xi1, yi1 = int(x0), int(y0), int(x1), int(y1)
        dx = abs(xi1 - xi0)
        sx = 1 if xi0 < xi1 else -1
        dy = -abs(yi1 - yi0)
        sy = 1 if yi0 < yi1 else -1
        err = dx + dy
        r, g, b = rgb
        while True:
            if 0 <= xi0 < w and 0 <= yi0 < h:
                i = (yi0 * w + xi0) * 3
                buf[i] = r
                buf[i + 1] = g
                buf[i + 2] = b
            if xi0 == xi1 and yi0 == yi1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                xi0 += sx
            if e2 <= dx:
                err += dx
                yi0 += sy

    @staticmethod
    def png_bytes(buf: bytearray, w: int, h: int) -> bytes:
        raw = bytearray()
        stride = w * 3
        for y in range(h):
            raw.append(0)
            raw += buf[y * stride:(y + 1) * stride]

        def chunk(tag: bytes, data: bytes) -> bytes:
            return (struct.pack(">I", len(data)) + tag + data +
                    struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))

        return (b"\x89PNG\r\n\x1a\n" +
                chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)) +
                chunk(b"IDAT", zlib.compress(bytes(raw), 6)) +
                chunk(b"IEND", b""))

    def render_raster(self, w: int, h: int):
        d = self.data
        assert d is not None
        buf = bytearray(bytes((13, 17, 23)) * (w * h))

        # graticule
        step = 0.05 if self.k > 2200 else 0.1 if self.k > 900 else 0.5 if self.k > 180 else 1.0 if self.k > 60 else 2.0
        lon0 = math.floor(self.center_lon / step) * step
        lat0 = math.floor(self.center_lat / step) * step
        grid = (20, 27, 35)
        for i in range(-40, 41):
            x, _ = self.project(lon0 + i * step, self.center_lat)
            if 0 <= x < w:
                self.line(buf, w, h, x, 0, x, h, grid)
            _, y = self.project(self.center_lon, lat0 + i * step)
            if 0 <= y < h:
                self.line(buf, w, h, 0, y, w, y, grid)

        # road network
        self.step = max(1, d.edge_count // 40000)
        ash = self.aspect()
        k = self.k
        clon, clat = self.center_lon, self.center_lat
        ox, oy = w * 0.5 + self.pan_x, h * 0.5 + self.pan_y
        lon_arr, lat_arr = d.lon, d.lat
        efrom, eto = d.efrom, d.eto
        net = (36, 64, 92)
        for i in range(0, d.edge_count, self.step):
            a = efrom[i]
            b = eto[i]
            x1 = (lon_arr[a] - clon) * k * ash + ox
            y1 = -(lat_arr[a] - clat) * k + oy
            x2 = (lon_arr[b] - clon) * k * ash + ox
            y2 = -(lat_arr[b] - clat) * k + oy
            if (x1 < -2 and x2 < -2) or (x1 > w + 2 and x2 > w + 2):
                continue
            if (y1 < -2 and y2 < -2) or (y1 > h + 2 and y2 > h + 2):
                continue
            self.line(buf, w, h, x1, y1, x2, y2, net)

        self.photo = tk.PhotoImage(data=self.png_bytes(buf, w, h))
        self.canvas.create_image(0, 0, anchor="nw", image=self.photo, tags="raster")

    def draw_route(self):
        if len(self.route) < 2:
            return
        pts: list[float] = []
        for lon, lat in self.route:
            x, y = self.project(lon, lat)
            pts.extend((x, y))
        self.canvas.create_line(*pts, fill=C_ACCENT, width=4, capstyle="round",
                                joinstyle="round", tags="route")
        self.canvas.create_line(*pts, fill="#ffe0a0", width=1, tags="route")

    def draw_markers(self):
        for point, color, label in ((self.start, C_ACCENT2, "A"), (self.goal, C_RED, "B")):
            if not point:
                continue
            x, y = self.project(point[0], point[1])
            r = 6
            self.canvas.create_oval(x - r, y - r, x + r, y + r, fill=color,
                                    outline="#0d1117", width=2, tags="marker")
            self.canvas.create_text(x, y - 14, text=label, fill=color,
                                    font=("Segoe UI", 9, "bold"), tags="marker")

    def draw_scalebar(self, w: int, h: int):
        km_per_px = 111.32 / max(self.k, 1e-6)
        target_px = 140
        raw = km_per_px * target_px
        nice = 1
        for cand in (5, 10, 20, 25, 50, 100, 200, 250, 500):
            if cand <= raw:
                nice = cand
        px = nice / km_per_px
        x0, y0 = 18, h - 52
        self.canvas.create_line(x0, y0, x0 + px, y0, fill=C_MUTED, width=2, tags="ui")
        self.canvas.create_line(x0, y0 - 4, x0, y0 + 4, fill=C_MUTED, width=1, tags="ui")
        self.canvas.create_line(x0 + px, y0 - 4, x0 + px, y0 + 4, fill=C_MUTED, width=1, tags="ui")
        self.canvas.create_text(x0 + px / 2, y0 - 12, text=f"{nice} km",
                                fill=C_MUTED, font=("Segoe UI", 8), tags="ui")

    def draw_hud(self, w: int, h: int):
        c = self.canvas
        if self.data:
            c.create_text(w - 14, 16, anchor="ne",
                          text=f"{self.data.node_count:,} nodes  |  {self.data.edge_count:,} edges"
                               f"  |  drawing {min(self.data.edge_count, self.data.edge_count // max(self.step,1) + 1):,}",
                          fill="#5c6b7a", font=("Segoe UI", 8), tags="ui")
        if self.route:
            c.create_text(w - 14, 34, anchor="ne",
                          text=f"route points: {len(self.route):,}",
                          fill=C_ACCENT, font=("Segoe UI", 8), tags="ui")

    # -- interaction ------------------------------------------------------
    def on_press(self, event):
        self.dragging = True
        self.drag_x = event.x
        self.drag_y = event.y

    def on_drag(self, event):
        if not self.dragging:
            return
        dx = event.x - self.drag_x
        dy = event.y - self.drag_y
        if abs(dx) < 1 and abs(dy) < 1:
            return
        self.drag_x = event.x
        self.drag_y = event.y
        self.pan_x += dx
        self.pan_y += dy
        for tag in ("raster", "route", "marker", "ui"):
            self.canvas.move(tag, dx, dy)
        self.app.on_view_changed()

    def on_release(self, event):
        if not self.dragging:
            return
        self.dragging = False
        if self.pick_mode.get() != "off" and self.data:
            lon, lat = self.unproject(event.x, event.y)
            self.app.pick_coordinate(lon, lat, self.pick_mode.get())
        self.redraw()

    def on_middle(self, event):
        self.fit_data()

    def on_wheel(self, event, direction=None):
        if not self.data:
            return
        if direction is None:
            delta = getattr(event, "delta", 0)
            direction = 1 if delta > 0 else (-1 if delta < 0 else 0)
        if getattr(event, "num", 0) == 4:
            direction = 1
        elif getattr(event, "num", 0) == 5:
            direction = -1
        if direction == 0:
            return
        factor = 1.2 if direction > 0 else 1 / 1.2
        wx, wy = self.unproject(event.x, event.y)
        self.k *= factor
        # keep the world point under the cursor fixed
        w, h = self.size()
        self.pan_x = event.x - (wx - self.center_lon) * self.k * self.aspect() - w * 0.5
        self.pan_y = event.y + (wy - self.center_lat) * self.k - h * 0.5
        self.redraw()

    def on_motion(self, event):
        if not self.data:
            return
        lon, lat = self.unproject(event.x, event.y)
        self.app.on_coord_hover(lat, lon)
        if self.dragging:
            self.app.on_coord_hover(lat, lon)

    def set_route(self, coords: list[tuple[float, float]]):
        self.route = coords
        self.redraw()


# ---------------------------------------------------------------------------
# Main application
# ---------------------------------------------------------------------------

class VectorRouterApp:
    def __init__(self, root: tk.Tk, engine: Engine, map_data: Optional[MapData] = None):
        self.root = root
        self.engine = engine
        self.busy = False
        self.ui_queue: queue.Queue = queue.Queue()
        self.busy_label = tk.StringVar(value="Ready")
        self.status_var = tk.StringVar(value="Ready")
        self.map_status = tk.StringVar(value="No map loaded")
        self.coord_var = tk.StringVar(value="--")
        self.score_var = tk.StringVar(value="--")

        root.title(f"{APP_TITLE}  v{APP_VERSION}")
        root.configure(background=C_BG)
        root.geometry("1440x900")
        root.minsize(1100, 700)

        self.style = ttk.Style(root)
        try:
            self.style.theme_use("clam")
        except tk.TclError:
            pass
        self.setup_style()

        self.container = ttk.Frame(root, style="App.TFrame")
        self.container.pack(fill="both", expand=True)
        self.container.grid_columnconfigure(1, weight=1)
        self.container.grid_rowconfigure(1, weight=1)

        self.build_header()
        self.build_sidebar()
        self.build_center()

        self.map = MapView(self.canvas, self)
        self.map.data = map_data
        if map_data is not None:
            self.root.after(50, self.map.fit_data)
            self.map_status.set(f"Loaded {map_data.node_count:,} nodes / {map_data.edge_count:,} edges")

        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.bind_shortcuts()
        self.root.after(80, self.poll_ui)
        self.log("info", f"Engine: {self.engine.exe if self.engine.exe else 'not found'}")
        self.log("info", f"DB: {self.engine.db}")
        self.log("info", f"Cache: {self.engine.vgraph} {'(present)' if self.engine.vgraph.is_file() else '(missing)'}")
        if not self.engine.available():
            self.log("error", "uk_vector_router executable not found; build it first.")

        if map_data is None and self.engine.vgraph.is_file():
            self.root.after(120, self.load_map_async)

    # -- styling ----------------------------------------------------------
    def setup_style(self):
        s = self.style
        s.configure("App.TFrame", background=C_BG)
        s.configure("Panel.TFrame", background=C_PANEL)
        s.configure("Card.TFrame", background=C_CARD)
        s.configure("Header.TFrame", background=C_PANEL)
        s.configure("TLabel", background=C_PANEL, foreground=C_FG, font=("Segoe UI", 9))
        s.configure("Card.TLabel", background=C_CARD, foreground=C_FG, font=("Segoe UI", 9))
        s.configure("Muted.TLabel", background=C_PANEL, foreground=C_MUTED, font=("Segoe UI", 8))
        s.configure("CardMuted.TLabel", background=C_CARD, foreground=C_MUTED, font=("Segoe UI", 8))
        s.configure("Title.TLabel", background=C_PANEL, foreground=C_FG, font=("Segoe UI", 15, "bold"))
        s.configure("Sub.TLabel", background=C_PANEL, foreground=C_MUTED, font=("Segoe UI", 8))
        s.configure("H2.TLabel", background=C_CARD, foreground=C_ACCENT, font=("Segoe UI", 9, "bold"))
        s.configure("TButton", background=C_CARD, foreground=C_FG, bordercolor=C_BORDER,
                    focuscolor=C_CARD, padding=(10, 6), font=("Segoe UI", 9), borderwidth=1)
        s.map("TButton", background=[("active", "#243244")], foreground=[("disabled", "#5c6b7a")])
        s.configure("Accent.TButton", background="#3a2a00", foreground=C_ACCENT,
                    bordercolor="#5a4300", padding=(10, 7), font=("Segoe UI", 9, "bold"))
        s.map("Accent.TButton", background=[("active", "#4a3600")])
        s.configure("Go.TButton", background="#0f2a17", foreground=C_ACCENT2,
                    bordercolor="#1c4a2a", padding=(10, 7), font=("Segoe UI", 9, "bold"))
        s.map("Go.TButton", background=[("active", "#143a20")])
        s.configure("TEntry", fieldbackground="#0f141a", foreground=C_FG,
                    insertcolor=C_FG, bordercolor=C_BORDER, lightcolor=C_BORDER,
                    darkcolor=C_BORDER, padding=4)
        s.configure("TCombobox", fieldbackground="#0f141a", background=C_CARD,
                    foreground=C_FG, arrowcolor=C_FG, bordercolor=C_BORDER,
                    lightcolor=C_BORDER, darkcolor=C_BORDER, padding=4)
        s.map("TCombobox", fieldbackground=[("readonly", "#0f141a")])
        s.configure("Treeview", background=C_CARD, fieldbackground=C_CARD,
                    foreground=C_FG, bordercolor=C_BORDER, rowheight=24,
                    font=("Segoe UI", 9))
        s.configure("Card.TRadiobutton", background=C_CARD, foreground=C_FG,
                    activebackground=C_CARD, activeforeground=C_ACCENT,
                    font=("Segoe UI", 9))
        s.map("Card.TRadiobutton", background=[("active", C_CARD)])
        s.configure("Treeview.Heading", background="#111820", foreground=C_MUTED,
                    font=("Segoe UI", 8, "bold"), borderwidth=0)
        s.map("Treeview", background=[("selected", "#243244")])
        s.configure("App.Horizontal.TProgressbar", background=C_ACCENT,
                    troughcolor="#0f141a", bordercolor=C_BORDER, lightcolor=C_ACCENT,
                    darkcolor=C_ACCENT)

    # -- widgets ----------------------------------------------------------
    def build_header(self):
        header = ttk.Frame(self.container, style="Header.TFrame", padding=(16, 12, 16, 10))
        header.grid(row=0, column=0, columnspan=2, sticky="ew")
        header.grid_columnconfigure(1, weight=1)

        left = ttk.Frame(header, style="Header.TFrame")
        left.grid(row=0, column=0, sticky="w")
        ttk.Label(left, text="UK VECTOR ROUTER", style="Title.TLabel").pack(anchor="w")
        ttk.Label(left, text="preprocessed 8D vector-corridor pathfinding  |  performance lab",
                  style="Sub.TLabel").pack(anchor="w")

        right = ttk.Frame(header, style="Header.TFrame")
        right.grid(row=0, column=2, sticky="e")
        self.status_dot = tk.Canvas(right, width=10, height=10, bg=C_PANEL, highlightthickness=0)
        self.status_dot.pack(side="left", padx=(0, 7), pady=(0, 2))
        self.dot_id = self.status_dot.create_oval(1, 1, 9, 9, fill=C_ACCENT2, outline="")
        ttk.Label(right, textvariable=self.busy_label, style="TLabel").pack(side="left")
        self.progress = ttk.Progressbar(right, mode="indeterminate", length=120, style="App.Horizontal.TProgressbar")
        self.progress.pack(side="left", padx=(14, 0))

    def build_sidebar(self):
        side = ttk.Frame(self.container, style="Panel.TFrame", width=350, padding=14)
        side.grid(row=1, column=0, sticky="ns")
        side.grid_propagate(False)
        for i in range(12):
            side.grid_rowconfigure(i, weight=0)

        # Data
        card = ttk.Frame(side, style="Card.TFrame", padding=12)
        card.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        card.grid_columnconfigure(0, weight=1)
        ttk.Label(card, text="DATA", style="H2.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(card, textvariable=self.map_status, style="CardMuted.TLabel").grid(row=1, column=0, sticky="w", pady=(2, 8))
        self.db_var = tk.StringVar(value=str(self.engine.db))
        self.vgraph_var = tk.StringVar(value=str(self.engine.vgraph))
        ttk.Label(card, text="SQLite DB", style="Card.TLabel").grid(row=2, column=0, sticky="w")
        dbrow = ttk.Frame(card, style="Card.TFrame")
        dbrow.grid(row=3, column=0, sticky="ew", pady=(2, 6))
        dbrow.grid_columnconfigure(0, weight=1)
        ttk.Entry(dbrow, textvariable=self.db_var, style="TEntry").grid(row=0, column=0, sticky="ew")
        ttk.Button(dbrow, text="...", width=3, command=self.choose_db).grid(row=0, column=1, padx=(6, 0))
        ttk.Label(card, text=".vgraph cache", style="Card.TLabel").grid(row=4, column=0, sticky="w")
        vrow = ttk.Frame(card, style="Card.TFrame")
        vrow.grid(row=5, column=0, sticky="ew", pady=(2, 8))
        vrow.grid_columnconfigure(0, weight=1)
        ttk.Entry(vrow, textvariable=self.vgraph_var, style="TEntry").grid(row=0, column=0, sticky="ew")
        ttk.Button(vrow, text="...", width=3, command=self.choose_vgraph).grid(row=0, column=1, padx=(6, 0))
        btns = ttk.Frame(card, style="Card.TFrame")
        btns.grid(row=6, column=0, sticky="ew")
        btns.grid_columnconfigure(0, weight=1)
        btns.grid_columnconfigure(1, weight=1)
        ttk.Button(btns, text="Import PBF / OSM", style="Accent.TButton",
                   command=self.choose_import).grid(row=0, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(btns, text="Prepare .vgraph", command=self.prepare_async).grid(row=0, column=1, sticky="ew", padx=(4, 0))
        ttk.Button(card, text="Reload map", command=self.load_map_async).grid(row=7, column=0, sticky="ew", pady=(6, 0))

        # Route
        card = ttk.Frame(side, style="Card.TFrame", padding=12)
        card.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        card.grid_columnconfigure(0, weight=1)
        card.grid_columnconfigure(1, weight=1)
        ttk.Label(card, text="ROUTE", style="H2.TLabel").grid(row=0, column=0, columnspan=2, sticky="w")
        self.start_lat = tk.StringVar(value="51.5074")
        self.start_lon = tk.StringVar(value="-0.1278")
        self.goal_lat = tk.StringVar(value="55.9533")
        self.goal_lon = tk.StringVar(value="-3.1883")
        self.profile_var = tk.StringVar(value="balanced")
        self.mode_var = tk.StringVar(value="bidir")
        self.pick_var = tk.StringVar(value="off")

        ttk.Label(card, text="Start  (lat, lon)", style="Card.TLabel").grid(row=1, column=0, sticky="w", pady=(4, 0))
        srow = ttk.Frame(card, style="Card.TFrame")
        srow.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(2, 4))
        srow.grid_columnconfigure(0, weight=1)
        srow.grid_columnconfigure(1, weight=1)
        ttk.Entry(srow, textvariable=self.start_lat, style="TEntry").grid(row=0, column=0, sticky="ew")
        ttk.Entry(srow, textvariable=self.start_lon, style="TEntry").grid(row=0, column=1, sticky="ew", padx=(5, 0))
        ttk.Label(card, text="Goal  (lat, lon)", style="Card.TLabel").grid(row=3, column=0, sticky="w")
        grow = ttk.Frame(card, style="Card.TFrame")
        grow.grid(row=4, column=0, columnspan=2, sticky="ew", pady=(2, 6))
        grow.grid_columnconfigure(0, weight=1)
        grow.grid_columnconfigure(1, weight=1)
        ttk.Entry(grow, textvariable=self.goal_lat, style="TEntry").grid(row=0, column=0, sticky="ew")
        ttk.Entry(grow, textvariable=self.goal_lon, style="TEntry").grid(row=0, column=1, sticky="ew", padx=(5, 0))

        ttk.Label(card, text="Profile", style="Card.TLabel").grid(row=5, column=0, sticky="w")
        ttk.Combobox(card, textvariable=self.profile_var, values=PROFILES, state="readonly",
                     style="TCombobox").grid(row=6, column=0, columnspan=2, sticky="ew", pady=(2, 6))
        ttk.Label(card, text="Pick on map", style="Card.TLabel").grid(row=7, column=0, sticky="w")
        pickrow = ttk.Frame(card, style="Card.TFrame")
        pickrow.grid(row=8, column=0, columnspan=2, sticky="w", pady=(2, 8))
        for val, text in (("off", "Off"), ("start", "Start"), ("goal", "Goal")):
            ttk.Radiobutton(pickrow, text=text, value=val, variable=self.pick_var,
                            style="Card.TRadiobutton", command=self.sync_pick).pack(side="left", padx=(0, 10))

        self.run_btn = ttk.Button(card, text="Run performance test", style="Go.TButton",
                                  command=self.run_benchmark)
        self.run_btn.grid(row=9, column=0, columnspan=2, sticky="ew")
        ttk.Button(card, text="Fit route", command=lambda: self.map.fit_route()).grid(row=10, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(card, text="Export GeoJSON", command=self.export_geojson).grid(row=10, column=1, sticky="ew", pady=(6, 0), padx=(4, 0))

        # Presets
        card = ttk.Frame(side, style="Card.TFrame", padding=12)
        card.grid(row=2, column=0, sticky="ew", pady=(0, 10))
        card.grid_columnconfigure(0, weight=1)
        card.grid_columnconfigure(1, weight=1)
        ttk.Label(card, text="QUICK PRESETS", style="H2.TLabel").grid(row=0, column=0, columnspan=2, sticky="w")
        self.city_var = tk.StringVar(value="London")
        ttk.Combobox(card, textvariable=self.city_var, values=list(CITIES.keys()), state="readonly",
                     style="TCombobox").grid(row=1, column=0, columnspan=2, sticky="ew", pady=(4, 6))
        ttk.Button(card, text="Set start", command=lambda: self.set_city("start")).grid(row=2, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(card, text="Set goal", command=lambda: self.set_city("goal")).grid(row=2, column=1, sticky="ew", padx=(4, 0))

        # Display
        card = ttk.Frame(side, style="Card.TFrame", padding=12)
        card.grid(row=3, column=0, sticky="ew")
        card.grid_columnconfigure(0, weight=1)
        ttk.Label(card, text="DISPLAY", style="H2.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(card, textvariable=self.coord_var, style="CardMuted.TLabel").grid(row=1, column=0, sticky="w", pady=(3, 0))
        ttk.Label(card, text="scroll = zoom  |  drag = pan  |  middle-click = fit",
                  style="CardMuted.TLabel").grid(row=2, column=0, sticky="w")
        ttk.Label(card, text="Score = 55% ETA vs 80 km/h baseline + 30% direct/road efficiency + 15% compute speed",
                  style="CardMuted.TLabel", wraplength=300, justify="left").grid(row=3, column=0, sticky="w", pady=(4, 0))

        # Side log
        card = ttk.Frame(side, style="Card.TFrame", padding=10)
        card.grid(row=4, column=0, sticky="nsew", pady=(10, 0))
        card.grid_columnconfigure(0, weight=1)
        card.grid_rowconfigure(1, weight=1)
        side.grid_rowconfigure(4, weight=1)
        ttk.Label(card, text="LOG", style="H2.TLabel").grid(row=0, column=0, sticky="w")
        self.log_text = ScrolledText(card, height=10, bg="#0f141a", fg=C_FG,
                                     insertbackground=C_FG, relief="flat", bd=0,
                                     font=("Consolas", 8), wrap="word")
        self.log_text.grid(row=1, column=0, sticky="nsew", pady=(4, 0))
        self.log_text.tag_configure("info", foreground=C_MUTED)
        self.log_text.tag_configure("ok", foreground=C_ACCENT2)
        self.log_text.tag_configure("err", foreground=C_RED)
        self.log_text.tag_configure("hi", foreground=C_ACCENT)

    def build_center(self):
        center = ttk.Frame(self.container, style="App.TFrame", padding=(10, 0, 12, 12))
        center.grid(row=1, column=1, sticky="nsew")
        center.grid_columnconfigure(0, weight=1)
        center.grid_rowconfigure(0, weight=1)
        center.grid_rowconfigure(1, weight=0)

        map_frame = ttk.Frame(center, style="Card.TFrame", padding=2)
        map_frame.grid(row=0, column=0, sticky="nsew")
        map_frame.grid_columnconfigure(0, weight=1)
        map_frame.grid_rowconfigure(0, weight=1)
        self.canvas = tk.Canvas(map_frame, bg=C_BG, highlightthickness=0, bd=0)
        self.canvas.grid(row=0, column=0, sticky="nsew")

        # results card
        results = ttk.Frame(center, style="Card.TFrame", padding=10)
        results.grid(row=1, column=0, sticky="ew", pady=(10, 0))
        results.grid_columnconfigure(0, weight=1)
        head = ttk.Frame(results, style="Card.TFrame")
        head.grid(row=0, column=0, sticky="ew")
        head.grid_columnconfigure(1, weight=1)
        ttk.Label(head, text="PERFORMANCE COMPARISON", style="H2.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(head, textvariable=self.score_var, style="CardMuted.TLabel").grid(row=0, column=1, sticky="e")

        columns = ("mode", "found", "dist", "eta", "hops", "exp", "cpu", "total", "score")
        self.tree = ttk.Treeview(results, columns=columns, show="headings", height=3, style="Treeview")
        headings = {
            "mode": ("Search", 120), "found": ("Found", 50), "dist": ("Dist km", 80),
            "eta": ("ETA min", 75), "hops": ("Hops", 60), "exp": ("Expanded", 80),
            "cpu": ("Compute ms", 90), "total": ("Total ms", 85), "score": ("Score", 60),
        }
        for col, (text, width) in headings.items():
            self.tree.heading(col, text=text)
            self.tree.column(col, width=width, anchor="center" if col != "mode" else "w")
        self.tree.tag_configure("best", foreground=C_ACCENT)
        self.tree.grid(row=1, column=0, sticky="ew", pady=(6, 0))

    # -- logging / status -------------------------------------------------
    def post(self, fn, *args):
        self.ui_queue.put((fn, args))

    def poll_ui(self):
        try:
            while True:
                fn, args = self.ui_queue.get_nowait()
                try:
                    fn(*args)
                except Exception as exc:
                    self.log("err", f"UI callback failed: {exc}")
        except queue.Empty:
            pass
        self.root.after(50, self.poll_ui)

    def log(self, kind: str, text: str):
        ts = time.strftime("%H:%M:%S")
        self.log_text.insert("end", f"{ts}  ", "info")
        self.log_text.insert("end", text + "\n", kind)
        self.log_text.see("end")

    def set_busy(self, busy: bool, label: str = "Ready"):
        self.busy = busy
        self.busy_label.set(label)
        self.status_dot.itemconfig(self.dot_id, fill=C_ACCENT if busy else C_ACCENT2)
        if busy:
            self.progress.start(12)
            self.run_btn.configure(state="disabled")
        else:
            self.progress.stop()
            self.run_btn.configure(state="normal")

    def on_view_changed(self):
        self.coord_var.set("view updated")

    def on_coord_hover(self, lat: float, lon: float):
        self.coord_var.set(f"{lat:.5f}, {lon:.5f}")

    def sync_pick(self):
        self.map.pick_mode.set(self.pick_var.get())

    def pick_coordinate(self, lon: float, lat: float, target: str):
        if target == "start":
            self.start_lat.set(f"{lat:.5f}")
            self.start_lon.set(f"{lon:.5f}")
            self.map.start = (lon, lat)
        elif target == "goal":
            self.goal_lat.set(f"{lat:.5f}")
            self.goal_lon.set(f"{lon:.5f}")
            self.map.goal = (lon, lat)
        self.map.redraw()
        self.log("info", f"picked {target}: {lat:.5f}, {lon:.5f}")

    # -- paths ------------------------------------------------------------
    def choose_db(self):
        p = filedialog.askopenfilename(title="Select SQLite DB", filetypes=[("SQLite", "*.db"), ("All", "*.*")])
        if p:
            self.db_var.set(p)
            self.engine.db = Path(p)

    def choose_vgraph(self):
        p = filedialog.askopenfilename(title="Select .vgraph", filetypes=[("VGraph", "*.vgraph"), ("All", "*.*")])
        if p:
            self.vgraph_var.set(p)
            self.engine.vgraph = Path(p)

    def choose_import(self):
        p = filedialog.askopenfilename(
            title="Select OSM PBF / XML",
            filetypes=[("OSM PBF", "*.pbf"), ("OSM XML", "*.osm"), ("All", "*.*")])
        if p:
            self.import_async(p)

    def set_city(self, which: str):
        lat, lon = CITIES[self.city_var.get()]
        if which == "start":
            self.start_lat.set(f"{lat:.5f}")
            self.start_lon.set(f"{lon:.5f}")
            self.map.start = (lon, lat)
        else:
            self.goal_lat.set(f"{lat:.5f}")
            self.goal_lon.set(f"{lon:.5f}")
            self.map.goal = (lon, lat)
        self.map.redraw()

    def sync_paths(self):
        self.engine.db = Path(self.db_var.get())
        self.engine.vgraph = Path(self.vgraph_var.get())

    # -- import / prepare -------------------------------------------------
    def import_async(self, path: str):
        if self.busy:
            return
        self.sync_paths()
        self.set_busy(True, "Importing")
        self.log("hi", f"importing {Path(path).name} ...")
        flag = "--import-pbf" if path.lower().endswith(".pbf") else "--import-osm"

        def worker():
            try:
                rc = self.engine.stream(["--db", str(self.engine.db), flag, path],
                                        lambda line: self.post(self.log, "info", line))
                msg = "import complete" if rc == 0 else f"import failed (rc={rc})"
                self.post(self.log, "ok" if rc == 0 else "err", msg)
            except Exception as exc:
                self.post(self.log, "err", str(exc))
            finally:
                self.post(self.set_busy, False, "Ready")

        threading.Thread(target=worker, daemon=True).start()

    def prepare_async(self):
        if self.busy:
            return
        self.sync_paths()
        self.set_busy(True, "Preprocessing")
        self.log("hi", "preprocessing graph -> .vgraph ...")

        def worker():
            try:
                rc = self.engine.stream(
                    ["--db", str(self.engine.db), "--preprocess", "--save-bin", str(self.engine.vgraph)],
                    lambda line: self.post(self.log, "info", line))
                msg = "preprocess complete" if rc == 0 else f"preprocess failed (rc={rc})"
                self.post(self.log, "ok" if rc == 0 else "err", msg)
                if rc == 0:
                    self.post(self.load_map_async)
            except Exception as exc:
                self.post(self.log, "err", str(exc))
            finally:
                self.post(self.set_busy, False, "Ready")

        threading.Thread(target=worker, daemon=True).start()

    def load_map_async(self):
        if self.busy:
            return
        self.sync_paths()
        if not self.engine.vgraph.is_file():
            self.log("err", "no .vgraph cache found; prepare one first")
            return
        self.set_busy(True, "Loading map")
        self.log("hi", f"loading {self.engine.vgraph.name} ...")

        def worker():
            try:
                data = MapData.load(self.engine.vgraph)
                self.post(self.on_map_loaded, data)
            except Exception as exc:
                self.post(self.log, "err", f"map load failed: {exc}")
            finally:
                self.post(self.set_busy, False, "Ready")

        threading.Thread(target=worker, daemon=True).start()

    def on_map_loaded(self, data: MapData):
        self.map.data = data
        self.map.route = []
        self.map.start = None
        self.map.goal = None
        self.map.fit_data()
        self.map_status.set(f"{data.node_count:,} nodes  |  {data.edge_count:,} edges")
        self.log("ok", f"map loaded: {data.node_count:,} nodes, {data.edge_count:,} edges")
        if data.raw_node_count:
            self.log("info", f"raw {data.raw_node_count:,} nodes -> {data.contracted:,} contracted, "
                             f"{data.peeled:,} peeled")

    # -- benchmark --------------------------------------------------------
    def read_coords(self) -> tuple[float, float, float, float]:
        return (float(self.start_lat.get()), float(self.start_lon.get()),
                float(self.goal_lat.get()), float(self.goal_lon.get()))

    def run_benchmark(self):
        if self.busy:
            return
        try:
            lat1, lon1, lat2, lon2 = self.read_coords()
        except ValueError:
            messagebox.showerror(APP_TITLE, "Start and goal must be decimal coordinates.")
            return
        if not self.engine.available():
            messagebox.showerror(APP_TITLE, "uk_vector_router executable not found.")
            return
        self.sync_paths()
        if not self.engine.vgraph.is_file():
            self.log("err", "no .vgraph cache; routing from the DB (slow). Prepare one first.")
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.score_var.set("running...")
        self.set_busy(True, "Routing")
        self.log("hi", f"benchmark {lat1:.5f},{lon1:.5f} -> {lat2:.5f},{lon2:.5f} "
                       f"({', '.join(m for m, _ in MODES)})")

        tmp = Path(tempfile.gettempdir()) / f"uvm_route_{int(time.time())}.geojson"
        results: list[RouteResult] = []
        profile = self.profile_var.get()

        def worker():
            try:
                for mode, _label in MODES:
                    gj = tmp if mode == "bidir" else None
                    res, text, _elapsed = self.engine.route(mode, profile,
                                                            lat1, lon1, lat2, lon2, str(gj) if gj else None)
                    if gj and gj.is_file():
                        try:
                            res.coords = parse_geojson(gj)
                        except Exception as exc:
                            self.post(self.log, "err", f"GeoJSON parse failed: {exc}")
                    results.append(res)
                    self.post(self.on_result, res)
                if tmp.is_file():
                    tmp.unlink()
            except Exception as exc:
                self.post(self.log, "err", str(exc))
            finally:
                self.post(self.on_benchmark_done, results)

        threading.Thread(target=worker, daemon=True).start()

    def on_result(self, res: RouteResult):
        mode_label = dict(MODES).get(res.mode, res.mode)
        if res.found:
            self.tree.insert("", "end", values=(
                mode_label, "YES", f"{res.dist_km:.1f}", f"{res.est_min:.0f}",
                res.hops, f"{res.expanded:,}", f"{res.compute_ms:.1f}",
                f"{res.total_ms:.0f}", f"{res.score:.0f}"))
            self.log("ok", f"{res.mode}: {res.dist_km:.1f} km, {res.est_min:.0f} min, "
                           f"{res.hops} hops, {res.expanded:,} expansions, "
                           f"{res.compute_ms:.1f} ms compute, score {res.score:.0f}")
        else:
            self.tree.insert("", "end", values=(mode_label, "NO", "-", "-", "-", "-", "-", "-", "0"))
            self.log("err", f"{res.mode}: {res.error or 'no route'}")
        if res.coords:
            self.map.start = (res.coords[0][0], res.coords[0][1])
            self.map.goal = (res.coords[-1][0], res.coords[-1][1])
            self.map.set_route(res.coords)

    def on_benchmark_done(self, results: list[RouteResult]):
        self.set_busy(False, "Ready")
        found = [r for r in results if r.found]
        if found:
            best = max(found, key=lambda r: r.score)
            self.score_var.set(f"best: {dict(MODES).get(best.mode, best.mode)}  |  score {best.score:.0f}")
            for i, item in enumerate(self.tree.get_children()):
                if i < len(results) and results[i] is best:
                    self.tree.item(item, tags=("best",))
        else:
            self.score_var.set("no route")
        if results and results[0].coords:
            self.map.set_route(results[0].coords)
        self.log("info", "benchmark finished")

    def export_geojson(self):
        if not self.map.route:
            messagebox.showinfo(APP_TITLE, "Run a route first.")
            return
        p = filedialog.asksaveasfilename(defaultextension=".geojson", filetypes=[("GeoJSON", "*.geojson")])
        if not p:
            return
        coords = [[lon, lat] for lon, lat in self.map.route]
        data = {"type": "FeatureCollection", "features": [{
            "type": "Feature", "properties": {"source": "uk_vector_router"},
            "geometry": {"type": "LineString", "coordinates": coords}}]}
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(data, fh)
        self.log("ok", f"exported {p}")

    # -- misc -------------------------------------------------------------
    def bind_shortcuts(self):
        self.root.bind("<Control-r>", lambda _e: self.run_benchmark())
        self.root.bind("<Control-o>", lambda _e: self.choose_import())
        self.root.bind("<F5>", lambda _e: self.prepare_async())
        self.root.bind("<Control-l>", lambda _e: self.load_map_async())

    def on_close(self):
        self.root.destroy()


# ---------------------------------------------------------------------------
# Tests / entry point
# ---------------------------------------------------------------------------

def self_test() -> int:
    sample = """
Start Node  : osm=107774 (51.50755, -0.12783)
Goal Node   : osm=3394046210 (55.95345, -3.18852)
Found       : YES
Road Dist   : 607.97 km (direct TM chord: 533.47 km, detour: 1.140x)
Est Time    : 373.3 min (22396 s)
Performance : 806 hops | 462519 nodes expanded | 198.000 ms
"""
    p = parse_route_output(sample)
    assert p["found"] and p["dist_km"] == 607.97 and p["direct_km"] == 533.47
    assert p["est_min"] == 373.3 and p["hops"] == 806 and p["expanded"] == 462519
    assert p["compute_ms"] == 198.0
    s = score_route(p, 533.47)
    assert 0 < s <= 100, s
    gj = {"type": "FeatureCollection", "features": [{"geometry": {"type": "LineString", "coordinates": [[0, 51], [1, 52]]}}]}
    with tempfile.NamedTemporaryFile("w", suffix=".geojson", delete=False) as fh:
        json.dump(gj, fh)
        name = fh.name
    assert parse_geojson(name) == [(0.0, 51.0), (1.0, 52.0)]
    os.unlink(name)
    print("self-test OK")
    return 0


def smoke_test(exe: Optional[Path], db: Path, vgraph: Path) -> int:
    root = tk.Tk()
    engine = Engine(exe, db, vgraph)
    app = VectorRouterApp(root, engine, map_data=synthetic_map())
    root.update()
    app.canvas.event_generate("<Configure>", width=900, height=620)
    root.update()
    app.map.center_lon, app.map.center_lat = -1.0, 57.0
    app.map.k = 14.0
    app.map.redraw()
    root.update()
    app.map.fit_data()
    root.update()
    root.destroy()
    print("smoke-test OK")
    return 0


def route_test(exe: Optional[Path], db: Path, vgraph: Path) -> int:
    engine = Engine(exe, db, vgraph)
    if not engine.available():
        print("engine not found")
        return 1
    res, text, _ = engine.route("bidir", "balanced", 51.5074, -0.1278, 55.9533, -3.1883)
    print(text)
    print(f"parsed: found={res.found} dist={res.dist_km:.1f}km eta={res.est_min:.0f}min "
          f"expanded={res.expanded} compute={res.compute_ms:.1f}ms score={res.score:.0f}")
    return 0 if res.found else 1


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=APP_TITLE)
    ap.add_argument("--exe")
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--vgraph", default=str(DEFAULT_VGRAPH))
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--smoke-test", action="store_true")
    ap.add_argument("--route-test", action="store_true")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()
    if args.smoke_test:
        return smoke_test(find_executable(args.exe), Path(args.db), Path(args.vgraph))
    if args.route_test:
        return route_test(find_executable(args.exe), Path(args.db), Path(args.vgraph))

    root = tk.Tk()
    VectorRouterApp(root, Engine(find_executable(args.exe), Path(args.db), Path(args.vgraph)))
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())