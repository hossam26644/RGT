"""
This code is adapted from https://github.com/MoncktonLab/MoncktonWaterfall
which is adapted from:  https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/waterfall.py
Credit to: Michael Wood and Nathalie Ridgillova
"""

import re
import hashlib
import colorsys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import OrderedDict
from matplotlib.colors import to_rgb
from pathlib import Path

# ── constants ────────────────────────────────────────────────────────────────
BLANK   = (1.0, 1.0, 1.0)
UNKNOWN = (0.75, 0.75, 0.75)

_DEFAULT_COLORS: dict[str, str] = {
    "CAG": "#DB6968",
    "CCG": "#4D97CD",
    "CAA": "#46C1BE",
    "CCA": "#E8C559",
    "CCT": "#B279B4",
    "CTG": "#DB6968",
    "CTC": "#46C1BE",
    "CAC": "#46C1BE",
    "CAT": "#E8C559",
}

def _hash_color(unit: str) -> tuple[float, float, float]:
    """Stable, process-independent color derived from repeat unit string."""
    h = int(hashlib.md5(unit.encode()).hexdigest(), 16)
    hue = (h % 3600) / 3600.0   # finer granularity than mod 360
    return colorsys.hsv_to_rgb(hue, 0.70, 0.85)

def _build_color_map(user_units: list[str]) -> OrderedDict[str, tuple]:
    """
    Merge default units with user-supplied units.
    Default units keep fixed colors; novel units get hash-derived colors.
    Order: defaults first (insertion-ordered), then novel user units.
    """
    cm: OrderedDict[str, tuple] = OrderedDict()
    for unit, hex_col in _DEFAULT_COLORS.items():
        cm[unit] = to_rgb(hex_col)
    for unit in user_units:
        if unit not in cm:
            cm[unit] = _hash_color(unit)
    return cm

# ── raster construction ───────────────────────────────────────────────────────
def _motif_raster(
    reads: list[str],
    color_map: OrderedDict[str, tuple],
) -> np.ndarray:
    """
    Build (n_reads × max_len × 3) float32 RGB raster.
    - Matched motif positions → motif color
    - Unmatched positions within read → UNKNOWN
    - Positions beyond read end → BLANK
    Motifs sorted longest-first to prevent short motif masking longer matches.
    """
    n      = len(reads)
    width  = max(len(r) for r in reads)
    raster = np.ones((n, width, 3), dtype=np.float32)  # BLANK everywhere

    motifs  = list(color_map.keys())
    motifs_sorted = sorted(motifs, key=len, reverse=True)
    pattern = re.compile('|'.join(f'({re.escape(m)})' for m in motifs_sorted))

    for i, seq in enumerate(reads):
        slen = len(seq)
        # mark entire read extent as UNKNOWN first
        raster[i, :slen] = UNKNOWN
        # paint matched motifs
        for match in pattern.finditer(seq):
            color = color_map[match.group()]
            raster[i, match.start():match.end()] = color

    return raster


# ── public API ────────────────────────────────────────────────────────────────
def draw_waterfall(
    reads:        list[str],
    settings:     dict,
    title: str,
    export_directory: str,
) -> np.ndarray:
    """
    Parameters
    ----------
    reads             : list of forward-strand sequences (str)
    settings          : user settings dict
    title             : string sampleID
    export_directory  : path string for output

    """
    if not reads:
        raise ValueError("reads list is empty")

    #sort reads by length
    reads = sorted(reads, key=len, reverse=True)
    user_units  = settings["repeat_units"]
    color_map   = _build_color_map(user_units)

    # retain only units that actually appear in at least one read (keeps legend clean)
    present = {u for seq in reads for u in color_map if u in seq}
    color_map = OrderedDict((u, c) for u, c in color_map.items() if u in present)

    raster = _motif_raster(reads, color_map)

    fig, ax = plt.subplots(figsize=(10, max(3, len(raster) * 0.15)))
    ax.imshow(raster, origin='lower', aspect='auto', interpolation='nearest')
    ax.set_xlabel('Position')
    ax.set_ylabel('Reads')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(title)

    patches = [
        mpatches.Patch(color=color, label=motif)
        for motif, color in color_map.items()
    ]
    ax.legend(handles=patches, bbox_to_anchor=(1.02, 0.6),
              loc='upper left', frameon=False, fontsize=7)

    plt.tight_layout()
    out = export_directory
    fig.savefig(out, dpi=400, format="png", bbox_inches='tight')
    plt.close(fig)

