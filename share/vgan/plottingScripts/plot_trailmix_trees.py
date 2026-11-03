#!/usr/bin/env python3
"""
plot_trailmix_trees.py  --  publication-style figures for TrailMix results.

Figures (per MCMC chain, plus a copy of the best-mixing chain as "chainBest"):

    {prefix}_chainN_figure.pdf       multi-panel summary: (a) tree with posterior
                                     placements, (b) mixture proportions,
                                     (c) placement posterior, (d) trace,
                                     (e) agreement between chains
    {prefix}_chainN_combined.pdf     the tree panel alone, all sources
    {prefix}_chainN_Source_K.pdf     the tree panel for one source
    {prefix}_chainN_proportions.pdf  posterior of the mixture proportions
and once per run:
    {prefix}_chains.pdf              do the chains agree? (proportions + placements)
    {prefix}_rpvg_*.pdf              RPVG's initial (pre-MCMC) estimate, as an
                                     aligned bar track on the tree + bar charts

How to read the tree
--------------------
Rectangular phylogram: right-angle branches, so branches cannot cross. Each
source's posterior is drawn ON the branch it landed on:

    * the swelling along the branch = posterior density of WHERE along the
      branch the source sits (thicker = more MCMC time spent there);
    * the dot = median position; its AREA = posterior probability of that
      branch (share of MCMC iterations spent on it);
    * the best-supported tip is set in bold with a color swatch; when the
      best placement is an internal branch its whole clade is shaded;
    * sources sharing a branch are stacked slightly apart so none hides another.

Usage
-----
    plot_trailmix_trees.py RUN_PREFIX
        Auto mode. Reads {RUN_PREFIX}summary.txt (written by vgan trailmix) and
        takes the database, chains and RPVG files from it. Draws
        haps_rerooted.treefile from the database directory when it exists
        (otherwise haps.treefile).

    plot_trailmix_trees.py TREE --multi --run-prefix RUN_PREFIX [--out-prefix P]
    plot_trailmix_trees.py TREE MCMC.tsv[.gz] [OUT_PREFIX]           (one chain)

Trees
-----
    --db DB          TrailMix database dir (--tm-files) or its haps.treefile --
                     the tree the MCMC placements refer to. TREE may be a
                     differently ROOTED display tree. Branches are matched by
                     the tip split they induce, never by name, so placements
                     land on the right branch whatever the rooting. (Node names
                     in a rerooted file are NOT reliable on the path between
                     the old and the new root, and the new root may sit in the
                     middle of a branch; both are handled.)
    --db-label TEXT  database name printed in the figure footer (default: directory name
                     without a trailing _nohh; "" omits it)
    --no-rerooted    auto mode: draw haps.treefile instead of haps_rerooted.treefile
    --layout L       phylogram (default; true branch lengths) | unit (every
                     branch length 1, for trees whose long branches hide the
                     interesting part) | cladogram (topology only, tips aligned)
    --break-long X   auto (default) | off | LENGTH. Over-long branches are drawn
                     shortened, marked //, but only when that makes the tree at
                     least 15% narrower.
    --keep LIST      comma list (or file, one per line) of tips / clades to keep;
                     everything else is pruned. Placed clades are always kept.
    --ladderize      sort clades by size.   --intern   print internal node names.

Output
------
    --png            also write 300 dpi PNGs
    --best-only      draw only the best chain (as chainBest) and the comparison
    --slot-order     keep vgan's Source_1, Source_2 ... slot numbers. By default sources
                     are numbered by decreasing mean proportion in the drawn chain
                     (Source 1 = largest), because the sampler may swap the slots
                     and this keeps colors and numbers the same in every panel,
                     in chains.pdf, and across runs.
    --min-prob P     do not draw branches below this posterior probability (default 0:
                     every branch the chain visited is drawn)
    --annotate       write probability / node name next to placement markers
    --pitch-mm MM    row height per tip (default 2.55)
    --rpvg-prefix P  prefix of {P}RPVGHaplotypeTranscripts.txt / RPVGHaplotypes.txt

A chain file that is still being written (truncated .gz) is read up to its last
complete row, so a run can be plotted while it is in progress.

Style: Arial/Helvetica-compatible sans serif (Liberation Sans / Nimbus Sans when
Arial is missing) at 5-7 pt, Nature column widths (89 / 120 / 183 mm), hairline
rules, fonts embedded as TrueType (editable text), color-blind-safe palette
(blue / orange / aqua), scale bar instead of a numeric axis.
"""

import sys
import os
import re
import glob
import shutil
import argparse
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from Bio import Phylo
import matplotlib as mpl

mpl.use("Agg")  # file output only; never needs a display
import matplotlib.pyplot as plt
from matplotlib import font_manager, patheffects
from matplotlib.collections import LineCollection
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
from matplotlib.textpath import TextPath
from scipy.ndimage import gaussian_filter1d

warnings.filterwarnings("ignore", message=".*Glyph.*missing.*")


# --------------------------------------------------------------------------- #
#  Palette
# --------------------------------------------------------------------------- #
# Categorical slots 1-3 of the validated 8-hue theme: every pair clears the
# color-blind (deutan/protan/tritan) and normal-vision separation checks in
# light mode. TrailMix supports at most three sources, so no hue is reused.
SOURCE_COLORS = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a"}
SOURCE_FALLBACK = "#898781"

# Text never wears a series color (the aqua slot is below 3:1 on white);
# identity is carried by swatches/markers next to ink-colored text.
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#898781"
RULE = "#c3c2b7"
GRID = "#e1e0d9"
SURFACE = "#ffffff"
TREE_INK = "#3f3e3c"
TREE_LEADER = "#c9c8bf"
HILITE_INTERNAL = "#9c4a3f"

# Sequential single-hue ramp for magnitudes (RPVG TPM / read counts).
_SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
        "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SEQ_CMAP = mpl.colors.LinearSegmentedColormap.from_list("seq_blue", _SEQ)

# Chains are an ordinal nuisance factor, not an identity: neutral grays.
CHAIN_GRAYS = ["#2b2a28", "#5b5a56", "#85847e", "#a9a8a1", "#c9c8bf", "#dedcd5"]


ANNOTATE = False     # set from --annotate
DB_LABEL = None      # set from --db-label (None: derive from the database directory name)


def src_color(i):
    return SOURCE_COLORS.get(int(i), SOURCE_FALLBACK)


def seq_color(v, vmax, floor=0.2):
    t = 1.0 if vmax <= 0 else min(max(v / vmax, 0.0), 1.0)
    return SEQ_CMAP(floor + (1.0 - floor) * t)


# --------------------------------------------------------------------------- #
#  Sizes / typography (Nature: 89 mm one column, 120 mm 1.5, 183 mm two)
# --------------------------------------------------------------------------- #
MM = 1.0 / 25.4
COL1, COL15, COL2 = 89.0, 120.0, 183.0      # mm
MAX_H = 247.0                               # mm, full-page ceiling

FS_TIP = 5.0        # tip labels
FS_TXT = 6.0        # annotations, tick labels
FS_AX = 6.5         # axis titles
FS_PANEL = 8.5      # panel letters (bold, lower case)
FONT = "DejaVu Sans"


def _pick_font():
    have = {f.name for f in font_manager.fontManager.ttflist}
    # Arial/Helvetica are what journals ask for; Liberation Sans and Nimbus
    # Sans are metric-compatible stand-ins that ship with most Linux hosts.
    for cand in ("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans",
                 "FreeSans", "DejaVu Sans"):
        if cand in have:
            return cand
    return "DejaVu Sans"


def setup_style():
    global FONT
    FONT = _pick_font()
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT],
        "font.size": FS_TXT,
        "axes.titlesize": FS_AX,
        "axes.labelsize": FS_AX,
        "axes.labelcolor": INK,
        "axes.labelpad": 2.0,
        "axes.linewidth": 0.5,
        "axes.edgecolor": RULE,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "text.color": INK,
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "xtick.labelsize": FS_TXT,
        "ytick.labelsize": FS_TXT,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.0,
        "ytick.major.size": 2.0,
        "xtick.major.pad": 2.0,
        "ytick.major.pad": 2.0,
        "legend.frameon": False,
        "legend.fontsize": FS_TXT,
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "lines.solid_capstyle": "butt",
    })


def text_w_mm(s, size=FS_TXT, weight="normal"):
    """Rendered width of a string in mm (so labels can be laid out exactly)."""
    if not s:
        return 0.0
    fp = FontProperties(family=FONT, size=size, weight=weight)
    return TextPath((0, 0), s, prop=fp).get_extents().width * 25.4 / 72.0


HALO = [patheffects.withStroke(linewidth=1.6, foreground=SURFACE)]


def new_fig(w_mm, h_mm):
    return plt.figure(figsize=(w_mm * MM, h_mm * MM))


def add_axes_mm(fig, left, bottom, width, height, **kw):
    """Axes placed by absolute mm from the figure's lower-left corner."""
    W, H = fig.get_size_inches() * 25.4
    return fig.add_axes([left / W, bottom / H, width / W, height / H], **kw)


def panel_letter(fig, letter, x_mm, y_mm_from_top):
    W, H = fig.get_size_inches() * 25.4
    fig.text(x_mm / W, 1.0 - y_mm_from_top / H, letter, fontsize=FS_PANEL,
             fontweight="bold", ha="left", va="top", color=INK)


def style_axes(ax, left=False):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(True)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(RULE)
        ax.spines[s].set_linewidth(0.5)
    ax.tick_params(colors=INK_2, width=0.5)


def save_fig(fig, path, png=False):
    fig.savefig(path)
    outs = [path]
    if png:
        p = re.sub(r"\.pdf$", ".png", path)
        fig.savefig(p, dpi=300)
        outs.append(p)
    plt.close(fig)
    for o in outs:
        print(f"  Saved: {o}")
    return path


# --------------------------------------------------------------------------- #
#  Tree model
# --------------------------------------------------------------------------- #

class Node:
    """Light tree node (so pruning/collapsing never fights Biopython)."""
    __slots__ = ("name", "length", "children", "parent", "x", "y", "broken")

    def __init__(self, name=None, length=0.0):
        self.name = name.strip() if isinstance(name, str) and name.strip() else None
        self.length = float(length) if length else 0.0
        self.children = []
        self.parent = None
        self.x = 0.0
        self.y = 0.0
        self.broken = False

    @property
    def is_tip(self):
        return not self.children

    def preorder(self):
        stack = [self]
        while stack:
            n = stack.pop()
            yield n
            stack.extend(reversed(n.children))

    def postorder(self):
        out = list(self.preorder())
        return reversed(out)

    def tips(self):
        return [n for n in self.preorder() if n.is_tip]


def _from_bio(clade, parent=None):
    n = Node(clade.name, clade.branch_length)
    n.parent = parent
    n.children = [_from_bio(c, n) for c in clade.clades]
    return n


def read_tree(path):
    return _from_bio(Phylo.read(path, "newick").root)


def tip_sets(root):
    """node -> frozenset of the tip names below it."""
    ts = {}
    for n in root.postorder():
        ts[n] = frozenset([n.name]) if n.is_tip else frozenset().union(*[ts[c] for c in n.children])
    return ts


def build_edge_map(root):
    """Placement coordinates on the plotted tree itself: a placement named X
    sits on X's own parent branch, t measured parent -> X.

    edge_map[name] = [(t0, t1, node, a0, a1), ...]: for t in [t0, t1] the
    point is on the parent branch of `node`, a fraction a0 + (t-t0)/(t1-t0)
    * (a1-a0) of the way from node's parent to node. clade_tips[name] = the
    tips below X in the DATABASE orientation."""
    ts = tip_sets(root)
    em, ct = {}, {}
    for n in root.preorder():
        if n.name:
            em[n.name] = [(0.0, 1.0, n, 0.0, 1.0)]
            ct[n.name] = ts[n]
    return em, ct


def map_db_edges(db_root, disp_root):
    """
    MCMC placements name branches of the DATABASE tree (haps.treefile). The
    display tree may be the same tree rooted elsewhere (haps_rerooted.treefile),
    in which case node names on the path between the two roots are NOT
    reliable (a name stays on its node while the branch above that node
    changes, and the new root may sit in the middle of a branch). So match
    branches by what they are -- the split of tips they induce -- not by name.

    Returns (edge_map, clade_tips), or None if the two trees do not have the
    same tips (caller falls back to name matching).
    """
    db_ts, dp_ts = tip_sets(db_root), tip_sets(disp_root)
    all_tips = frozenset(t.name for t in disp_root.tips())
    if all_tips != frozenset(t.name for t in db_root.tips()):
        return None
    ref = min(all_tips)

    def key(S):                     # canonical side of an unrooted split
        return S if ref not in S else all_tips - S

    disp_edges = {}
    for d in disp_root.preorder():
        if d.parent is not None:
            disp_edges.setdefault(key(dp_ts[d]), []).append(d)

    # internal names of the display file are unreliable: drop and re-derive
    for d in disp_root.preorder():
        if not d.is_tip:
            d.name = None

    em, ct = {}, {}
    n_unmapped = 0
    for c in db_root.preorder():
        if not c.name:
            continue
        S = db_ts[c]
        ct[c.name] = S
        if c.parent is None:
            continue
        lst = disp_edges.get(key(S))
        if not lst:
            n_unmapped += 1
            continue
        if len(lst) == 1:
            d = lst[0]
            same = dp_ts[d] == S
            em[c.name] = [(0.0, 1.0, d, 0.0, 1.0) if same else (0.0, 1.0, d, 1.0, 0.0)]
            if same and not c.is_tip:
                d.name = c.name            # safe to label: same node, same side
        elif len(lst) == 2 and lst[0].parent is lst[1].parent:
            A, B = lst                     # the display root sits on this branch
            child, par = (A, B) if dp_ts[A] == S else (B, A)
            lc, lp = max(child.length, 0.0), max(par.length, 0.0)
            frac = lp / (lc + lp) if (lc + lp) > 0 else 0.5
            segs = []
            if frac > 0:                   # DB parent end -> display root, along `par`, reversed
                segs.append((0.0, frac, par, 1.0, 0.0))
            if frac < 1:                   # display root -> DB child end, along `child`
                segs.append((frac, 1.0, child, 0.0, 1.0))
            em[c.name] = segs
        else:
            n_unmapped += 1
    if n_unmapped:
        print(f"  Note: {n_unmapped} database branch(es) have no counterpart in the display tree.",
              file=sys.stderr)
    return em, ct


def recover_labels(display_root, orig_root):
    """Fallback when the display tree and the database tree have different
    tips: stamp lost internal names back by identical tip set."""
    ts_d, ts_o = tip_sets(display_root), tip_sets(orig_root)
    by_tips = {ts_d[n]: n for n in display_root.preorder()}
    names = {n.name for n in display_root.preorder() if n.name}
    got = []
    for n in orig_root.preorder():
        if not n.name or n.name in names:
            continue
        tgt = by_tips.get(ts_o[n])
        if tgt is not None and not tgt.name:
            tgt.name = n.name
            got.append(n.name)
    if got:
        print(f"Recovered {len(got)} internal label(s) by tip set: {', '.join(sorted(got))}",
              file=sys.stderr)


def resolve_db_tree(db_arg):
    cand = os.path.join(db_arg, "haps.treefile") if os.path.isdir(db_arg) else db_arg
    if not os.path.isfile(cand):
        print(f"Warning: --db tree not found: {cand}", file=sys.stderr)
        return None
    return cand


def _collapse_unary(root, edge_map):
    """Merge single-child internal nodes (left behind by pruning) and keep
    every placement coordinate valid: a position on a merged-away branch is
    re-expressed as a position on the merged branch."""
    for n in list(root.postorder()):
        if n.is_tip or len(n.children) != 1 or n is root:
            continue
        y = n.children[0]
        tot = n.length + y.length
        a = (n.length / tot) if tot > 0 else 0.5
        for name, segs in list(edge_map.items()):
            new = []
            for (t0, t1, tgt, c0, c1) in segs:
                if tgt is n:
                    new.append((t0, t1, y, c0 * a, c1 * a))
                elif tgt is y:
                    new.append((t0, t1, y, a + (1 - a) * c0, a + (1 - a) * c1))
                else:
                    new.append((t0, t1, tgt, c0, c1))
            edge_map[name] = new
        y.length = tot
        p = n.parent
        p.children[p.children.index(n)] = y
        y.parent = p
    while len(root.children) == 1:          # a unary root: promote its only child
        y = root.children[0]
        for name, segs in list(edge_map.items()):
            edge_map[name] = [sg for sg in segs if sg[2] is not y]
        y.parent = None
        y.length = 0.0
        root = y
    return root


def prune_tree(root, keep_tips, edge_map):
    def rec(n):
        n.children = [c for c in n.children if rec(c)]
        for c in n.children:
            c.parent = n
        return bool(n.children) or (n.is_tip and n.name in keep_tips)

    rec(root)
    root = _collapse_unary(root, edge_map)
    alive = {id(n) for n in root.preorder()}
    for name, segs in list(edge_map.items()):
        segs = [sg for sg in segs if id(sg[2]) in alive and sg[2].parent is not None]
        if segs:
            edge_map[name] = segs
        else:
            edge_map.pop(name)
    return root


def ladderize(root):
    size = {}
    for n in root.postorder():
        size[id(n)] = 1 if n.is_tip else sum(size[id(c)] for c in n.children)
    for n in root.preorder():
        n.children.sort(key=lambda c: size[id(c)])


@dataclass
class Layout:
    root: Node
    tips: list
    nodes: list
    by_name: dict
    edge_map: dict               # name -> [(t0, t1, node, a0, a1), ...]  (see build_edge_map)
    clade_tips: dict             # name -> tips below it (database orientation)
    mode: str
    xmax: float                  # x of the right-most tip (tree width)
    cap: float = None
    n_broken: int = 0


def _auto_cap(root):
    """Truncate over-long branches only when it genuinely buys resolution:
    the layout must get at least 15% narrower. Otherwise draw true lengths."""
    lens = np.array([n.length for n in root.preorder() if n is not root and n.length > 0])
    if lens.size < 8:
        return None
    cap = max(3.0 * np.quantile(lens, 0.90), 1e-12)
    if lens.max() <= cap:
        return None
    full = max(_true_depth(root).values())
    cut = max(d for d in _cut_depth(root, cap).values())
    return cap if cut < 0.85 * full else None


def _cut_depth(root, cap):
    d = {root: 0.0}
    for n in root.preorder():
        for c in n.children:
            L = c.length
            if L > cap:
                L = cap + 0.08 * (L - cap)
            d[c] = d[n] + L
    return d


def _true_depth(root):
    d = {root: 0.0}
    for n in root.preorder():
        for c in n.children:
            d[c] = d[n] + c.length
    return d


def compute_layout(root, edge_map, clade_tips, mode="phylogram", break_long="auto"):
    tips = root.tips()
    for i, t in enumerate(tips):
        t.y = float(i)
    for n in root.postorder():
        if not n.is_tip:
            n.y = 0.5 * (n.children[0].y + n.children[-1].y)

    cap, n_broken = None, 0
    if mode == "phylogram":
        if break_long == "auto":
            cap = _auto_cap(root)
        elif break_long not in ("off", None):
            cap = float(break_long)
        root.x = 0.0
        for n in root.preorder():
            for c in n.children:
                L = c.length
                c.broken = False
                if cap is not None and L > cap:
                    L = cap + 0.08 * (L - cap)
                    c.broken = True
                    n_broken += 1
                c.x = n.x + L
    elif mode == "unit":
        root.x = 0.0
        for n in root.preorder():
            for c in n.children:
                c.x = n.x + 1.0
                c.broken = False
    elif mode == "cladogram":
        h = {}
        for n in root.postorder():
            h[n] = 0 if n.is_tip else 1 + max(h[c] for c in n.children)
        H = h[root]
        for n in root.preorder():
            n.x = float(H - h[n])
            n.broken = False
    else:
        raise ValueError(mode)

    nodes = list(root.preorder())
    by_name = {n.name: n for n in nodes if n.name}
    xmax = max(t.x for t in tips)
    return Layout(root, tips, nodes, by_name, edge_map, clade_tips, mode, xmax, cap, n_broken)


def _seg_at(lay, name, t):
    segs = lay.edge_map[name]
    for sg in segs:
        if sg[0] <= t <= sg[1]:
            return sg
    return segs[-1] if t > segs[-1][1] else segs[0]


def _seg_x(sg, tt_arr):
    """x for fractions along a segment (tt measured parent -> node)."""
    node = sg[2]
    xp = node.parent.x
    return xp + (node.x - xp) * tt_arr


def place_xy(lay, name, t):
    """(x, y) of position t (scalar) on the branch carrying placement `name`."""
    t0, t1, node, a0, a1 = _seg_at(lay, name, t)
    u = 0.0 if t1 <= t0 else (t - t0) / (t1 - t0)
    return float(_seg_x((t0, t1, node, a0, a1), a0 + u * (a1 - a0))), node.y


# --------------------------------------------------------------------------- #
#  Posterior summaries
# --------------------------------------------------------------------------- #

@dataclass
class BranchPost:
    name: str
    p: float
    t: np.ndarray
    mapped: bool = True


@dataclass
class SourcePost:
    idx: int                           # slot number as written by vgan (1-based)
    label_idx: int                     # number shown on the figure
    branches: list = field(default_factory=list)   # by decreasing p
    prop: np.ndarray = None
    n: int = 0

    @property
    def top(self):
        return self.branches[0] if self.branches else None

    @property
    def color(self):
        return src_color(self.label_idx)

    @property
    def name(self):
        return self.top.name if self.top else f"source {self.label_idx}"


def _match_tree_name(name, known):
    """MCMC rows carry the GBWT path name (e.g. 'HV2a1*1'); the tree tip is the sanitized
    name ('HV2a1_1'). vgan replaces + ' * @ ( ) by '_' and drops a trailing .1/.2, so try
    the same steps before declaring a placement 'absent from the plotted tree'."""
    if name in known:
        return name
    c = re.sub(r"[+'*@()]", "_", name)
    if c in known:
        return c
    c2 = re.sub(r"\.(1|2)$", "", c)
    return c2 if c2 in known else name


def load_mcmc(path):
    """Read a chain file. A chain that is still being written (truncated gzip)
    is read up to its last complete row, with a warning, rather than failing."""
    comp = "gzip" if str(path).endswith(".gz") else None
    try:
        return pd.read_csv(path, sep="\t", compression=comp)
    except (EOFError, OSError) as err:
        if comp is None:
            raise
        import gzip
        import io
        lines = []
        with gzip.open(path, "rt") as fh:
            try:
                for ln in fh:
                    lines.append(ln)
            except (EOFError, OSError):
                pass
        if lines and not lines[-1].endswith("\n"):
            lines.pop()
        if len(lines) < 3:
            raise
        print(f"  Warning: {path} is truncated ({err}); using its {len(lines) - 1:,} complete "
              f"rows. Is the run still in progress?", file=sys.stderr)
        return pd.read_csv(io.StringIO("".join(lines)), sep="\t")


def summarize_chain(mcmc, lay, sort_sources=False):
    slots = [i for i in (1, 2, 3)
             if f"Source_{i}" in mcmc.columns and f"branch_position_derived_{i}" in mcmc.columns]
    posts = []
    for i in slots:
        names = mcmc[f"Source_{i}"].astype(str).str.strip()
        # was: names = mcmc[...].astype(str).str.strip().to_numpy()  -- placements on 'HV2a1*1' (tree tip 'HV2a1_1') were
        # treated as absent from the tree and silently dropped
        _alias = {u: _match_tree_name(u, lay.edge_map) for u in names.unique()}
        names = names.map(_alias).to_numpy()
        t = pd.to_numeric(mcmc[f"branch_position_derived_{i}"], errors="coerce").to_numpy()
        ok = ~np.isnan(t)
        names, t = names[ok], t[ok]
        n = len(names)
        sp = SourcePost(idx=i, label_idx=i, n=n)
        pc = f"proportion_{i}"
        if pc in mcmc.columns:
            sp.prop = pd.to_numeric(mcmc[pc], errors="coerce").dropna().to_numpy()
        if n:
            order = pd.Series(np.arange(n)).groupby(names).apply(lambda s: s.to_numpy())
            for nm, idxs in order.items():
                sp.branches.append(BranchPost(nm, len(idxs) / n, t[idxs],
                                              mapped=nm in lay.edge_map and bool(lay.edge_map[nm])))
            sp.branches.sort(key=lambda b: -b.p)
        posts.append(sp)

    if sort_sources and len(posts) > 1:
        mean = lambda s: float(np.mean(s.prop)) if s.prop is not None and len(s.prop) else 0.0
        for rank, sp in enumerate(sorted(posts, key=lambda s: -mean(s)), start=1):
            sp.label_idx = rank
        posts.sort(key=lambda s: s.label_idx)

    lost = {}
    for sp in posts:
        for b in sp.branches:
            if not b.mapped:
                lost[b.name] = lost.get(b.name, 0.0) + b.p * sp.n
    tot_n = sum(sp.n for sp in posts) or 1
    if lost and sum(lost.values()) / tot_n >= 0.005:
        tot = sum(sp.n for sp in posts) or 1
        items = ", ".join(f"{k} ({v / tot:.1%})" for k, v in sorted(lost.items(), key=lambda kv: -kv[1])[:6])
        print(f"  Warning: placements on branches absent from the plotted tree "
              f"are not drawn: {items}.\n  (Wrong database tree for this run? Use --db "
              f"with the --tm-files directory the run used.)", file=sys.stderr)
    return posts


def ci(v, lo=2.5, hi=97.5):
    return float(np.percentile(v, lo)), float(np.percentile(v, hi))


# --------------------------------------------------------------------------- #
#  The tree panel
# --------------------------------------------------------------------------- #

@dataclass
class TreeGeom:
    ax: object
    lay: Layout
    pitch_mm: float
    ax_w_mm: float
    label_w_mm: float
    s: float = 1.0                     # mm per data unit (x)
    lab_x: float = 0.0                 # x of the label column (data)
    extra_rows: float = 0.0

    def dx(self, mm):                  # mm -> data units along x
        return mm / self.s


LABEL_GAP_MM = 2.2        # tip end -> label text
LEADER_GAP_MM = 0.5


def tip_label_width_mm(lay, extra=0.0):
    return max((text_w_mm(t.name or "", FS_TIP, "bold") for t in lay.tips), default=10.0) + extra


def setup_tree_axes(ax, lay, ax_w_mm, pitch_mm, label_w_mm, extra_rows, track_w_mm=0.0):
    """Fix data limits so that 1 row == pitch_mm and the label column fits."""
    W = lay.xmax if lay.xmax > 0 else 1.0
    stem = 0.025 * W
    xmin = -stem
    usable = ax_w_mm - label_w_mm - LABEL_GAP_MM - 0.8 - track_w_mm
    s = usable / (W - xmin)
    xmax = xmin + ax_w_mm / s
    n = len(lay.tips)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(n - 0.5 + extra_rows, -0.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.patch.set_alpha(0.0)
    g = TreeGeom(ax, lay, pitch_mm, ax_w_mm, label_w_mm, s=s,
                 lab_x=W + LABEL_GAP_MM / s, extra_rows=extra_rows)
    return g


def draw_tree_skeleton(g, internal_labels=False):
    ax, lay = g.ax, g.lay
    segs = []
    for n in lay.nodes:
        p = n.parent
        if p is None:
            segs.append([(-0.025 * lay.xmax, n.y), (n.x, n.y)])     # root stem
            continue
        segs.append([(p.x, p.y), (p.x, n.y), (n.x, n.y)])
    ax.add_collection(LineCollection(segs, colors=TREE_INK, linewidths=0.55,
                                     capstyle="butt", joinstyle="miter", zorder=2))

    # // marks on truncated long branches
    for n in lay.nodes:
        if n.broken and n.parent is not None:
            xm = 0.5 * (n.parent.x + n.x)
            d = g.dx(0.55)
            dy = 0.34
            ax.add_patch(Rectangle((xm - d * 1.5, n.y - 0.45), 3 * d, 0.9,
                                   facecolor=SURFACE, edgecolor="none", zorder=3))
            for sh in (-0.9, 0.9):
                ax.plot([xm + sh * d - d, xm + sh * d + d], [n.y + dy, n.y - dy],
                        color=TREE_INK, lw=0.5, zorder=4, solid_capstyle="butt")

    # aligned tip labels with dotted leader lines
    for t in lay.tips:
        if t.x < lay.xmax - g.dx(0.8):
            ax.plot([t.x + g.dx(LEADER_GAP_MM), g.lab_x - g.dx(0.1)], [t.y, t.y],
                    color=TREE_LEADER, lw=0.35, dashes=(1, 1.6), zorder=1)

    if internal_labels:
        for n in lay.nodes:
            if n.name and not n.is_tip:
                ax.text(n.x - g.dx(0.6), n.y - 0.2, n.name, fontsize=FS_TIP - 0.6,
                        ha="right", va="bottom", color=HILITE_INTERNAL, zorder=5,
                        path_effects=HALO)


def draw_tip_labels(g, hilite):
    """hilite: {tip name: color}. Highlighted tips get a swatch + bold ink."""
    ax = g.ax
    for t in g.lay.tips:
        c = hilite.get(t.name)
        x = g.lab_x + g.dx(1.5)
        if c:
            ax.add_patch(Rectangle((g.lab_x + g.dx(0.1), t.y - 0.3), g.dx(1.0), 0.6,
                                   facecolor=c, edgecolor="none", zorder=4))
        ax.text(x, t.y, t.name, fontsize=FS_TIP, va="center", ha="left",
                color=INK if c else INK_2, fontweight="bold" if c else "normal",
                zorder=5)


def nice_bar(target):
    if target <= 0:
        return 1.0
    mag = 10 ** np.floor(np.log10(target))
    for m in (1, 2, 5, 10):
        if m * mag >= target:
            return m * mag
    return mag


def draw_scale_bar(g):
    """Scale bar (left) and, right-aligned in the same strip, nothing else:
    keeps the bottom of the tree uncluttered. Returns the row it used."""
    lay, ax = g.lay, g.ax
    n = len(lay.tips)
    y = n - 0.5 + 1.3
    if lay.mode == "phylogram":
        bar = nice_bar(lay.xmax * 0.15)
        ax.plot([0, bar], [y, y], color=INK, lw=0.9, solid_capstyle="butt", clip_on=False)
        ax.text(bar + g.dx(1.5), y, f"{bar:g} substitutions / site", ha="left", va="center",
                fontsize=FS_TXT - 0.5, color=INK_2)
    else:
        what = "branch length = 1 (topology only)" if lay.mode == "unit" else "topology only, tips aligned"
        ax.text(0, y, what, ha="left", va="center", fontsize=FS_TXT - 0.5, color=INK_2)
    if lay.n_broken:
        ax.text(0, y + 1.7, "//  branch drawn shortened", ha="left", va="center",
                fontsize=FS_TXT - 0.5, color=INK_3)
    return y


# ----- placements ---------------------------------------------------------- #

def _hist_profile(t, nb=36):
    cnt, _ = np.histogram(np.clip(t, 0, 1), bins=nb, range=(0.0, 1.0))
    cnt = cnt.astype(float)
    k = np.array([1.0, 2.0, 3.0, 2.0, 1.0])
    k /= k.sum()
    pad = np.concatenate([np.full(2, cnt[0]), cnt, np.full(2, cnt[-1])])
    sm = np.convolve(pad, k, mode="valid")
    return sm  # same length as cnt (preserves mass, soft edges)


def _marker_area(p):
    """Marker area (pt^2) for a branch visited with posterior probability p."""
    return 2.5 + 36.5 * np.sqrt(p)


def draw_placements(g, posts, min_prob=0.0, hmax_rows=0.40, annotate=False):
    ax, lay = g.ax, g.lay
    shown = {}
    for sp in posts:
        shown[sp.idx] = [b for b in sp.branches if b.mapped and b.p >= min_prob]
    # Sources that sit on the same branch are stacked a little apart (dodged)
    # so no source hides another.
    users = {}
    for sp in posts:
        for b in shown[sp.idx]:
            users.setdefault(b.name, []).append(sp.idx)

    def dodge(name, idx):
        u = users[name]
        n = len(u)
        return (u.index(idx) - (n - 1) / 2.0) * 0.30 if n > 1 else 0.0, \
            hmax_rows / (1.0 + 0.45 * (n - 1))

    for sp in posts:
        branches = shown[sp.idx]
        if not branches:
            continue
        # spindle thickness is normalized per source: the busiest bin of the
        # busiest branch is full height, everything else is relative to it
        profiles = {b.name: _hist_profile(b.t) / max(len(b.t), 1) * b.p for b in branches}
        gmax = max(prof.max() for prof in profiles.values())
        for b in branches:
            prof = profiles[b.name]
            nb = len(prof)
            cen = (np.arange(nb) + 0.5) / nb
            dy, hmax = dodge(b.name, sp.idx)
            for sg in lay.edge_map[b.name]:
                t0, t1, node, a0, a1 = sg
                if t1 <= t0:
                    continue
                tg = np.linspace(t0, t1, 121)
                h = np.interp(tg, cen, prof) / gmax * hmax
                xs = _seg_x(sg, a0 + (tg - t0) / (t1 - t0) * (a1 - a0))
                y = node.y + dy
                poly = np.vstack([np.column_stack([xs, y - h]),
                                  np.column_stack([xs[::-1], (y + h)[::-1]])])
                ax.add_patch(Polygon(poly, closed=True, facecolor=sp.color, edgecolor="none",
                                     alpha=0.45, zorder=6, linewidth=0))
        # markers last, so they sit on top of every spindle
        for b in branches:
            xm, ym = place_xy(lay, b.name, float(np.median(b.t)))
            dy, _ = dodge(b.name, sp.idx)
            ax.scatter([xm], [ym + dy], s=_marker_area(b.p), color=sp.color,
                       edgecolors="white", linewidths=0.45, zorder=8)
        if annotate:
            for b in branches[:3]:
                if b.p < 0.10:
                    continue
                xm, ym = place_xy(lay, b.name, float(np.median(b.t)))
                txt = f"{b.p:.2f}" if b.p < 0.995 else ""
                if b.name in lay.clade_tips and len(lay.clade_tips[b.name]) > 1:
                    txt = b.name + ("  " + txt if txt else "")
                if txt:
                    ax.text(xm, ym - 0.62, txt, fontsize=FS_TIP - 0.4, ha="center",
                            va="bottom", color=INK_2, zorder=9, path_effects=HALO)


def _runs(ys):
    """Contiguous runs of integer row numbers."""
    ys = sorted(set(int(round(v)) for v in ys))
    out, start = [], None
    for i, v in enumerate(ys):
        if start is None:
            start = prev = v
        elif v == prev + 1:
            prev = v
        else:
            out.append((start, prev))
            start = prev = v
    if start is not None:
        out.append((start, prev))
    return out


def clade_bands(g, posts, min_p=0.5):
    """Light band behind the rows of a clade when a source's best placement
    is an internal branch (otherwise the placement is just a dot on a short
    stem and the reader cannot tell which tips it implicates)."""
    ax, lay = g.ax, g.lay
    row_of = {t.name: t.y for t in lay.tips}
    below = tip_sets(lay.root)          # display node -> tips below it, in the DRAWN orientation
    for sp in posts:
        b = sp.top
        if not b or not b.mapped or b.p < min_p:
            continue
        sg = _seg_at(lay, b.name, float(np.median(b.t)))
        # The band must cover the tips BELOW the dot as drawn. lay.clade_tips holds the clade in the orientation of the
        # DATABASE tree, which for branches between the database root (Node1) and the display root is the complement of
        # what lies below the dot in the rerooted display (e.g. a Node44 placement shaded the 7 archaic tips instead of
        # the modern-human clade the dot sits above). Taking the tips below the display node of the segment is
        # correct in both orientations.
        # OLD: tips = [t for t in lay.clade_tips.get(b.name, ()) if t in row_of]
        tips = [t for t in below.get(sg[2], ()) if t in row_of]
        if len(tips) < 2 or len(tips) > 0.8 * len(lay.tips):
            continue
        x0 = sg[2].parent.x
        for r0, r1 in _runs([row_of[t] for t in tips]):
            ax.add_patch(Rectangle((x0, r0 - 0.5), g.lab_x + g.dx(g.label_w_mm + 2.0) - x0,
                                   r1 - r0 + 1.0, facecolor=sp.color, alpha=0.08,
                                   edgecolor="none", zorder=0))


def hilite_map(posts, min_p=0.05):
    out = {}
    for sp in posts:
        for b in sp.branches[:2]:
            if b.mapped and b.p >= min_p and b.name not in out:
                out[b.name] = sp.color
    return out


def draw_key(g, posts, y=None):
    """Marker-size key, right-aligned under the label column."""
    ax, lay = g.ax, g.lay
    n = len(lay.tips)
    y = (n - 0.5 + 1.3) if y is None else y
    fs = FS_TXT - 0.7
    head = "Branch posterior probability"
    entries = [(1.0, "1"), (0.5, "0.5"), (0.1, "0.1")]
    w_mm = text_w_mm(head, fs) + sum(3.4 + text_w_mm(l, fs) + 2.2 for _, l in entries)
    xr = ax.get_xlim()[1]
    x = xr - g.dx(w_mm + 0.5)
    ax.text(x, y, head, fontsize=fs, color=INK_2, ha="left", va="center")
    x += g.dx(text_w_mm(head, fs) + 3.0)
    for p, lab in entries:
        ax.scatter([x], [y], s=_marker_area(p), color=INK_3, edgecolors="white",
                   linewidths=0.45, zorder=8, clip_on=False)
        ax.text(x + g.dx(1.6), y, lab, fontsize=fs, color=INK_2, ha="left", va="center")
        x += g.dx(3.4 + text_w_mm(lab, fs) + 2.2 - 0.6)


# ----- legend row above a tree --------------------------------------------- #

def header_legend(fig, posts, x_mm, y_mm_from_top, max_w_mm):
    """One row of  [swatch] Source K  NAME  entries, wrapped if needed."""
    W, H = fig.get_size_inches() * 25.4
    ax = fig.add_axes([0, 0, 1, 1], facecolor="none")
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    x, y = x_mm, y_mm_from_top
    for sp in posts:
        b = sp.top
        lab1 = f"Source {sp.label_idx}"
        lab2 = b.name if b else ""
        w = 3.2 + text_w_mm(lab1, FS_TXT, "bold") + 1.6 + text_w_mm(lab2, FS_TXT) + 4.5
        if x + w > x_mm + max_w_mm and x > x_mm:
            x, y = x_mm, y + 3.6
        ax.add_patch(Rectangle((x, y - 1.2), 2.2, 2.4, facecolor=sp.color, edgecolor="none"))
        ax.text(x + 3.2, y, lab1, fontsize=FS_TXT, fontweight="bold", va="center", ha="left")
        ax.text(x + 3.2 + text_w_mm(lab1, FS_TXT, "bold") + 1.6, y, lab2, fontsize=FS_TXT,
                va="center", ha="left", color=INK_2)
        x += w
    return y


def footer(fig, meta, x_mm=2.0, y_mm=2.0):
    """Provenance line (sample, chain, iterations, database); wraps onto
    more lines when the figure is narrow instead of running off the edge."""
    if not meta:
        return
    W, H = fig.get_size_inches() * 25.4
    fs = FS_TXT - 1.5
    parts = meta.split("  \u00b7  ")
    lines, cur = [], ""
    for part in parts:
        trial = (cur + "  \u00b7  " + part) if cur else part
        if cur and text_w_mm(trial, fs) > W - 2 * x_mm:
            lines.append(cur)
            cur = part
        else:
            cur = trial
    lines.append(cur)
    for i, ln in enumerate(reversed(lines)):
        fig.text(x_mm / W, (y_mm + i * 2.6) / H, ln, fontsize=fs, color=INK_3, ha="left",
                 va="bottom")


# ----- full figure builders ------------------------------------------------ #

def fig_tree(lay, posts, path, which=None, meta="", pitch_mm=2.55, internal=False,
             min_prob=0.0, png=False, width_mm=COL15):
    sel = [sp for sp in posts if which is None or sp.idx in which]
    if not sel:
        return None
    n = len(lay.tips)
    extra = 4.0
    top_mm, bot_mm = 11.0, 7.0
    ax_h = (n + extra) * pitch_mm
    H = top_mm + ax_h + bot_mm
    fig = new_fig(width_mm, H)
    label_w = tip_label_width_mm(lay)
    ax = add_axes_mm(fig, 1.5, bot_mm, width_mm - 3.0, ax_h)
    g = setup_tree_axes(ax, lay, width_mm - 3.0, pitch_mm, label_w, extra)
    clade_bands(g, sel)
    draw_tree_skeleton(g, internal)
    draw_placements(g, sel, min_prob=min_prob, annotate=ANNOTATE)
    draw_tip_labels(g, hilite_map(sel))
    y_bar = draw_scale_bar(g)
    draw_key(g, sel, y_bar)
    header_legend(fig, sel, 3.0, 4.5, width_mm - 6.0)
    footer(fig, meta, 3.0, 2.0)
    return save_fig(fig, path, png)


# --------------------------------------------------------------------------- #
#  Proportions, placement bars, traces
# --------------------------------------------------------------------------- #

def _kde_reflected(vals, grid):
    """Gaussian KDE of samples on [0,1], reflected at both boundaries (no mass
    leaks out). The bandwidth is Scott's rule applied to the SAMPLES (sd * n^-1/5),
    so the curve narrows as the posterior narrows. (Computing it on a mirrored
    copy of the data, whose spread is ~0.6 whatever the posterior, pins the
    kernel at ~0.09 and hides the narrowing.) Binned on the plotting grid and
    smoothed with a reflecting Gaussian filter; integrates to 1."""
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 5:
        return None
    nb = len(grid)
    step = 1.0 / nb
    h = float(np.std(vals, ddof=1)) * len(vals) ** (-0.2)
    h = max(h, 1.5 * step)                       # never finer than the grid
    cnt, _ = np.histogram(np.clip(vals, 0.0, 1.0), bins=nb, range=(0.0, 1.0))
    sm = gaussian_filter1d(cnt.astype(float), sigma=h / step, mode="reflect")
    return sm / (len(vals) * step)


def draw_proportions(ax, posts):
    nb = 4000                                    # fine enough for sd ~0.002 posteriors
    grid = (np.arange(nb) + 0.5) / nb            # bin centers
    dens = {}
    for sp in posts:
        if sp.prop is None or len(sp.prop) < 5 or np.std(sp.prop) < 1e-9:
            continue
        d = _kde_reflected(sp.prop, grid)
        if d is not None:
            dens[sp.idx] = d
    if not dens:
        return False
    shown = [sp for sp in posts if sp.idx in dens]
    ymax = max(d.max() for d in dens.values())
    ax_w_mm = ax.get_window_extent().width / ax.figure.dpi * 25.4
    ax_h_mm = ax.get_window_extent().height / ax.figure.dpi * 25.4

    # vertical budget (mm): headroom for the legend rows, curves, CI bars
    row_mm, head_mm = 3.3, 3.3 * len(shown) + 1.5
    ci_mm = 1.9
    foot_mm = 1.5 + ci_mm * len(shown)
    curve_mm = max(ax_h_mm - head_mm - foot_mm, 8.0)
    per_mm = ymax / curve_mm                       # data units per mm
    ax.set_ylim(-foot_mm * per_mm, (curve_mm + head_mm) * per_mm)
    ax.set_xlim(0, 1)

    for r, sp in enumerate(shown):
        d, v = dens[sp.idx], sp.prop
        lo, hi = ci(v)
        med = float(np.median(v))
        y_ci = -(1.4 + ci_mm * r) * per_mm
        ax.fill_between(grid, d, color=sp.color, alpha=0.20, linewidth=0, zorder=2)
        ax.plot(grid, d, color=sp.color, lw=0.9, zorder=3)
        ax.plot([lo, hi], [y_ci, y_ci], color=sp.color, lw=2.0, solid_capstyle="round", zorder=3)
        ax.scatter([med], [y_ci], s=9, color="white", edgecolors=sp.color, linewidths=0.8, zorder=4)
        # legend row: swatch, name (ink), median (95% CI) -- above every curve
        yr = 1.0 - (0.9 + row_mm * (r + 0.5)) / ax_h_mm
        tr = ax.transAxes
        ax.add_patch(Rectangle((0.0, yr - 0.9 / ax_h_mm), 1.7 / ax_w_mm, 1.8 / ax_h_mm,
                               transform=tr, facecolor=sp.color, edgecolor="none", clip_on=False))
        ax.text(3.0 / ax_w_mm, yr, sp.name, transform=tr, fontsize=FS_TXT, fontweight="bold",
                ha="left", va="center", color=INK)
        ax.text((3.0 + text_w_mm(sp.name, FS_TXT, "bold") + 2.0) / ax_w_mm, yr,
                f"{med:.2f}  ({lo:.2f}\u2013{hi:.2f})", transform=tr, fontsize=FS_TXT - 0.5,
                ha="left", va="center", color=INK_2)
    ax.set_yticks([])
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax.set_xlabel("Mixture proportion")
    ax.set_ylabel("Posterior density")
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    return True


def fig_proportions(posts, path, meta="", png=False):
    if not any(sp.prop is not None and len(sp.prop) > 5 and np.std(sp.prop) > 1e-9 for sp in posts):
        return None
    H = 64.0
    fig = new_fig(COL1, H)
    ax = add_axes_mm(fig, 8.0, 12.0, COL1 - 12.0, 42.0)
    if not draw_proportions(ax, posts):
        plt.close(fig)
        return None
    fig.text(3.0 / COL1, 1 - 4.0 / H,
             "Thick bar: 95% credible interval  \u00b7  open dot: median",
             fontsize=FS_TXT - 0.7, color=INK_3, ha="left", va="center")
    footer(fig, meta, 3.0, 1.5)
    for sp in posts:
        if sp.prop is not None and len(sp.prop) > 1:
            lo, hi = ci(sp.prop)
            print(f"proportion_{sp.idx}: mean={sp.prop.mean():.6f}, 95% CI=[{lo:.6f}, {hi:.6f}], "
                  f"min={sp.prop.min():.6f}, max={sp.prop.max():.6f}, n={len(sp.prop)}")
    return save_fig(fig, path, png)


def draw_placement_bars(ax, posts, top=4, min_p=0.005, label_mm=20.0):
    rows = []     # (label, p, color, kind)
    for sp in posts:
        bs = [b for b in sp.branches if b.p >= min_p][:top]
        rows.append(("hdr", f"Source {sp.label_idx}", sp.color))
        for b in bs:
            rows.append(("bar", b.name, sp.color, b.p))
        if len(sp.branches) > top:
            rest = sum(b.p for b in sp.branches[top:])
            if rest >= 0.005:
                rows.append(("more", f"{len(sp.branches) - top} other branches", sp.color, rest))
        rows.append(("gap",))
    if rows and rows[-1][0] == "gap":
        rows.pop()
    ypos, labels, vals, cols, hdr_y = [], [], [], [], []
    y = 0.0
    ax.set_xlim(0, 1.0)
    for r in rows:
        if r[0] == "hdr":
            axw = ax.get_window_extent().width / ax.figure.dpi * 25.4
            x0 = -(label_mm + 1.0) / axw
            hdr_y.append((y, r[1], r[2], x0, axw))
            y += 0.85
        elif r[0] in ("bar", "more"):
            ax.barh(y, r[3], height=0.62, color=r[2] if r[0] == "bar" else INK_3,
                    alpha=1.0 if r[0] == "bar" else 0.45, linewidth=0, zorder=3)
            ax.text(r[3] + 0.015, y, f"{r[3]:.2f}", fontsize=FS_TXT - 0.7, ha="left",
                    va="center", color=INK_2, zorder=4)
            ypos.append(y)
            labels.append(r[1])
            y += 1.0
        else:
            y += 0.55
    ax.set_ylim(y - 0.4, -0.7)
    for yy, txt, col, x0, axw in hdr_y:
        tr = mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData)
        ax.add_patch(Rectangle((x0, yy - 0.24), 1.7 / axw, 0.48, transform=tr, facecolor=col,
                               edgecolor="none", clip_on=False))
        ax.text(x0 + 3.0 / axw, yy, txt, fontsize=FS_TXT, fontweight="bold", ha="left",
                va="center", color=INK, transform=tr, clip_on=False)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=FS_TXT)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_xlabel("Posterior probability of branch")
    ax.xaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    return y


def placement_rows(posts, top=4):
    # rows (in units of one bar) the placement panel needs
    n = 0
    for sp in posts:
        n += 0.85 + min(top, len(sp.branches)) + (1 if len(sp.branches) > top else 0) + 0.55
    return n


def draw_trace(ax, posts, mcmc):
    n_it = len(mcmc)
    stride = max(1, n_it // 1500)
    it = np.arange(n_it)[::stride]
    drawn = False
    for sp in posts:
        if sp.prop is not None and len(sp.prop) == n_it and np.std(sp.prop) > 1e-9:
            ax.plot(it, sp.prop[::stride], color=sp.color, lw=0.45, alpha=0.95, zorder=3)
            drawn = True
    if drawn:
        ax.set_ylim(0, 1)
        ax.set_ylabel("Mixture proportion")
    elif "Loglikelihood" in mcmc.columns:
        ll = pd.to_numeric(mcmc["Loglikelihood"], errors="coerce").to_numpy()
        ax.plot(it, ll[::stride], color=INK_2, lw=0.45)
        ax.set_ylabel("Log-likelihood")
        drawn = True
    ax.set_xlim(0, n_it)
    ax.set_xlabel("MCMC iteration (after burn-in)")
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(
        lambda v, _: f"{v / 1000:g}k" if v >= 1000 else f"{v:g}"))
    ax.yaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    style_axes(ax, left=True)
    return drawn


def fig_composite(lay, posts, mcmc, path, meta="", pitch_mm=2.55, internal=False,
                  min_prob=0.0, png=False, chain_data=None, best_idx=None):
    n = len(lay.tips)
    extra = 4.0
    top_mm, bot_mm = 11.0, 7.0
    ax_h = (n + extra) * pitch_mm
    W = COL2
    tree_w = 108.0
    x0 = 118.0
    wcol = W - x0 - 3.0
    XL, TOP, gap = 11.0, 6.0, 4.0              # below-axes room, above-axes room, block gap

    has_prop = any(sp.prop is not None and len(sp.prop) > 5 and np.std(sp.prop) > 1e-9
                   for sp in posts)
    show_chains = bool(chain_data) and len(chain_data) > 1 and \
        any(len(p_) > 1 for p_, _ in chain_data.values())

    def plan(compact):
        top_n = 3 if compact else 4
        row_mm = 2.7 if compact else 3.1
        blocks = []                              # [kind, height_mm, stretchable, min_height]
        if has_prop:
            blocks.append(["prop", TOP + 30.0 + XL, True, TOP + 24.0 + XL])
        blocks.append(["place", TOP + placement_rows(posts, top_n) * row_mm + XL, False, 0])
        blocks.append(["trace", TOP + 24.0 + XL, True, TOP + 16.0 + XL])
        if show_chains:
            blocks.append(["chains", TOP + (4.4 if compact else 5.2) * len(chain_data) + XL, False, 0])
        return blocks, top_n, row_mm

    H0 = max(top_mm + ax_h + bot_mm, 125.0)
    avail = lambda H: H - 1.0 - bot_mm - 2.0
    blocks, top_n, row_mm = plan(False)
    natural = lambda bl: sum(b[1] for b in bl) + gap * (len(bl) - 1)
    if natural(blocks) > avail(H0):
        blocks, top_n, row_mm = plan(True)
        for b in blocks:                         # squeeze stretchables to their minimum
            if b[2] and natural(blocks) > avail(H0):
                b[1] = max(b[3], b[1] - (natural(blocks) - avail(H0)))
    H = H0
    if natural(blocks) > avail(H0):               # still too tall: grow the page (cap 247 mm)
        H = min(natural(blocks) + 1.0 + bot_mm + 2.0, MAX_H)
    spare = avail(H) - natural(blocks)
    if spare > 0:
        stretch = [b for b in blocks if b[2]]
        add = min(spare / max(len(stretch), 1), 40.0)
        for b in stretch:
            b[1] += add

    fig = new_fig(W, H)
    ax_t = add_axes_mm(fig, 1.5, H - top_mm - ax_h, tree_w, ax_h)
    label_w = tip_label_width_mm(lay)
    g = setup_tree_axes(ax_t, lay, tree_w, pitch_mm, label_w, extra)
    clade_bands(g, posts)
    draw_tree_skeleton(g, internal)
    draw_placements(g, posts, min_prob=min_prob, annotate=ANNOTATE)
    draw_tip_labels(g, hilite_map(posts))
    y_bar = draw_scale_bar(g)
    draw_key(g, posts, y_bar)
    header_legend(fig, posts, 5.0, 4.5, tree_w - 6.0)
    panel_letter(fig, "a", 0.8, 1.0)

    letters = iter("bcdef")
    y = 1.0
    for kind, h_blk, _, _ in blocks:
        h_ax = h_blk - TOP - XL
        if kind == "prop":
            axb = add_axes_mm(fig, x0 + 9.0, H - (y + TOP + h_ax), wcol - 11.0, h_ax)
            draw_proportions(axb, posts)
            title, xl = "Posterior mixture proportions", x0 + 9.0
        elif kind == "place":
            label_mm = min(max([text_w_mm(b.name, FS_TXT) for sp in posts for b in sp.branches[:top_n]]
                               + [14.0]) + 6.0, 30.0)
            axc = add_axes_mm(fig, x0 + label_mm + 8.0, H - (y + TOP + h_ax),
                              wcol - label_mm - 16.0, h_ax)
            draw_placement_bars(axc, posts, top=top_n, label_mm=label_mm)
            title, xl = "Placement posterior", x0 + 9.0
        elif kind == "trace":
            axd = add_axes_mm(fig, x0 + 11.0, H - (y + TOP + h_ax), wcol - 13.0, h_ax)
            if not draw_trace(axd, posts, mcmc):
                axd.remove()
                continue
            title, xl = "Trace", x0 + 11.0
        else:
            axe = add_axes_mm(fig, x0 + 14.0, H - (y + TOP + h_ax), wcol - 18.0, h_ax)
            draw_chain_bars(axe, chain_data, best_idx)
            title, xl = "Agreement between chains", x0 + 9.0
        panel_letter(fig, next(letters), x0, y)
        fig.text(xl / W, 1 - (y + 2.6) / H, title, fontsize=FS_TXT, fontweight="bold",
                 ha="left", va="center")
        y += h_blk + gap
    footer(fig, meta, 3.0, 2.0)
    return save_fig(fig, path, png)


# --------------------------------------------------------------------------- #
#  Chain comparison
# --------------------------------------------------------------------------- #

def chain_rows(chain_data):
    """Per chain: sources ordered by mean proportion (largest first), so a
    color means the same thing in every row even when the sampler swapped
    the slot numbers between chains."""
    rows = {}
    for c, (posts, _) in chain_data.items():
        have = [sp for sp in posts if sp.prop is not None and len(sp.prop)]
        if len(have) > 1:
            have.sort(key=lambda s_: -float(np.mean(s_.prop)))
        items = []
        for rank, sp in enumerate(have, start=1):
            lo, hi = ci(sp.prop)
            tb = sp.top
            items.append(dict(rank=rank, mean=float(np.mean(sp.prop)), lo=lo, hi=hi,
                              name=tb.name if tb else "", p=tb.p if tb else 0.0))
        rows[c] = items
    return rows


def draw_chain_bars(ax, chain_data, best_idx, label_size=None):
    rows = chain_rows(chain_data)
    chains = sorted(rows)
    ax.set_xlim(0, 1)
    ax.set_ylim(len(chains) - 0.35, -0.65)
    for i, c in enumerate(chains):
        left = 0.0
        for it in rows[c]:
            col = src_color(it["rank"])
            ax.barh(i, it["mean"], left=left, height=0.62, color=col, linewidth=0, zorder=3)
            if it["mean"] > 0.07:
                ax.text(left + it["mean"] / 2, i, f"{it['mean']:.2f}", fontsize=FS_TXT - 0.7,
                        color="white", ha="center", va="center", zorder=5, fontweight="bold")
            left += it["mean"]
        best = c == best_idx
        ax.text(-0.025, i, f"chain {c}", fontsize=FS_TXT, ha="right", va="center",
                fontweight="bold" if best else "normal", color=INK,
                transform=ax.get_yaxis_transform(), clip_on=False)
    ax.set_yticks([])
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_xlabel("Posterior mean mixture proportion")
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    return rows


def fig_chains(chain_data, best_idx, path, meta="", png=False):
    """Do the chains agree? One row per chain: stacked posterior-mean
    proportions plus the best-supported branch of each source."""
    chains = sorted(chain_data)
    row = 6.5
    H = 12.0 + row * len(chains) + 17.0
    fig = new_fig(COL2, H)
    ax = add_axes_mm(fig, 22.0, 14.0, 78.0, row * len(chains))
    axt = add_axes_mm(fig, 108.0, 14.0, 72.0, row * len(chains))
    axt.set_xlim(0, 1)
    axt.set_ylim(len(chains) - 0.35, -0.65)
    axt.axis("off")
    rows = draw_chain_bars(ax, chain_data, best_idx)
    for i, c in enumerate(chains):
        its = rows[c]
        for j, it in enumerate(its):
            yy = i + ((-0.2 if j == 0 else 0.2) if len(its) > 1 else 0.0)
            col = src_color(it["rank"])
            axt.add_patch(Rectangle((0.0, yy - 0.09), 0.03, 0.18, facecolor=col, edgecolor="none"))
            axt.text(0.05, yy, it["name"] + (f"   {it['p']:.2f}" if it["p"] < 0.995 else ""),
                     fontsize=FS_TXT - 0.3, ha="left", va="center", color=INK)
    fig.text(22.0 / COL2, 1 - 5.0 / H, "Mixture proportion (sources ordered largest first)",
             fontsize=FS_TXT, color=INK, ha="left", va="center", fontweight="bold")
    fig.text(108.0 / COL2, 1 - 5.0 / H, "Best-supported branch per source (posterior probability)",
             fontsize=FS_TXT, color=INK, ha="left", va="center", fontweight="bold")
    footer(fig, meta, 3.0, 2.0)
    return save_fig(fig, path, png)


# --------------------------------------------------------------------------- #
#  RPVG (initial, pre-MCMC) estimates
# --------------------------------------------------------------------------- #

_GBWT_NAME_RE = re.compile(r"^_gbwt_ref_(.+)_\d+_\d+$")


def strip_gbwt_wrapper(name):
    m = _GBWT_NAME_RE.match(name)
    return m.group(1) if m else name


def load_rpvg_haplotype_transcripts(path):
    df = pd.read_csv(path, sep="\t")
    df["Name_stripped"] = df["Name"].astype(str).map(strip_gbwt_wrapper)
    for col in ("ReadCount", "TPM", "HaplotypeProbability"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_rpvg_haplotype_combos(path):
    df = pd.read_csv(path, sep="\t")
    name_cols = sorted((c for c in df.columns if re.match(r"^Name_\d+$", c)),
                       key=lambda c: int(c.split("_")[1]))
    if not name_cols:
        return pd.DataFrame()

    def lab(row):
        return " + ".join(strip_gbwt_wrapper(str(row[c])) for c in name_cols if pd.notna(row[c]))

    return pd.DataFrame({
        "Combo": df.apply(lab, axis=1),
        "HaplotypingProbability": pd.to_numeric(
            df.get("HaplotypingProbability", pd.Series(dtype=float)), errors="coerce"),
    })


def fig_rpvg_tree(lay, rdf, path, value_col="TPM", meta="", pitch_mm=2.55, internal=False,
                  png=False):
    sub = rdf[pd.to_numeric(rdf[value_col], errors="coerce").fillna(0) > 0]
    vals = {r["Name_stripped"]: float(r[value_col]) for _, r in sub.iterrows()
            if r["Name_stripped"] in lay.by_name and lay.by_name[r["Name_stripped"]].is_tip}
    if not vals:
        print(f"  None of the nonzero RPVG {value_col} haplotypes matched a tree tip; skipping.",
              file=sys.stderr)
        return None
    vmax = max(vals.values())
    n = len(lay.tips)
    extra = 4.0
    top_mm, bot_mm = 11.0, 7.0
    ax_h = (n + extra) * pitch_mm
    H = top_mm + ax_h + bot_mm
    W = COL15 + 28.0
    fig = new_fig(W, H)
    tree_w = COL15 - 3.0
    label_w = tip_label_width_mm(lay)
    ax = add_axes_mm(fig, 1.5, bot_mm, tree_w, ax_h)
    g = setup_tree_axes(ax, lay, tree_w, pitch_mm, label_w, extra)
    draw_tree_skeleton(g, internal)
    top_names = {k: SOURCE_FALLBACK for k, v in vals.items() if v >= 0.05 * vmax}
    draw_tip_labels(g, {k: seq_color(v, vmax) for k, v in vals.items() if k in top_names})
    draw_scale_bar(g)
    # aligned bar track: exactly the tip rows of the tree
    axb = add_axes_mm(fig, 1.5 + tree_w + 1.0, bot_mm + extra * pitch_mm,
                      W - tree_w - 5.0, n * pitch_mm)
    axb.set_ylim(n - 0.5, -0.5)
    axb.set_xlim(0, vmax * 1.08)
    for t in lay.tips:
        v = vals.get(t.name, 0.0)
        if v > 0:
            axb.barh(t.y, v, height=0.72, color=seq_color(v, vmax), linewidth=0, zorder=3)
    axb.set_yticks([])
    axb.xaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
    axb.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(
        lambda v, _: f"{v / 1000:g}k" if abs(v) >= 1000 else f"{v:g}"))
    axb.set_xlabel(f"RPVG initial {value_col}")
    axb.xaxis.grid(True, color=GRID, lw=0.4)
    axb.set_axisbelow(True)
    style_axes(axb)
    axb.spines["left"].set_visible(False)
    header_y = 4.5
    fig.text(3.0 / W, 1 - header_y / H, "RPVG initial estimate (before MCMC refinement)",
             fontsize=FS_TXT, fontweight="bold", color=INK, ha="left", va="center")
    footer(fig, meta, 3.0, 2.0)
    return save_fig(fig, path, png)


def fig_rpvg_bar(rdf, path, value_col="TPM", top_n=25, meta="", png=False):
    sub = rdf[pd.to_numeric(rdf[value_col], errors="coerce").fillna(0) > 0]
    if sub.empty:
        return None
    sub = sub.sort_values(value_col, ascending=False).head(top_n)
    labels = sub["Name_stripped"].tolist()[::-1]
    values = sub[value_col].tolist()[::-1]
    vmax = max(values)
    n = len(labels)
    H = 14.0 + n * 3.0 + 6.0
    lab_mm = max(text_w_mm(s, FS_TXT) for s in labels) + 3.0
    fig = new_fig(COL1, H)
    ax = add_axes_mm(fig, lab_mm + 3.0, 10.0, COL1 - lab_mm - 8.0, n * 3.0 + 1.0)
    ax.barh(range(n), values, height=0.68, color=[seq_color(v, vmax) for v in values],
            linewidth=0, zorder=3)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=FS_TXT - 0.5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(f"RPVG initial {value_col}")
    ax.xaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    ax.set_ylim(-0.7, n - 0.3)
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    fig.text(3.0 / COL1, 1 - 4.0 / H, f"RPVG initial per-haplotype {value_col} (top {n})",
             fontsize=FS_TXT, fontweight="bold", color=INK, ha="left", va="center")
    footer(fig, meta, 3.0, 1.5)
    return save_fig(fig, path, png)


def fig_rpvg_combos(cdf, path, top_n=15, meta="", png=False):
    sub = cdf.dropna(subset=["HaplotypingProbability"]) if not cdf.empty else cdf
    if sub.empty:
        return None
    sub = sub.sort_values("HaplotypingProbability", ascending=False).head(top_n)
    labels = sub["Combo"].tolist()[::-1]
    values = sub["HaplotypingProbability"].tolist()[::-1]
    vmax = max(values)
    n = len(labels)
    lab_mm = min(max(text_w_mm(s, FS_TXT) for s in labels) + 3.0, 80.0)
    H = 14.0 + n * 3.0 + 6.0
    fig = new_fig(COL15, H)
    ax = add_axes_mm(fig, lab_mm + 3.0, 10.0, COL15 - lab_mm - 8.0, n * 3.0 + 1.0)
    ax.barh(range(n), values, height=0.68, color=[seq_color(v, vmax) for v in values],
            linewidth=0, zorder=3)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=FS_TXT - 0.5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Haplotyping probability")
    ax.xaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    ax.set_ylim(-0.7, n - 0.3)
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    fig.text(3.0 / COL15, 1 - 4.0 / H, f"RPVG initial haplotype combinations (top {n})",
             fontsize=FS_TXT, fontweight="bold", color=INK, ha="left", va="center")
    footer(fig, meta, 3.0, 1.5)
    return save_fig(fig, path, png)


def run_rpvg_plots(lay, rpvg_prefix, out_prefix, meta, opts):
    if not rpvg_prefix:
        return {}
    out = {}
    ht = f"{rpvg_prefix}RPVGHaplotypeTranscripts.txt"
    hp = f"{rpvg_prefix}RPVGHaplotypes.txt"
    print(f"\n--- RPVG initial estimates (prefix: {rpvg_prefix}) ---")
    if os.path.isfile(ht):
        rdf = load_rpvg_haplotype_transcripts(ht)
        for col in ("TPM", "ReadCount"):
            if col not in rdf.columns:
                continue
            p = fig_rpvg_tree(lay, rdf, f"{out_prefix}_rpvg_transcripts_tree_{col}.pdf", col,
                              meta, opts.pitch_mm, opts.intern, opts.png)
            if p:
                out[f"rpvg_tree_{col}"] = p
            p = fig_rpvg_bar(rdf, f"{out_prefix}_rpvg_transcripts_bar_{col}.pdf", col,
                             meta=meta, png=opts.png)
            if p:
                out[f"rpvg_bar_{col}"] = p
    else:
        print(f"  Not found: {ht}  (skipping per-haplotype RPVG plots)", file=sys.stderr)
    if os.path.isfile(hp):
        p = fig_rpvg_combos(load_rpvg_haplotype_combos(hp), f"{out_prefix}_rpvg_combos_bar.pdf",
                            meta=meta, png=opts.png)
        if p:
            out["rpvg_combos"] = p
    else:
        print(f"  Not found: {hp}  (skipping haplotype-combination RPVG plot)", file=sys.stderr)
    return out


# --------------------------------------------------------------------------- #
#  Run discovery, chain selection
# --------------------------------------------------------------------------- #

def load_summary(prefix):
    path = f"{prefix}summary.txt"
    if not os.path.isfile(path):
        return None, path
    info = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            k, _, v = line.partition(":")
            info[k.strip()] = v.strip()
    return info, path


def discover_run_files(run_prefix):
    branch_path = f"{run_prefix}BranchEstimate.txt"
    diag_path = f"{run_prefix}Diagnostics.txt"
    files = glob.glob(f"{run_prefix}Result*_chain*.mcmc.gz") or glob.glob(f"{run_prefix}Result*_chain*.mcmc")

    def idx(p):
        m = re.search(r"chain(\d+)\.mcmc", p)
        return int(m.group(1)) if m else 999

    return branch_path, diag_path, sorted(files, key=idx)


def _parse_branch_estimate(path):
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            f = line.split("\t")
            if f[0].strip() == "Source":
                header = [x.strip().replace(" ", "_") for x in f]
                continue
            if header is None:
                continue
            if len(f) < len(header):
                f += [""] * (len(header) - len(f))
            rows.append(dict(zip(header, [x.strip() for x in f])))
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    for c in df.columns:
        if c != "Source":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _ess(x):
    """Effective sample size of a trace (FFT autocorrelation, Geyer's initial positive sequence).
    NaN for a constant trace."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 10 or np.std(x) < 1e-12:
        return float("nan")
    x = x - x.mean()
    f = np.fft.rfft(x, 2 * n)
    acov = np.fft.irfft(f * np.conj(f))[:n] / n
    rho = acov / acov[0]
    tau = 1.0
    for t in range(1, n - 1, 2):
        pair = rho[t] + rho[t + 1]
        if pair < 0:
            break
        tau += 2.0 * pair
    return n / max(tau, 1.0)


def proportion_ess(mcmc):
    """min over the sources of the ESS of the mixture proportion trace (None if there is none)."""
    vals = []
    for c in mcmc.columns:
        if re.fullmatch(r"proportion_\d+", str(c)):
            e = _ess(pd.to_numeric(mcmc[c], errors="coerce").to_numpy())
            if np.isfinite(e):
                vals.append(e)
    return min(vals) if vals else None


def select_best_chain(branch_path, chain_files, theta_ess=None):
    """Best chain = highest min(Pi ESS, Placement ESS, proportion ESS) over its sources
    (a chain is only as good as its worst-mixing parameter); ties broken by the lowest
    mean autocorrelation. theta_ess = {chain index: min ESS of the proportion traces}.
    The proportions used to be ignored, so a chain whose proportions mixed badly (or
    whose sources swapped labels part-way) could be chosen."""
    def default():
        print("  Defaulting to chain 0.", file=sys.stderr)
        idxs = [int(re.search(r"chain(\d+)\.mcmc", f).group(1)) for f in chain_files
                if re.search(r"chain(\d+)\.mcmc", f)]
        return 0, (idxs[0] if idxs else 0)

    if not os.path.isfile(branch_path):
        print(f"Warning: BranchEstimate file not found: {branch_path}", file=sys.stderr)
        return default()
    df = _parse_branch_estimate(branch_path)
    if df.empty:
        print("Warning: BranchEstimate file parsed to empty table.", file=sys.stderr)
        return default()
    ess = next((c for c in df.columns if "Pi_ESS" in c and "Placement" not in c), None)
    plc = next((c for c in df.columns if "Placement_ESS" in c), None)
    acc = next((c for c in df.columns if "Autocorrelation" in c), None)
    chc = next((c for c in df.columns if c.lower() == "chain"), None)
    if ess is None or chc is None:
        print(f"Warning: could not find Pi_ESS or Chain column in {branch_path}.", file=sys.stderr)
        return default()
    agg = {}
    for ci_, grp in df.groupby(chc):
        a, b = grp[ess].min(), (grp[plc].min() if plc else grp[ess].min())
        th = (theta_ess or {}).get(int(ci_))
        agg[int(ci_)] = {"score": min(a, b) if th is None else min(a, b, th), "pi": a, "plc": b,
                         "th": th, "ac": grp[acc].mean() if acc else 1.0}
    print("\nChain selection summary (aggregated across sources):")
    print(f"  {'Chain':>6}  {'min(Pi_ESS)':>12}  {'min(Plc_ESS)':>13}  {'min(Theta_ESS)':>15}  {'mean_Autocorr':>14}  {'score':>10}")
    for c in sorted(agg):
        a = agg[c]
        th = "n/a" if a["th"] is None else f"{a['th']:.1f}"
        print(f"  {c:>6}  {a['pi']:>12.1f}  {a['plc']:>13.1f}  {th:>15}  {a['ac']:>14.4f}  {a['score']:>10.1f}")
    best = max(agg, key=lambda c: (agg[c]["score"], -agg[c]["ac"]))
    b = agg[best]
    th_s = "n/a" if b["th"] is None else f"{b['th']:.1f}"
    print(f"\nBest chain: chain{best}  (score={b['score']:.1f}, min_Pi_ESS={b['pi']:.1f}, "
          f"min_Plc_ESS={b['plc']:.1f}, min_Theta_ESS={th_s}, mean_Autocorr={b['ac']:.4f})")
    idxs = [int(re.search(r"chain(\d+)\.mcmc", f).group(1)) for f in chain_files
            if re.search(r"chain(\d+)\.mcmc", f)]
    if best in idxs:
        return idxs.index(best), best
    print(f"Warning: chain{best} not found in chain file list; defaulting to position 0.",
          file=sys.stderr)
    return 0, best


# --------------------------------------------------------------------------- #
#  Orchestration
# --------------------------------------------------------------------------- #

def meta_string(summary, chain_label, mcmc_rows, tree_note=""):
    bits = []
    if summary:
        if summary.get("sample"):
            bits.append(summary["sample"])
    if chain_label:
        bits.append(chain_label)
    if summary:
        if summary.get("k"):
            bits.append(f"k = {summary['k']}")
        if summary.get("iterations"):
            try:
                bits.append(f"{int(summary['iterations']):,} iterations")
            except ValueError:
                pass
        if DB_LABEL is None and summary.get("database_dir"):
            # The directory name is an internal build name; a trailing '_nohh' only says how the database was built
            # (without the removed sequence) and is not meant for publication, so it is not printed.
            nm = os.path.basename(summary["database_dir"].rstrip("/"))
            bits.append("database: " + re.sub(r"_nohh$", "", nm))
        elif DB_LABEL:
            bits.append("database: " + DB_LABEL)
    if mcmc_rows:
        bits.append(f"{mcmc_rows:,} posterior samples")
    if tree_note:
        bits.append(tree_note)
    return "  ·  ".join(bits)


def run_single_chain(lay, mcmc_path, out_prefix, chain_label, opts, summary=None, draw=True,
                     loaded=None, chain_data=None, best_idx=None):
    """Summarize one chain and (unless draw=False) write all of its figures."""
    if loaded is None:
        print(f"\n--- Processing {chain_label or 'chain'}  ({mcmc_path}) ---")
        mcmc = load_mcmc(mcmc_path)
        posts = summarize_chain(mcmc, lay, opts.sort_sources)
    else:
        posts, mcmc = loaded
        print(f"\n--- Drawing {chain_label or 'chain'} ---")
    meta = meta_string(summary, chain_label, len(mcmc))
    out = {}
    if not draw:
        return out, posts, mcmc

    p = fig_proportions(posts, f"{out_prefix}_proportions.pdf", meta, opts.png)
    if p:
        out["proportions"] = p
    if not posts or not any(sp.branches for sp in posts):
        print("  Warning: no branch position data found.", file=sys.stderr)
        return out, posts, mcmc
    for sp in posts:
        p = fig_tree(lay, posts, f"{out_prefix}_Source_{sp.label_idx}.pdf", which={sp.idx}, meta=meta,
                     pitch_mm=opts.pitch_mm, internal=opts.intern, min_prob=opts.min_prob,
                     png=opts.png)
        if p:
            out[f"source_{sp.label_idx}"] = p
    p = fig_tree(lay, posts, f"{out_prefix}_combined.pdf", meta=meta, pitch_mm=opts.pitch_mm,
                 internal=opts.intern, min_prob=opts.min_prob, png=opts.png)
    if p:
        out["combined"] = p
    p = fig_composite(lay, posts, mcmc, f"{out_prefix}_figure.pdf", meta=meta,
                      pitch_mm=opts.pitch_mm, internal=opts.intern, min_prob=opts.min_prob,
                      png=opts.png, chain_data=chain_data, best_idx=best_idx)
    if p:
        out["figure"] = p
    return out, posts, mcmc


def read_keep(arg):
    if not arg:
        return None
    if os.path.isfile(arg):
        with open(arg) as fh:
            return {ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")}
    return {x.strip() for x in arg.split(",") if x.strip()}


def apply_keep(root, edge_map, clade_tips, keep, placed_names):
    """Tips/clades to retain: the user's list plus every placed clade."""
    by_name = {n.name: n for n in root.preorder() if n.name}
    want = set()
    for nm in list(keep) + list(placed_names):
        if nm in clade_tips:
            want.update(clade_tips[nm])
        elif nm in by_name:
            want.update(t.name for t in by_name[nm].tips())
    if not want:
        print("Warning: --keep matched nothing; keeping the full tree.", file=sys.stderr)
        return root
    return prune_tree(root, want, edge_map)


def peek_placed(chain_files):
    """Names the MCMC actually placed sources on (so --keep never hides them)."""
    names = set()
    for f in chain_files:
        try:
            df = pd.read_csv(f, sep="\t", compression="gzip" if f.endswith(".gz") else None)
        except Exception:
            continue
        for i in (1, 2, 3):
            c = f"Source_{i}"
            if c in df.columns:
                vc = df[c].astype(str).str.strip().value_counts(normalize=True)
                names.update(vc[vc >= 0.01].index)
    return names


def main():
    ap = argparse.ArgumentParser(
        description="Publication-style figures for TrailMix: MCMC placement posterior on "
                    "the phylogeny, mixture proportions, chain comparison, RPVG initial estimate.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tree", nargs="?", default=None,
                    help="Newick tree. Given ALONE it is taken as a vgan trailmix output prefix "
                         "and everything is read from {prefix}summary.txt.")
    ap.add_argument("mcmc_file", nargs="?", default=None,
                    help="MCMC file (.tsv/.tsv.gz) -- single-chain mode")
    ap.add_argument("out_prefix", nargs="?", default=None,
                    help="Output prefix (single-chain mode); overridden by --out-prefix")
    ap.add_argument("--multi", action="store_true", help="auto-discover chains from --run-prefix")
    ap.add_argument("--run-prefix", metavar="PREFIX")
    ap.add_argument("--out-prefix", metavar="PREFIX", dest="out_prefix_flag")
    ap.add_argument("--intern", action="store_true", help="print internal node names")
    ap.add_argument("--orig-tree", metavar="TREE", help="same as --db (a tree file)")
    ap.add_argument("--db", metavar="DB_PATH",
                    help="TrailMix database dir (--tm-files) or its haps.treefile")
    ap.add_argument("--rpvg-prefix", metavar="PREFIX")
    ap.add_argument("--layout", choices=("phylogram", "unit", "cladogram"), default="phylogram")
    ap.add_argument("--break-long", default="auto", metavar="auto|off|LENGTH")
    ap.add_argument("--keep", metavar="LIST|FILE")
    ap.add_argument("--ladderize", action="store_true")
    ap.add_argument("--sort-sources", dest="sort_sources", action="store_true",
                    help="(default; kept for old command lines)")
    ap.add_argument("--slot-order", dest="sort_sources", action="store_false",
                    help="number sources by vgan's slot instead of by proportion")
    ap.set_defaults(sort_sources=True)
    ap.add_argument("--min-prob", type=float, default=0.0)
    ap.add_argument("--pitch-mm", type=float, default=2.55, help="row height per tip, mm")
    ap.add_argument("--png", action="store_true", help="also write 300 dpi PNGs")
    ap.add_argument("--db-label", default=None, metavar="TEXT",
                    help="database name printed in the figure footer (default: the database directory name without a "
                         "trailing _nohh; an empty string omits it)")
    ap.add_argument("--no-rerooted", action="store_true",
                    help="auto mode: draw the database tree itself instead of haps_rerooted.treefile")
    ap.add_argument("--annotate", action="store_true",
                    help="write branch probability / node name next to placement markers")
    ap.add_argument("--best-only", action="store_true",
                    help="only draw the best chain (as chainBest) plus the chain comparison; "
                         "skips the per-chain figures")
    opts = ap.parse_args()

    setup_style()
    global ANNOTATE, DB_LABEL
    ANNOTATE = opts.annotate
    DB_LABEL = opts.db_label

    summary = None
    auto = (opts.tree is not None and opts.mcmc_file is None and opts.out_prefix is None
            and not opts.multi and not opts.run_prefix)
    if auto:
        run_prefix = opts.tree
        summary, spath = load_summary(run_prefix)
        if summary is None:
            sys.exit(f"Error: no tree/MCMC file given, and no summary.txt found at {spath}.\n"
                     f"Pass a tree + mcmc file (or --multi --run-prefix), or rerun vgan trailmix "
                     f"with a build that writes summary.txt.")
        db_tree = summary.get("database_tree")
        if not db_tree or not os.path.isfile(db_tree):
            sys.exit(f"Error: {spath} has no usable 'database_tree' entry (got: {db_tree!r}).")
        print(f"Auto-configured from {spath}")
        for k in ("database_tree", "database_dir", "k", "chains"):
            print(f"  {k + ':':16s}{summary.get(k, '?')}")
        # The database tree (haps.treefile) is what the MCMC placements refer to;
        # the display tree (haps_rerooted.treefile, same directory) is the same
        # tree rooted for drawing. Draw the rerooted one when it exists and
        # recover the internal-node names rerooting drops from the database tree.
        shown = db_tree
        rr = os.path.join(os.path.dirname(db_tree), "haps_rerooted.treefile")
        if not opts.no_rerooted and os.path.isfile(rr):
            shown = rr
            print(f"  display_tree:   {rr}")
        if not opts.db and not opts.orig_tree:
            opts.db = db_tree
        opts.tree, opts.multi, opts.run_prefix = shown, True, run_prefix
        opts.out_prefix_flag = opts.out_prefix_flag or run_prefix
        opts.rpvg_prefix = opts.rpvg_prefix or run_prefix
    elif opts.tree is None:
        sys.exit("Error: provide a tree file, or a single output-prefix argument (see --help).")
    elif opts.run_prefix:
        summary, _ = load_summary(opts.run_prefix)

    root = read_tree(opts.tree)
    ref = opts.orig_tree
    if opts.db:
        ref = resolve_db_tree(opts.db) or ref
    edge_map = clade_tips = None
    if ref and os.path.abspath(ref) != os.path.abspath(opts.tree):
        db_root = read_tree(ref)
        res = map_db_edges(db_root, root)
        if res is None:
            print("Warning: display tree and database tree have different tips; "
                  "matching branches by name instead.", file=sys.stderr)
            recover_labels(root, db_root)
        else:
            edge_map, clade_tips = res
            print("Database branches mapped onto the display tree by tip split "
                  "(rooting-independent).")
    if edge_map is None:
        edge_map, clade_tips = build_edge_map(root)

    # Chain files (needed up front: --keep must never prune a placed clade)
    if opts.multi:
        if not opts.run_prefix:
            sys.exit("Error: --multi requires --run-prefix")
        branch_path, diag_path, chain_files = discover_run_files(opts.run_prefix)
        if not chain_files:
            sys.exit(f"Error: no chain files found matching {opts.run_prefix}Result*_chain*.mcmc.gz")
    else:
        if not opts.mcmc_file:
            sys.exit("Error: provide an MCMC file, or use --multi --run-prefix")
        chain_files = [opts.mcmc_file]

    keep = read_keep(opts.keep)
    if keep:
        root = apply_keep(root, edge_map, clade_tips, keep, peek_placed(chain_files))
    if opts.ladderize:
        ladderize(root)
    lay = compute_layout(root, edge_map, clade_tips, opts.layout, opts.break_long)
    print(f"Tree: {len(lay.tips)} tips, layout={opts.layout}"
          + (f", {lay.n_broken} long branch(es) truncated" if lay.n_broken else ""))

    if opts.multi:
        run_pfx = opts.run_prefix
        out_pfx = opts.out_prefix_flag or run_pfx
        rpvg_pfx = opts.rpvg_prefix or run_pfx
        print(f"Run prefix:       {run_pfx}\nBranchEstimate:   {branch_path}\n"
              f"Diagnostics:      {diag_path}\nChain files:      {len(chain_files)}\n"
              f"Output prefix:    {out_pfx}\nRPVG prefix:      {rpvg_pfx}")
        # pass 1: load + summarize every chain (the composite figure compares them)
        chain_data, idx_of, theta_ess = {}, {}, {}
        for pos, cf in enumerate(chain_files):
            m = re.search(r"chain(\d+)\.mcmc", cf)
            ci_ = int(m.group(1)) if m else pos
            idx_of[pos] = ci_
            print(f"Reading chain {ci_}: {cf}")
            mc = load_mcmc(cf)
            chain_data[ci_] = (summarize_chain(mc, lay, opts.sort_sources), mc)
            th = proportion_ess(mc)
            if th is not None:
                theta_ess[ci_] = th
        list_pos, best_idx = select_best_chain(branch_path, chain_files, theta_ess)
        best_actual = idx_of[list_pos]

        # pass 2: draw
        all_out = {}
        for pos, cf in enumerate(chain_files):
            ci_ = idx_of[pos]
            if opts.best_only and ci_ != best_actual:
                continue
            tag = "chainBest" if opts.best_only else f"chain{ci_}"
            lab = f"chain {ci_}" + (" (best)" if ci_ == best_actual else "")
            out, _, _ = run_single_chain(lay, cf, f"{out_pfx}_{tag}", lab, opts, summary,
                                         loaded=chain_data[ci_], chain_data=chain_data,
                                         best_idx=best_actual)
            all_out[pos] = (ci_, out)

        best_paths = all_out[list_pos][1] if list_pos in all_out else {}
        if not opts.best_only:
            print(f"\n--- Copying chain{best_actual} outputs \u2192 chainBest ---")
            for key, src in best_paths.items():
                dst = src.replace(f"_chain{best_actual}_", "_chainBest_")
                shutil.copy2(src, dst)
                if opts.png and os.path.isfile(re.sub(r"\.pdf$", ".png", src)):
                    shutil.copy2(re.sub(r"\.pdf$", ".png", src), re.sub(r"\.pdf$", ".png", dst))
                print(f"  {src}  \u2192  {dst}")

        if len(chain_data) > 1 and any(len(p) > 1 for p, _ in chain_data.values()):
            print("\n--- Chain comparison ---")
            fig_chains(chain_data, best_actual, f"{out_pfx}_chains.pdf",
                       meta_string(summary, "", 0), opts.png)

        run_rpvg_plots(lay, rpvg_pfx, out_pfx, meta_string(summary, "", 0), opts)
    else:
        out_pfx = opts.out_prefix_flag or opts.out_prefix or "tree_source"
        print(f"Tree file:     {opts.tree}\nMCMC file:     {opts.mcmc_file}\nOutput prefix: {out_pfx}")
        run_single_chain(lay, opts.mcmc_file, out_pfx, "", opts, summary)
        if opts.rpvg_prefix:
            run_rpvg_plots(lay, opts.rpvg_prefix, out_pfx, "", opts)

    print("\nDone.")


if __name__ == "__main__":
    main()
