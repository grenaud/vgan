#!/usr/bin/env bash
# Build vgan against vg_newest, PINNED to a specific vg commit.
#
# vg is actively developed and its default branch moves often - a plain
# `git clone https://github.com/vgteam/vg.git` gets whatever is current
# HEAD *today*, not the commit vgan's source patches, src/Makefile link
# flags, and reference-database formats (HashGraph, SnarlDistanceIndex v5,
# etc.) were actually verified against. To avoid re-discovering/re-fixing
# vg-API breakage every time upstream changes something, this script pins
# to one exact commit and leaves both vg_newest/ and dep/vg/ in DETACHED
# HEAD state there (not tracking origin/master) - so an accidental `git
# pull` in either directory fails loudly instead of silently drifting.
#
# To move to a newer vg on purpose, update PINNED_VG_COMMIT below only
# after re-verifying euka/soibean/trailmix end-to-end (see MIGRATION.txt).
PINNED_VG_COMMIT="6193310e3"   # v1.75.0-258-g6193310e3 "Spike"
VG_REPO_URL="https://github.com/vgteam/vg.git"
#
# vgan's dep/ and lib/ have already been wired by Claude:
#   dep/vg     -> ./vg_newest                     (the new vg)
#   dep/spimap -> /home/garen18/repos/vgan/dep/spimap   (prebuilt, reused)
#   dep/rpvg   -> /home/garen18/repos/vgan/dep/rpvg     (prebuilt, reused)
#   lib/libgab -> /home/garen18/repos/vgan/lib/libgab   (prebuilt, reused)
# and src/{libgabmade,vgmade,rpvgmade,spimapmade,*filesmade} markers are touched
# so `make` will NOT try to re-download/re-clone any dependency.
#
# Run this from the vgan repo root: /home/garen18/projects/vgan/vgan
set -euo pipefail
cd "$(dirname "$0")"
ROOT="$PWD"
VG="$ROOT/vg_newest"

# --- Step 0: clone + pin vg_newest, if it doesn't already exist --------------
if [ ! -d "$VG/.git" ]; then
    echo ">>> Cloning vg (full history needed to reach the pinned commit)"
    git clone "$VG_REPO_URL" "$VG"
fi
echo ">>> Pinning vg_newest to $PINNED_VG_COMMIT (detached HEAD)"
git -C "$VG" fetch --unshallow 2>/dev/null || git -C "$VG" fetch origin
git -C "$VG" checkout "$PINNED_VG_COMMIT"
git -C "$VG" submodule update --init --recursive

# --- Step 1: system dependencies (needs sudo; only you can do this) -----------
# This installs protobuf, boost, jansson, meson, cairo, etc. Run once:
#     sudo make -C vg_newest get-deps
# (Uncomment to run here if you have sudo.)
# sudo make -C "$VG" get-deps

# --- Step 2: build vg_newest, MEMORY-SAFELY -----------------------------------
# This box has 7.2 GiB RAM + 4 GiB swap. The full `make` that links bin/vg is
# what OOM-kills. vgan does its own final link, so we only need libvg.a, a few
# subcommand objects, the jemalloc allocator config object, and the deps libs.
# Use -j2 to stay under the memory ceiling.
cd "$VG"
echo ">>> Building vg deps + libvg.a (-j2, this is the long part)"
make -j2 lib/libvg.a
echo ">>> Building the subcommand + config objects vgan links against"
make -j2 \
    obj/subcommand/subcommand.o \
    obj/subcommand/giraffe_main.o \
    obj/subcommand/gamsort_main.o \
    obj/subcommand/gbwt_main.o \
    obj/subcommand/filter_main.o \
    obj/subcommand/view_main.o \
    obj/config/allocator_config_jemalloc.o

# --- Step 3: build vgan -------------------------------------------------------
cd "$ROOT/src"
echo ">>> Building vgan"
make -j2 vgan
echo ">>> Done. Binary at: $ROOT/bin/vgan"
