#!/usr/bin/env bash
# build_deps.sh -- fetch and build everything vgan links against, then build vgan. Run from the vgan checkout.
#
#   ./build_deps.sh            (JOBS=8 NICE="nice -n 19" ./build_deps.sh   to change the defaults)
#
# What it does (each step is skipped when its result already exists):
#   1 libgab   lib/libgab        clone + make; gzstream/version is renamed (it shadows the C++20 <version> header)
#   2 spimap   dep/spimap        clone + make
#   3 vg       vg_newest         clone, PINNED to one commit, build lib/libvg.a + the subcommand objects vgan links;
#                                dep/vg is a symlink to it. (Not bin/vg: it is not needed and needs a lot of RAM.)
#   4 rpvg     dep/rpvg          github.com/grenaud/rpvg, tag vgan-1 (upstream jonassibbesen/rpvg 46ba0bb + the changes made for
#                                vgan: rpvg_main() entry point, --vgan-temp-dir, ...; see VGAN.md there), compiled to lib/librpvg.a
#   5 SAFARI   dep/safari_vg/bin/vg   the released static binary TrailMix runs as a subprocess to map reads
#   6 vgan     bin/vgan          make -C src vgan
# Needs no sudo if the system packages (README) are installed; for a missing libdw-dev / libelf-dev / nasm /
# meson / ninja see "No sudo" in the README section.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOBS="${JOBS:-$(( $(nproc) / 2 > 0 ? $(nproc) / 2 : 1 ))}"
NICE="${NICE-nice -n 19}"
PINNED_VG_COMMIT="6193310e3"       # v1.75.0-258-g6193310e3; vgan, its link flags and the index formats were verified against this one
RPVG_URL="${RPVG_URL:-https://github.com/grenaud/rpvg.git}"
RPVG_REF="${RPVG_REF:-vgan-1}"     # tag of the vgan version of rpvg (based on upstream 46ba0bb, 2022-11-09)
SAFARI_URL="${SAFARI_URL:-https://github.com/grenaud/SAFARI/releases/download/v1.0.1/vg}"
VG_DIR="${VG_DIR:-$ROOT/vg_newest}"
cd "$ROOT"
say() { echo; echo ">>> $*"; }

say "checking tools"
missing=()
for t in g++ make git cmake pkg-config protoc bc autoconf automake libtool meson ninja curl patch; do
    command -v "$t" >/dev/null 2>&1 || missing+=("$t")
done
pkg-config --exists libdw libelf 2>/dev/null || missing+=("libdw-dev/libelf-dev")
if [ ${#missing[@]} -gt 0 ]; then
    echo "missing: ${missing[*]}"
    echo "Install the packages listed in the README (sudo apt-get install ...), or see 'No sudo' there."
    exit 1
fi

# 1 libgab
if [ ! -f lib/libgab/libgab.a ]; then
    say "libgab"
    [ -d lib/libgab/.git ] || { rm -rf lib/libgab; git clone --recursive https://github.com/grenaud/libgab.git lib/libgab; }
    $NICE make -C lib/libgab -j"$JOBS"
fi
[ -f lib/libgab/gzstream/version ] && mv lib/libgab/gzstream/version lib/libgab/gzstream/version.txt

# 2 spimap
if [ ! -f dep/spimap/src/newick.o ]; then
    say "spimap"
    mkdir -p dep
    [ -d dep/spimap/.git ] || { rm -rf dep/spimap; git clone --recursive https://github.com/mdrasmus/spimap.git dep/spimap; }
    $NICE make -C dep/spimap -j"$JOBS"
fi

# 3 vg (pinned)
if [ ! -f "$VG_DIR/lib/libvg.a" ]; then
    say "vg $PINNED_VG_COMMIT (the long step)"
    [ -d "$VG_DIR/.git" ] || git clone https://github.com/vgteam/vg.git "$VG_DIR"
    git -C "$VG_DIR" checkout "$PINNED_VG_COMMIT"
    # vg's cairo and pixman submodules use git://anongit.freedesktop.org, which is blocked by many firewalls (the clone then
    # hangs without an error). The gitlab mirrors have the same commits and use https.
    git -C "$VG_DIR" submodule init deps/cairo deps/pixman
    git -C "$VG_DIR" config submodule.deps/cairo.url https://gitlab.freedesktop.org/cairo/cairo.git
    git -C "$VG_DIR" config submodule.deps/pixman.url https://gitlab.freedesktop.org/pixman/pixman.git
    git -C "$VG_DIR" submodule update --init --recursive
    $NICE make -C "$VG_DIR" -j"$JOBS" lib/libvg.a
fi
$NICE make -C "$VG_DIR" -j"$JOBS" \
    obj/subcommand/subcommand.o obj/subcommand/giraffe_main.o obj/subcommand/gamsort_main.o obj/subcommand/gbwt_main.o \
    obj/subcommand/filter_main.o obj/subcommand/view_main.o obj/subcommand/options.o obj/config/allocator_config_jemalloc.o
mkdir -p dep && ln -sfn "$VG_DIR" dep/vg

# 4 rpvg
if [ ! -f dep/rpvg/lib/librpvg.a ]; then
    say "rpvg $RPVG_URL @ $RPVG_REF"
    if [ ! -d dep/rpvg/.git ]; then rm -rf dep/rpvg; git clone --recursive --branch "$RPVG_REF" "$RPVG_URL" dep/rpvg; fi
    git -C dep/rpvg checkout "$RPVG_REF"
    git -C dep/rpvg submodule update --init --recursive
    mkdir -p dep/rpvg/lib dep/rpvg/obj
    V="$VG_DIR"
    INC="-Isrc -Ideps/eigen -Ideps/sparsepp -Ideps/cxxopts/include -I$V/include -I$V/src -I$V/deps/libvgio/include \
 -I$V/deps/libhandlegraph/src/include -I$V/deps/sdsl-lite/include -I$V/deps/gbwt/include -I$V/deps/gbwtgraph/include \
 -I$V/deps/xg/src -I$V/deps/libbdsg/include"
    # the 10 core files; rpvg_main() itself is in the header src/rpvg_main.hpp, which vgan's run_trailmix.cpp compiles in
    ( cd dep/rpvg
      for f in fragment_length_dist paths_index alignment_path alignment_path_finder path_clusters read_path_probabilities \
               path_estimator path_posterior_estimator path_abundance_estimator threaded_output_writer; do
          $NICE g++ -std=c++17 -O3 -fopenmp -fPIC -w -DGIT_COMMIT=\"$(git rev-parse --short HEAD)\" $INC -c "src/$f.cpp" -o "obj/$f.o" &
      done
      wait
      ar rcs lib/librpvg.a obj/*.o )
fi

# 5 SAFARI
if [ ! -x dep/safari_vg/bin/vg ]; then
    say "SAFARI binary"
    mkdir -p dep/safari_vg/bin
    if [ -n "${SAFARI_BIN:-}" ]; then cp "$SAFARI_BIN" dep/safari_vg/bin/vg; else curl -L -o dep/safari_vg/bin/vg "$SAFARI_URL"; fi
    chmod +x dep/safari_vg/bin/vg
fi

# 6 vgan
say "vgan"
$NICE make -C src -j"$JOBS" vgan
say "done: $ROOT/bin/vgan ($(bin/vgan version 2>&1 | head -1)). Databases: see the README ('Database files')."
