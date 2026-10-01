#!/usr/bin/env bash
###############################################################################
# migrate_refs.sh
#
# Migrate vgan reference databases (HaploCart / euka / soibean) from the OLD
# vg (~v1.44, ODGI-based) format to the NEW vg v1.75 format, so vgan built
# against vg_newest can load and map against them.
#
# Two things change between old and new vg data:
#   1. The in-memory graph .og is ODGI format, which new bdsg CANNOT read.
#      -> Convert it to HashGraph with the OLD vg (old vg can still read ODGI).
#   2. The giraffe distance index (.dist) is an old SnarlDistanceIndex version
#      (v1022; new vg wants v5) and the minimizer (.min) is stale.
#      -> Rebuild .dist and .min (and, for euka/soibean, the .gbz) with NEW vg.
#
# VERIFIED end-to-end for HaploCart on this project (rCRS -> H2a2a1).
# euka/soibean steps use the same principles + `vg autoindex`; they are
# UNVERIFIED here only because their DBs (~10GB, 6.6GB distance index) exceed
# the 7GB test machine. Run them on a box with plenty of RAM (>=32GB rec.).
#
# Usage:
#   ./migrate_refs.sh haplocart
#   ./migrate_refs.sh euka
#   ./migrate_refs.sh soibean
###############################################################################
set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG -- edit these paths for your machine
# ---------------------------------------------------------------------------
# OLD vg binary that can READ ODGI (~v1.44). On the dev box this was:
OLD_VG="${OLD_VG:-/home/garen18/repos/vgan/dep/vg/bin/vg}"
# NEW vg v1.75 binary (build it, or use the vgan vg_newest tree's bin/vg).
# If you only built libvg.a (not bin/vg), build bin/vg or the minimal driver.
NEW_VG="${NEW_VG:-$PWD/vg_newest/bin/vg}"
# vgan share dir root (where euka_dir/ hcfiles/ soibean_dir/ live)
SHARE="${SHARE:-$PWD/share/vgan}"
FTP="ftp://ftp.healthtech.dtu.dk/public"
THREADS="${THREADS:-$(nproc)}"

log(){ echo -e "\n>>> $*" >&2; }
need(){ command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
need wget; [ -x "$OLD_VG" ] || { echo "OLD_VG not executable: $OLD_VG" >&2; exit 1; }
[ -x "$NEW_VG" ] || { echo "NEW_VG not executable: $NEW_VG" >&2; exit 1; }

# ODGI .og  ->  HashGraph (overwrite in place, keep .odgi.bak)
og_to_hashgraph(){
  local og="$1"
  log "Converting $(basename "$og") ODGI -> HashGraph (old vg)"
  "$OLD_VG" convert -a "$og" > "$og.hashgraph"
  mv "$og" "$og.odgi.bak"
  mv "$og.hashgraph" "$og"
  # sanity: HashGraph magic is 0x284d4f38
  local m; m=$(xxd -l4 -p "$og")
  [ "$m" = "284d4f38" ] || { echo "WARNING: unexpected magic $m (expected 284d4f38)"; }
}

# ---------------------------------------------------------------------------
# HaploCart  (ships a combined graph.giraffe.gbz -> simplest case; VERIFIED)
# ---------------------------------------------------------------------------
migrate_haplocart(){
  local D="$SHARE/hcfiles"; mkdir -p "$D"; cd "$D"
  log "HaploCart: downloading files (skip the huge graph.xg, not needed)"
  for f in graph.giraffe.gbz graph.dist graph.og \
           k31_w11.min k17_w18.min path_supports parsed_pangenome_mapping \
           parents.txt children.txt graph_paths mappability.tsv; do
    wget -q -c "$FTP/haplocart/hcfiles/$f" -O "$f"
  done
  og_to_hashgraph "$D/graph.og"
  log "Rebuilding distance index (from the gbz)"
  "$NEW_VG" index -j graph.dist graph.giraffe.gbz
  log "Rebuilding minimizers"
  "$NEW_VG" minimizer -k 31 -w 11 -d graph.dist -o k31_w11.min graph.giraffe.gbz
  "$NEW_VG" minimizer -k 17 -w 18 -d graph.dist -o k17_w18.min graph.giraffe.gbz
  log "HaploCart done. First run auto-builds graph.shortread.withzip.min (slow, once)."
}

# ---------------------------------------------------------------------------
# euka / soibean  (ship .gg + .gbwt + .og separately, NOT a combined gbz)
#
# Cleanest, consistent path with new vg: rebuild the giraffe indexes from the
# shipped GFA via `vg autoindex`, which emits a consistent PREFIX.giraffe.gbz +
# PREFIX.dist + PREFIX.min. Then convert the .og to HashGraph for the
# classification step (Euka/soibean readPathHandleGraph loads DB.og).
#
# IMPORTANT: vgan's euka/soibean currently invoke giraffe with the multi-file
# style (-g DB.gg -H DB.gbwt -x DB.og -d DB.dist -m DB.min). Because autoindex
# rebuilds the graph, keep everything CONSISTENT by switching that invocation
# to the GBZ style (like HaploCart). See MIGRATION.txt "CODE CHANGE" section.
# ---------------------------------------------------------------------------
migrate_multifile(){
  local NAME="$1" DIR="$2" PREFIX="$3"   # e.g. euka euka_dir euka_db
  local D="$SHARE/$DIR"; mkdir -p "$D"; cd "$D"
  log "$NAME: downloading graph + haplotypes + metadata (skip old .dist/.min/.ry)"
  # .gfa = base graph (used to rebuild indexes); .og = graph for classification;
  # .gbwt/.gg kept for reference; clade/bins/etc are the metadata vgan needs.
  for f in ${PREFIX}.gfa ${PREFIX}.og ${PREFIX}.gbwt ${PREFIX}.gg \
           ${PREFIX}.clade ${PREFIX}.bins; do
    wget -q -c "$FTP/${NAME}_files/$f" -O "$f" || echo "  (optional/missing: $f)"
  done
  # euka needs an extra path-supports file:
  [ "$NAME" = euka ] && wget -q -c "$FTP/euka_files/euka_db_graph_path_supports" \
        -O euka_db_graph_path_supports || true
  # soibean needs baseFreq + tree_dir:
  if [ "$NAME" = soibean ]; then
    wget -q -c "$FTP/soibean_files/soibean_db.baseFreq" -O soibean_db.baseFreq || true
    wget -q -r -nH --cut-dirs=2 -np "$FTP/soibean_files/tree_dir/" -P . || true
  fi

  log "$NAME: rebuilding giraffe indexes from GFA (needs LOTS of RAM/time)"
  # Produces: ${PREFIX}.giraffe.gbz  ${PREFIX}.dist  ${PREFIX}.min
  "$NEW_VG" autoindex --workflow giraffe -g "${PREFIX}.gfa" -p "${PREFIX}" -t "$THREADS"

  # The minimizer HaploCart-style names (k31/k17) are not produced by autoindex;
  # if vgan expects ${PREFIX}.min it is already there. If it expects specific
  # k/w minimizers, build them explicitly from the new gbz instead:
  #   "$NEW_VG" minimizer -k 31 -w 11 -d ${PREFIX}.dist -o ${PREFIX}.min ${PREFIX}.giraffe.gbz

  og_to_hashgraph "$D/${PREFIX}.og"
  log "$NAME done. Point vgan's giraffe call at ${PREFIX}.giraffe.gbz (-Z). See MIGRATION.txt."
}

case "${1:-}" in
  haplocart) migrate_haplocart ;;
  euka)      migrate_multifile euka    euka_dir    euka_db ;;
  soibean)   migrate_multifile soibean soibean_dir soibean_db ;;
  *) echo "usage: $0 {haplocart|euka|soibean}"; exit 1 ;;
esac
log "ALL DONE for: $1"
