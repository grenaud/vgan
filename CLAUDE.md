# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What vgan is

vgan (v3.0.0) is a C++ suite of pangenomics tools that runs **on top of vg** (github.com/vgteam/vg). It is not a standalone program: it statically links vg's `libvg.a` plus several of vg's compiled subcommand objects, and it invokes vg's mapper (`giraffe`) and other subcommands **in-process** by calling `vg::subcommand::Subcommand::get(argc, argv)` and running them. This makes vgan tightly coupled to one specific vg version and its on-disk index formats.

Subcommands (dispatched in `src/vgan.cpp`): `haplocart` (human mtDNA haplogroup), `euka` (bilaterian abundance from ancient eDNA), `soibean` (species ID / source estimation), `trailmix` (ancient hominin mtDNA mixture deconvolution / phylogenetic placement), `tempeh` (CDX pangenome coordinate index build/inspect + GAM coverage, a subprocess wrapper — see below), `gam2prof` (deamination profile from a GAM), `duprm` (PCR-dedup a GAM), `keelime`, `version`.

## Build

The build is driven by `src/Makefile` (not a root Makefile). It expects four dependencies wired under `dep/` and `lib/`:
- `dep/vg` → the vg source tree, **already built** (needs `libvg.a`, the deps `.a` libs, `obj/subcommand/*.o`, `obj/config/allocator_config_jemalloc.o`, `lib/libjemalloc.a`). Currently a symlink to `./vg_newest` (vg v1.75.0).
- `dep/spimap`, `dep/rpvg`, `lib/libgab` → prebuilt (symlinked to a sibling checkout here).

The Makefile's `libgabmade`/`vgmade` targets that auto-download deps from `ftp.healthtech.dtu.dk` are **dead** (FTP recipes obsolete). The setup is bypassed via marker files: `src/{libgabmade,vgmade,rpvgmade,spimapmade,*filesmade}` are touched so `make` skips re-downloading. If these markers are missing, `make` will try (and fail) to re-fetch deps.

```bash
make -C src vgan     # dynamic build  -> bin/vgan
make -C src static   # fully static build -> bin/vgan (-static -static-libstdc++ -static-libgcc)
make -C src test     # boost unit-test binary -> bin/test
make -C src clean     # rm *.o and bin/*
```

Both `vgan` and `static` targets write to `bin/vgan` — build one at a time. Compiler is `g++ -std=c++2a -fopenmp`. `FLAGS_DYNAMIC`/`FLAGS_STATIC` in `src/Makefile` mirror vg v1.75's own link recipe (from `make -pn` in vg: `LD_LIB_FLAGS` etc.); when the vg version changes, re-derive these.

### Building the vg dependency (memory-sensitive)

vgan does its own final link, so you only need vg's `libvg.a` + specific objects, **not** the full `bin/vg`. Building `bin/vg` is the memory-heavy step that OOMs low-RAM machines — avoid it. `build_with_vg_newest.sh` builds just what vgan links, memory-safely (`-j2`). vg itself needs system deps via `sudo make -C vg_newest get-deps` (protobuf, boost, meson, etc.).

## Run

Subcommands need their reference database dir:
```bash
./bin/vgan haplocart --hc-files share/vgan/hcfiles/ -fq1 reads.fq.gz -t 4      # or -f consensus.fa
./bin/vgan euka      --euka_dir  share/vgan/euka_dir/  -fq1 reads.fq.gz -t 4
./bin/vgan soibean   --soibean_dir share/vgan/soibean_dir/ ...
```
Test inputs live in `test/input_files/` (e.g. `rCRS.fa` → expected haplogroup `H2a2a1`; `euka/*.fq.gz`; `soibean/*.fq.gz`).

## Architecture notes (the non-obvious parts)

- **Graph loading uses `bdsg::HashGraph`.** The reference graph files are named `*.og` but must be **HashGraph** format for vg v1.75 (ODGI was removed from bdsg). `readPathHandleGraph.cpp` / `readOG_Euka.h` `deserialize()` them by value. vgan only *reads* graphs (get_handle / get_sequence / node ids / paths); it never mutates or serializes them.

- **Mapping pipeline (per subcommand).** Reads → `map_giraffe*.cpp` builds a giraffe argv and runs it in-process (stdout redirected to a FIFO). HaploCart/soibean use a 3-fork FIFO chain `map_giraffe → filter → gamsort → readGAM`; the analysis code then reads the GAM back. HaploCart `-f` consensus input is chunked into ~165 bp reads by `fa2fq` before mapping. Both euka and HaploCart use the combined `-Z <prefix>.giraffe.gbz -d <prefix>.dist -m <prefix>.min` giraffe style (euka/soibean's `map_giraffe_Euka.cpp` was switched from the old multi-file `-g .gg -H .gbwt -x .og` style during the v1.75 migration, see `MIGRATION.txt`).

- **soibean always classifies against a per-taxon subgraph, not the full DB.** `soibean_db.*` is only ever chunked (`vg chunk -r <range> -x soibean_db.og`) into a small per-taxon graph (e.g. `Ursidae.og/.gbwt/.giraffe.gbz/.dist/.min`, built by `share/vgan/soibean_dir/make_graph_files.sh`); `--dbprefix <taxon>` is mandatory and soibean.cpp throws without it. Never build a full `soibean_db.giraffe.gbz`/`.dist`/`.min` expecting to classify against it directly — only `.og` and `.gbwt` are read directly (by the shared `Euka::readPathHandleGraph`), for chunking.

- **vg subcommands beyond giraffe (index/minimizer/autoindex/chunk/mod/convert) don't need a full `bin/vg`.** vgan already proves you only need `libvg.a` + specific `obj/subcommand/*.o` objects; a minimal driver that calls `vg::subcommand::Subcommand::get()` (see `MIGRATION.txt`'s "vgtool workaround") avoids needing to build vg's cairo/pixman-dependent full CLI, which can fail for reasons unrelated to vgan (meson/toolchain version mismatches).

- **vg is tightly version-coupled, and pinned to one exact commit on purpose.** vgan links vg internals, so both the code and the reference-data index formats must match one vg version — and vg is actively developed, so a plain `git clone` gets whatever's on the default branch *today*, not a reproducible version. `vg_newest/` and `dep/vg/` are therefore kept in **detached HEAD** at commit `6193310e3` (`v1.75.0-258-g6193310e3`), not tracking `origin/master` — a stray `git pull` in either will refuse rather than silently drifting forward. `build_with_vg_newest.sh` re-pins to this same commit automatically (cloning fresh if needed). Upgrading vg on purpose means: (1) bump `PINNED_VG_COMMIT` in `build_with_vg_newest.sh`, (2) re-sweep source for renamed vg APIs, (3) re-derive `src/Makefile` link flags, (4) **regenerate the reference databases** (`.og`→HashGraph, and rebuild `.dist`/`.min` — the giraffe distance index has a hard version check), (5) re-run the euka/soibean/trailmix/tempeh end-to-end tests before trusting it. See `MIGRATION.txt`'s "VG IS PINNED TO ONE EXACT COMMIT" section and `migrate_refs.sh` for the data migration recipe.

- **Giraffe index cache gotcha.** On first run giraffe auto-builds `<prefix>.shortread.withzip.min` + `.shortread.zipcodes` next to the graph (slow). An interrupted build leaves them corrupt → every later run reports "no reads mapped". Recover with `rm share/vgan/<dir>/*.shortread.*` then rerun uninterrupted.

- **Known upstream bug:** `bdsg .../snarl_distance_index.cpp: Assertion is_chain(parent_chain) failed` can fire (in **both** static and dynamic builds) on some longer reads; it intermittently breaks HaploCart's `-f` consensus mode. Short-read (`-fq1`) input is unaffected. This is a vg/bdsg issue, not vgan.

- **TrailMix maps reads via a SAFARI subprocess, not in-process.** Unlike every other subcommand, `Trailmix::map_giraffe()` (`src/trailmix_functions.h`) does not call `vg::subcommand::Subcommand::get()` in-process — it forks and `execv()`s a separate prebuilt binary, `dep/safari_vg/bin/vg` (a damage-aware giraffe fork, github.com/grenaud/SAFARI, built from vg v1.44), redirecting its stdout GAM stream into the usual FIFO. This is a deliberate process-boundary choice: SAFARI's changes touch vg's core mapper internals in ways that don't line up with v1.75's zipcode-based seed clustering, so porting it in-process was judged much higher-risk than shelling out. TrailMix also needs `dep/rpvg` (`-lrpvg`, linked before `libvg.a` so libvg.a resolves its symbols) for haplotype-aware path abundance, and its own reference-data directories `share/vgan/tmfiles/` (full) and `share/vgan/publication_tmfiles/` (small test set, prefix `pub.graph`) — these must stay in **separate** directories because each has its own matching `haps.treefile`; mixing them pairs the wrong tree with the wrong graph and crashes MCMC. See `MIGRATION.txt` section D) for the full port history and bugs fixed.

- **`tempeh` is not a vg-linked port at all — it's a pure subprocess wrapper around a separate CMake project.** `dep/cdx/` is a full clone of github.com/JolanBoucher/cdx (a CDX pangenome coordinate-index builder + GAM coverage tool), built entirely by its **own CMake system** (`cmake -S dep/cdx -B dep/cdx/build && cmake --build dep/cdx/build`), completely independent of `src/Makefile` and never linked into `bin/vgan`. The merged dispatcher executable's CMake target was renamed from upstream's `cdx` to `tempeh` in `dep/cdx/CMakeLists.txt` (output binary is `dep/cdx/build/tempeh`), and every user-facing CLI string (usage/help text) was patched to match, so `vgan tempeh --help` reads as one consistent tool rather than a program still calling itself "cdx" — this is a real, intentional divergence from upstream to keep in mind on any future `git pull` in `dep/cdx`. (The on-disk *index format* is still called "CDX" / `.cdx` — that's a separate, unrelated concept from the executable's name and was left alone.) `run_tempeh_subprocess()` in `src/vgan.cpp` just `fork()`+`execv()`s `dep/cdx/build/tempeh` with argv forwarded unmodified — it picks build/inspect/coverage mode itself from the input file's binary signature and has its own complete CLI11-based `--help`, so vgan does no argument handling of its own. No FIFO needed (unlike SAFARI): it reads/writes real files directly, so stdio is just inherited. It fetches its own separate copies of libbdsg/GBWT/GBWTGraph/sdsl-lite/libvgio/Abseil via CMake FetchContent (deliberately not shared with vg_newest's copies, since it's a different process); the one dependency reused from vg_newest is its already-built HTSlib 1.19.1, exposed to its CMake via an isolated `dep/cdx_pkgconfig_hints/` directory containing only a symlinked `htslib.pc` — **do not** point `PKG_CONFIG_PATH` at the whole `vg_newest/lib/pkgconfig/` directory, since its `cairo.pc` has hardcoded paths from a different machine and breaks the cairo lookup. See `MIGRATION.txt` section E) for the full build recipe and gotchas.
