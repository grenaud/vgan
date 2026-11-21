SHELL := /bin/bash
SPIMAPDIR := ../dep/spimap/
VGDIR := ../dep/vg/
RPVGDIR := ../dep/rpvg/

# Directory where the binaries are placed
BIN_DIR := ../bin/

VGAN_PREREQS = assembly.o soibean.o soibean_functions.o compute_init_vec.o MCMC.o haplocart_functions.o map_giraffe_Euka.o damage.o Euka.o gam2prof.o vgan.o Clade.o update_likelihood.o ../lib/libgab/libgab.o HaploCart.o get_posterior.o process_mapping.o fa2fq.o get_p_obs_base.o \
        read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o TrailMix.o run_trailmix.o ../lib/libgab/libgab.a \
        ../lib/libgab/gzstream/gzstream.o ${VGDIR}obj/subcommand/subcommand.o ${VGDIR}obj/subcommand/giraffe_main.o  \
        ${VGDIR}obj/subcommand/gamsort_main.o ${VGDIR}obj/subcommand/gbwt_main.o ${VGDIR}obj/subcommand/subcommand.o \
        ${VGDIR}obj/config/allocator_config_jemalloc.o ${VGDIR}obj/subcommand/filter_main.o ${VGDIR}obj/subcommand/view_main.o

TEST_PREREQS =  test.o ${VGAN_PREREQS}

# libvg.a and vg lib dir are now carried in FLAGS_DYNAMIC (see below). The old rpvg-bundled
# gbwt/sdsl libs are intentionally NOT linked: they are built against old vg and clash with v1.75.
CLUSTER_LINKER_FLAGS =

# https://stackoverflow.com/questions/2527496/how-can-i-write-a-makefile-to-auto-detect-and-parallelize-the-build-with-gnu-mak
NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)

# Only take half as many processors as available
NPROCS := $(shell echo "$(NPROCS)/2"|bc)

CXX      = g++ -w -std=c++2a -w -Wextra -Wshadow -Wpedantic -s ## -g


CLUSTER_FLAGS = -I ${RPVGDIR}deps/gbwt/include/ -I ${RPVGDIR}deps/sdsl-lite/include/ -I${RPVGDIR}deps/sparsepp/ \
                -I${RPVGDIR}src/tests/ -I${RPVGDIR}deps/libvgio/deps/libhandlegraph/src/include/ \
                -I${RPVGDIR}deps/eigen/ -I${RPVGDIR}deps/cxxopts/include/ -I${VGDIR}src/io -I${VGDIR}src/algorithms -I${VGDIR}deps/sdsl-lite/include/sdsl/ \
                -I${RPVGDIR}src/

CXXFLAGS =  -fopenmp -fconcepts -Wno-unused-variable -Wno-char-subscripts -Wno-reorder -Wno-unknown-pragmas -lm -O3 -lz  -I. -I../lib/libgab/ \
            -I../lib/libgab/gzstream/ -I${VGDIR}include/ -I${VGDIR}deps/progress_bar/ -I${VGDIR}src/io \
            -I/usr/lib -I${VGDIR}src/algorithms ${CLUSTER_FLAGS} -I${VGDIR}src -I${VGDIR}deps/libvgio/include/ -c

# Fully-static build mirrors vg v1.75's static target: same libs as FLAGS_DYNAMIC plus
# -static -static-libstdc++ -static-libgcc (no -Bstatic/-Bdynamic markers, no rpath/-rdynamic).
# All needed static archives are present (vg provides libcairo.a/libpixman-1.a; system provides
# the rest incl. libc.a, libbrotli*.a, libboost_program_options.a).
FLAGS_STATIC = ${FLAGS_DYNAMIC_OBJS} \
               -static -static-libstdc++ -static-libgcc -L../dep/vg/lib -fopenmp \
               -L../dep/rpvg/lib -lrpvg \
               ../dep/vg/lib/libvg.a \
               -lvcflib -lwfa2 -ltabixpp -lgssw -lssw -lsublinearLS -lncurses -lgcsa2 -lgbwtgraph -lgbwt -lkff \
               -ldivsufsort -ldivsufsort64 -lraptor2 -lpinchesandcacti -l3edgeconnected -lsonlib -lstructures -lbdsg \
               -lxg -lsdsl -lhandlegraph -lzstd -lcrypto -ldl -ldw -lelf -latomic "-lboost_program_options" -llzma -lbz2 \
               -lvgio -lhts -ldeflate -lz -lbz2 -llzma -lprotobuf -ljansson \
               ../dep/vg/lib/libjemalloc.a \
               -lcairo -lfontconfig -lbz2 -lbrotlicommon -lexpat -lfreetype -lpng16 -lbrotlidec -lpixman-1 \
               -pthread -lz -lpthread -lm

# Object files (vgan + spimap + libgab + vg subcommand/config objects)
FLAGS_DYNAMIC_OBJS = ${SPIMAPDIR}src/treevis.o ${SPIMAPDIR}src/newick.o ${SPIMAPDIR}src/logging.o ${SPIMAPDIR}src/parsing.o haplocart_functions.o \
                 map_giraffe_Euka.o rmdup.o damage.o Euka.o vgan.o Clade.o update_likelihood.o HaploCart.o get_posterior.o \
                 process_mapping.o fa2fq.o MCMC.o compute_init_vec.o gam2prof.o soibean.o soibean_functions.o assembly.o \
                 get_p_obs_base.o read_fasta.o map_giraffe.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o \
                 TrailMix.o run_trailmix.o ../lib/libgab/libgab.a \
                ../lib/libgab/gzstream/gzstream.o ../dep/vg/obj/subcommand/subcommand.o ../dep/vg/obj/subcommand/giraffe_main.o \
                ../dep/vg/obj/subcommand/gamsort_main.o ../dep/vg/obj/subcommand/gbwt_main.o ../dep/vg/obj/config/allocator_config_jemalloc.o \
                ../dep/vg/obj/subcommand/filter_main.o ../dep/vg/obj/subcommand/options.o

# Library flags mirror vg v1.75's own link recipe (make -pn: LD_LIB_FLAGS / LD_STATIC_LIB_FLAGS /
# LD_STATIC_LIB_DEPS / LD_EXE_LIB_FLAGS / LD_CAIRO_LIB_FLAGS). vg now uses system libdw (backward-cpp),
# builds headless cairo (no X11/xcb), and no longer ships libfml/elfutils.
FLAGS_DYNAMIC = ${FLAGS_DYNAMIC_OBJS} \
                -L../dep/vg/lib -fopenmp -Wl,-rpath,$(abspath ../dep/vg/lib) -rdynamic \
                -L../dep/rpvg/lib -lrpvg \
                ../dep/vg/lib/libvg.a \
                -lvcflib -lwfa2 -ltabixpp -lgssw -lssw -lsublinearLS -lncurses -lgcsa2 -lgbwtgraph -lgbwt -lkff \
                -ldivsufsort -ldivsufsort64 -lraptor2 -lpinchesandcacti -l3edgeconnected -lsonlib -lstructures -lbdsg \
                -lxg -lsdsl -lhandlegraph -lzstd -lcrypto -ldl -ldw -lelf -latomic "-lboost_program_options" -llzma -lbz2 \
                -Wl,-Bstatic -lvgio -lhts -ldeflate -lz -lbz2 -llzma -lprotobuf -ljansson -Wl,-Bdynamic \
                -lpthread -lm ../dep/vg/lib/libjemalloc.a \
                -lcairo -lfontconfig -lbz2 -lbrotlicommon -lexpat -lfreetype -lpng16 -lbrotlidec -lpixman-1 \
                -pthread -lz -lpthread -lm \

FLAGS_TEST =    test.o $(filter-out vgan.o, ${FLAGS_DYNAMIC}) -lboost_unit_test_framework


 all:	vgan hcfiles eukafiles soibeanfiles

.PHONY:	all clean autocomp libgab static hcfiles eukafiles soibeanfiles

autocomp:
	#source ../tools/autocomp.bash

#.NOTPARALLEL: hcfiles libgab

# The database files are .tar.gz archives on Google Drive; each unpacks into a folder of the same name in ../share/vgan/
GDRIVE = https://drive.usercontent.google.com/download

hcfiles hcfilesmade:
	wget -O ../share/vgan/hcfiles.tar.gz "${GDRIVE}?id=1gPue0IKdEz3UpY2KMbAHnXZXKvEcCdmQ&export=download&confirm=t" && \
	tar -xzf ../share/vgan/hcfiles.tar.gz -C ../share/vgan/ && rm -f ../share/vgan/hcfiles.tar.gz && touch hcfilesmade

eukafiles eukafilesmade:
	wget -O ../share/vgan/euka_dir.tar.gz "${GDRIVE}?id=1cHDZujMQaMiJxZHhSD56BqlziJWthhPk&export=download&confirm=t" && \
	tar -xzf ../share/vgan/euka_dir.tar.gz -C ../share/vgan/ && rm -f ../share/vgan/euka_dir.tar.gz && touch eukafilesmade

soibeanfiles soibeanfilesmade:
	rm -rf ../share/vgan/soibean_dir/tree_dir
	wget -O ../share/vgan/soibean_dir.tar.gz "${GDRIVE}?id=1Hzjs4Bc8V71OYp8vY7waI2Q2x4vi8fSP&export=download&confirm=t" && \
	tar -xzf ../share/vgan/soibean_dir.tar.gz -C ../share/vgan/ && rm -f ../share/vgan/soibean_dir.tar.gz && touch soibeanfilesmade

# Dependencies (libgab, spimap, rpvg, vg) are fetched and built by ../build_deps.sh. make never downloads them;
# if one is missing it stops here and says what to run.
DEP_FILES = ../lib/libgab/libgab.a ../dep/rpvg/lib/librpvg.a ../dep/vg/lib/libvg.a ../dep/spimap/src/newick.o

${DEP_FILES}:
	@echo "error: $@ is missing. Run ../build_deps.sh first (see ../README.md, 'Building vgan from source')." >&2; false

static:	${DEP_FILES} ${VGAN_PREREQS}
	${CXX} -o ${BIN_DIR}vgan $(FLAGS_STATIC)

vgan: ${DEP_FILES} ${VGAN_PREREQS}
	${MAKE} autocomp && ${CXX} ${CLUSTER_FLAGS} -o ${BIN_DIR}vgan ${FLAGS_DYNAMIC} ${CLUSTER_LINKER_FLAGS}

test: ${TEST_PREREQS}
	${MAKE} autocomp && ${CXX} -o ${BIN_DIR}test ${FLAGS_TEST} ${CLUSTER_LINKER_FLAGS}

clean:
	rm -f *.o && rm -f ../dep/src/spimap/*.o && rm -f ../bin/*
