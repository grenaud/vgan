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

#vg vgmade:
#	mkdir ${VGDIR}
#	wget -nc -l0 --recursive --no-parent -P ${VGDIR} ftp://ftp.healthtech.dtu.dk:/public/vg/ && \
#	mv -fv ${VGDIR}ftp.healthtech.dtu.dk/public/vg/* ${VGDIR} && rm -rf ${VGDIR}ftp.healthtech.dtu.dk && touch vgmade

#${VGDIR}:
#	rm -rf ${VGDIR} && git clone --recurse-submodules -j ${NPROCS} https://github.com/vgteam/vg.git ${VGDIR}

#${VGDIR}include/vg/vg.pb.h vgmade: ${VGDIR}
#	(cd ${VGDIR} && git checkout de51eb534394bb93a406baa809ac079ed2223eeb && . ./source_me.sh && ${MAKE} -j ${NPROCS}) && touch vgmade

hcfiles hcfilesmade: vgmade
	wget -nc -l1 --recursive --no-parent -P ../share/vgan/ ftp://ftp.healthtech.dtu.dk:/public/haplocart/hcfiles/ && \
        mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/haplocart/hcfiles/* ../share/vgan/hcfiles && rm -rf ../share/vgan/ftp.healthtech.dtu.dk && touch hcfilesmade

eukafiles eukafilesmade: vgmade
	wget -nc -l0 --recursive --no-parent -P ../share/vgan/ ftp://ftp.healthtech.dtu.dk:/public/euka_files/ && \
	mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/euka_files/* ../share/vgan/euka_dir && rm -rf ../share/vgan/ftp.healthtech.dtu.dk && touch eukafilesmade

soibeanfiles soibeanfilesmade: vgmade
	rm -rf ../share/vgan/soibean_dir/tree_dir
	wget -nc -l0 --recursive --no-parent -P ../share/vgan/ ftp://ftp.healthtech.dtu.dk:/public/soibean_files/ && \
	mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/soibean_files/* ../share/vgan/soibean_dir && rm -rf ../share/vgan/ftp.healthtech.dtu.dk && touch soibeanfilesmade

../lib/libgab/libgab.a: libgabmade
	echo ""


libgabmade:

	# Removing and cloning spimap repository
	rm -rf ${SPIMAPDIR}
	git clone --recursive https://github.com/mdrasmus/spimap.git ../dep/spimap/
	(cd ${SPIMAPDIR} && make -j ${NPROCS}) && touch spimapmade

	# Setting up RPVG directory and downloading its contents
	mkdir -p ${RPVGDIR}
	wget -nc -l0 --recursive --no-parent -P ${RPVGDIR} ftp://ftp.healthtech.dtu.dk:/public/rpvg/
	mv -fv ${RPVGDIR}ftp.healthtech.dtu.dk/public/rpvg/* ${RPVGDIR}
	rm -rf ${RPVGDIR}ftp.healthtech.dtu.dk && touch rpvgmade

	# Removing and cloning libgab repository
	rm -rf ../lib/libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../lib/libgab/
	(cd ../lib/libgab/ && make -j ${NPROCS} && make -j ${NPROCS} -C gzstream/) && touch libgabmade

	mkdir ${VGDIR}
	wget -nc -l0 --recursive --no-parent -P ${VGDIR} ftp://ftp.healthtech.dtu.dk:/public/vg/ && \
	mv -fv ${VGDIR}ftp.healthtech.dtu.dk/public/vg/* ${VGDIR} && rm -rf ${VGDIR}ftp.healthtech.dtu.dk && touch vgmade

	${MAKE} vgan -j ${NPROCS}

static:	libgabmade rpvgmade spimapmade ${VGAN_PREREQS}
	${CXX} -o ${BIN_DIR}vgan $(FLAGS_STATIC)

vgan: libgabmade ${VGAN_PREREQS}
	${MAKE} autocomp && ${CXX} ${CLUSTER_FLAGS} -o ${BIN_DIR}vgan ${FLAGS_DYNAMIC} ${CLUSTER_LINKER_FLAGS}

test: ${TEST_PREREQS}
	${MAKE} autocomp && ${CXX} -o ${BIN_DIR}test ${FLAGS_TEST} ${CLUSTER_LINKER_FLAGS}

clean:
	rm -f *.o && rm -f ../dep/src/spimap/*.o && rm -f ../bin/*
