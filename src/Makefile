SHELL := /bin/bash
SPIMAPDIR := ../dep/spimap/
VGDIR := ../dep/vg/
RPVGDIR := ../dep/rpvg/

# Directory where the binaries are placed
BIN_DIR := ../bin/

VGAN_PREREQS = soibean.o soibean_functions.o compute_init_vec.o MCMC.o haplocart_functions.o map_giraffe_Euka.o damage.o Euka.o gam2prof.o vgan.o Clade.o update_likelihood.o ../lib/libgab/libgab.o HaploCart.o get_posterior.o process_mapping.o fa2fq.o get_p_obs_base.o \
        read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.a \
        ../lib/libgab/gzstream/gzstream.o ${VGDIR}obj/subcommand/subcommand.o ${VGDIR}obj/subcommand/giraffe_main.o  \
        ${VGDIR}obj/subcommand/gamsort_main.o ${VGDIR}obj/subcommand/gbwt_main.o ${VGDIR}obj/subcommand/subcommand.o \
        ${VGDIR}obj/config/allocator_config_jemalloc.o ${VGDIR}obj/subcommand/filter_main.o ${VGDIR}obj/subcommand/view_main.o

TEST_PREREQS =  test.o ${VGAN_PREREQS}

CLUSTER_LINKER_FLAGS = -L${RPVGDIR}lib/ -L${RPVGDIR}deps/sdsl-lite/lib/ -L../dep/vg/lib ../dep/vg/lib/libvg.a \
                       ${RPVGDIR}deps/gbwt/lib/libgbwt.a

# https://stackoverflow.com/questions/2527496/how-can-i-write-a-makefile-to-auto-detect-and-parallelize-the-build-with-gnu-mak
NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)

# Only take half as many processors as available
NPROCS := $(shell echo "$(NPROCS)/2"|bc)

CXX      = g++ -w -g -std=c++2a -w -Wextra -Wshadow -Wpedantic


CLUSTER_FLAGS = -I ${RPVGDIR}deps/gbwt/include/ -I ${RPVGDIR}deps/sdsl-lite/include/ -I${RPVGDIR}deps/sparsepp/ \
                -I${RPVGDIR}src/tests/ -I${RPVGDIR}deps/libvgio/deps/libhandlegraph/src/include/ \
                -I${RPVGDIR}deps/eigen/ -I${RPVGDIR}deps/cxxopts/include/ -I${VGDIR}src/io -I${VGDIR}src/algorithms -I${VGDIR}deps/sdsl-lite/include/sdsl/ \
                -I${RPVGDIR}src/

CXXFLAGS =  -fopenmp -fconcepts -Wno-unused-variable -Wno-char-subscripts -Wno-reorder -Wno-unknown-pragmas -lm -O3 -lz  -I. -I../lib/libgab/ \
            -I../lib/libgab/gzstream/ -I${VGDIR}include/ -I${VGDIR}deps/progress_bar/ -I${VGDIR}src/io \
            -I/usr/lib -I${VGDIR}src/algorithms ${CLUSTER_FLAGS} -I${VGDIR}src -I${VGDIR}deps/libvgio/include/ -c

FLAGS_STATIC = ${SPIMAPDIR}src/treevis.o ${SPIMAPDIR}src/newick.o ${SPIMAPDIR}src/logging.o ${SPIMAPDIR}src/parsing.o soibean.o soibean_functions.o haplocart_functions.o MCMC.o compute_init_vec.o gam2prof.o map_giraffe_Euka.o damage.o Euka.o vgan.o Clade.o update_likelihood.o HaploCart.o get_posterior.o \
               process_mapping.o  fa2fq.o get_p_obs_base.o read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o \
               AlignmentInfo.o ../lib/libgab/libgab.a ../lib/libgab/gzstream/gzstream.o ../dep/vg/obj/subcommand/subcommand.o \
               ../dep/vg/obj/subcommand/giraffe_main.o ../dep/vg/obj/subcommand/gamsort_main.o ../dep/vg/obj/subcommand/gbwt_main.o \
               ../dep/vg/obj/config/allocator_config_jemalloc.o ../dep/vg/obj/subcommand/filter_main.o ../dep/vg/lib/libjemalloc.a \
               ../dep/vg//lib/libtabixpp.a ../dep/vg//lib/libdw.a ../dep/vg//lib/libdwfl.a ../dep/vg/lib/libvg.a ../dep/vg/lib/libxg.a \
               ../dep/vg/lib/libvgio.a ../dep/vg/deps/elfutils/libdw/libdw.a ../dep/vg/deps/elfutils/libdwfl/libdwfl.a ../dep/vg//lib/libhts.a \
               ../dep/vg//lib/libdeflate.a -static -static-libstdc++ -static-libgcc -ldl -L../dep/boost_1_61_0 -L../dep/vg/lib -fopenmp -lvg  -lstructures \
               -ltabixpp -lvcflib -lgssw -lssw -lboost_unit_test_framework -lncurses -lgcsa2 -lgbwtgraph -lgbwt -ldivsufsort -ldivsufsort64 -lraptor2 \
               -lpinchesandcacti -l3edgeconnected -lsonlib -lfml -lbdsg -lxg -lsdsl -lhandlegraph -lcairo -lgobject-2.0 -lffi -lglib-2.0 -lpcre -lpixman-1 \
               -lfontconfig -luuid -lexpat -lfreetype -lpng16 -lxcb-shm -lxcb-render -lXext -lX11 -lxcb -lXau -lXdmcp -ljansson -latomic -lsublinearLS \
               "-Wl,-rpath,"../dep/vg/"lib" -rdynamic -ldwfl -ldwelf -lelf -lebl "-lboost_program_options""" -ldl -lzstd -lbz2 -llzma -lprotobuf \
               -pthread -lz -lpthread -lm

FLAGS_DYNAMIC = ${SPIMAPDIR}src/treevis.o ${SPIMAPDIR}src/newick.o ${SPIMAPDIR}src/logging.o ${SPIMAPDIR}src/parsing.o haplocart_functions.o \
                 map_giraffe_Euka.o rmdup.o damage.o Euka.o vgan.o Clade.o update_likelihood.o HaploCart.o get_posterior.o \
                 process_mapping.o fa2fq.o MCMC.o compute_init_vec.o gam2prof.o soibean.o soibean_functions.o \
                 get_p_obs_base.o read_fasta.o map_giraffe.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.a \
                ../lib/libgab/gzstream/gzstream.o ../dep/vg/obj/subcommand/subcommand.o ../dep/vg/obj/subcommand/giraffe_main.o \
                ../dep/vg/obj/subcommand/gamsort_main.o ../dep/vg/obj/subcommand/gbwt_main.o ../dep/vg/obj/config/allocator_config_jemalloc.o \
                ../dep/vg/obj/subcommand/filter_main.o ../dep/vg/lib/libjemalloc.a ../dep/vg//lib/libtabixpp.a ../dep/vg//lib/libdw.a ../dep/vg//lib/libdwfl.a \
                ../dep/vg/lib/libvg.a ../dep/vg/lib/libxg.a ../dep/vg/lib/libvgio.a ../dep/vg/deps/elfutils/libdw/libdw.a ../dep/vg/deps/elfutils/libdwfl/libdwfl.a \
                ../dep/vg//lib/libhts.a ../dep/vg//lib/libdeflate.a -ldl -L../dep/vg/lib -fopenmp -Wl,-Bstatic -lvg -lstructures  -ltabixpp \
                -lvcflib -lgssw -lssw -lboost_unit_test_framework -lncurses -lgcsa2 -lgbwtgraph -lgbwt -ldivsufsort -ldivsufsort64 -lraptor2 -lpinchesandcacti \
                -l3edgeconnected -lsonlib -lfml -lbdsg -lxg -lsdsl -lhandlegraph -lcairo -lgobject-2.0 -lffi -lglib-2.0 -lpcre -lpixman-1 -lfontconfig -luuid -lexpat \
                -lfreetype -lpng16 -lxcb-shm -lxcb-render -lXext -lX11 -lxcb -lXau -lXdmcp -ljansson -latomic -lsublinearLS "-Wl,-rpath,"../dep/vg/"lib" -rdynamic \
                -ldwfl -ldwelf -lelf -lebl "-lboost_program_options""" -lzstd -lbz2 -llzma -lprotobuf -pthread -lz -Wl,-Bdynamic -lpthread -lm -lrpvg \

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
	mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/euka_files/* ../share/vgan/euka_dir && rm -rf ../share/ftp.healthtech.dtu.dk && touch eukafilesmade

soibeanfiles soibeanfilesmade: vgmade
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
