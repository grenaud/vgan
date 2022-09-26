SHELL := /bin/bash

# Directory where the binaries are placed
BIN_DIR := ../bin/

VGAN_PREREQS = haplocart_functions.o vgan.o update_likelihood.o ../lib/libgab/libgab.o HaploCart.o get_posterior.o process_mapping.o fa2fq.o get_p_obs_base.o \
        read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.a \
        ../lib/libgab/gzstream/gzstream.o ${VGDIR}obj/subcommand/subcommand.o ${VGDIR}obj/subcommand/giraffe_main.o  \
        ${VGDIR}obj/subcommand/gamsort_main.o ${VGDIR}obj/subcommand/gbwt_main.o ${VGDIR}obj/subcommand/subcommand.o \
        ${VGDIR}obj/config/allocator_config_jemalloc.o ${VGDIR}obj/subcommand/filter_main.o ${VGDIR}obj/subcommand/view_main.o

# https://stackoverflow.com/questions/2527496/how-can-i-write-a-makefile-to-auto-detect-and-parallelize-the-build-with-gnu-mak
NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)

# Only take half as many processors as available
NPROCS := $(shell echo "$(NPROCS)/2"|bc)

CXX      = g++ -s -std=c++2a -w
CXXFLAGS =  -fopenmp -fconcepts -Wno-unused-variable -Wno-char-subscripts -Wno-reorder -Wno-unknown-pragmas -lm -O3 -lz  -I. -I../lib/libgab/ \
            -I../lib/libgab/gzstream/ -I${VGDIR}include/ -I${VGDIR}deps/libvgio/include/  -I${VGDIR}src/ -I${VGDIR}deps/progress_bar/ -I${VGDIR}src/io \
            -I/usr/lib -I${VGDIR}src/algorithms -c

VGDIR := ../dep/vg/

FLAGS_STATIC = haplocart_functions.o vgan.o update_likelihood.o HaploCart.o get_posterior.o \
               process_mapping.o fa2fq.o get_p_obs_base.o read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o \
               AlignmentInfo.o ../lib/libgab/libgab.a ../lib/libgab/gzstream/gzstream.o ../dep/vg/obj/subcommand/subcommand.o \
               ../dep/vg/obj/subcommand/giraffe_main.o ../dep/vg/obj/subcommand/gamsort_main.o ../dep/vg/obj/subcommand/gbwt_main.o \
               ../dep/vg/obj/config/allocator_config_jemalloc.o ../dep/vg/obj/subcommand/filter_main.o ../dep/vg/lib/libjemalloc.a \
               ../dep/vg//lib/libtabixpp.a ../dep/vg//lib/libdw.a ../dep/vg//lib/libdwfl.a ../dep/vg/lib/libvg.a ../dep/vg/lib/libxg.a \
               ../dep/vg/lib/libvgio.a ../dep/vg/deps/elfutils/libdw/libdw.a ../dep/vg/deps/elfutils/libdwfl/libdwfl.a ../dep/vg//lib/libhts.a \
               ../dep/vg//lib/libdeflate.a -static -static-libstdc++ -static-libgcc -ldl -L../dep/boost_1_61_0 -L../dep/vg/lib -fopenmp -lvg -lstructures \
               -ltabixpp -lvcflib -lgssw -lssw -lboost_unit_test_framework -lncurses -lgcsa2 -lgbwtgraph -lgbwt -ldivsufsort -ldivsufsort64 -lraptor2 \
               -lpinchesandcacti -l3edgeconnected -lsonlib -lfml -lbdsg -lxg -lsdsl -lhandlegraph -lcairo -lgobject-2.0 -lffi -lglib-2.0 -lpcre -lpixman-1 \
               -lfontconfig -luuid -lexpat -lfreetype -lpng16 -lxcb-shm -lxcb-render -lXext -lX11 -lxcb -lXau -lXdmcp -ljansson -latomic -lsublinearLS \
               "-Wl,-rpath,"../dep/vg/"lib" -rdynamic -ldwfl -ldwelf -lelf -lebl "-lboost_program_options""" -ldl -lzstd -lbz2 -llzma -lprotobuf \
               -pthread -lz -lpthread -lm

FLAGS_DYNAMIC = -L../dep/boost_1_61_0 haplocart_functions.o vgan.o update_likelihood.o HaploCart.o get_posterior.o \
                process_mapping.o fa2fq.o \
                get_p_obs_base.o read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.a \
                ../lib/libgab/gzstream/gzstream.o ../dep/vg/obj/subcommand/subcommand.o ../dep/vg/obj/subcommand/giraffe_main.o \
                ../dep/vg/obj/subcommand/gamsort_main.o ../dep/vg/obj/subcommand/gbwt_main.o ../dep/vg/obj/config/allocator_config_jemalloc.o \
                ../dep/vg/obj/subcommand/filter_main.o ../dep/vg/lib/libjemalloc.a ../dep/vg//lib/libtabixpp.a ../dep/vg//lib/libdw.a ../dep/vg//lib/libdwfl.a \
                ../dep/vg/lib/libvg.a ../dep/vg/lib/libxg.a ../dep/vg/lib/libvgio.a ../dep/vg/deps/elfutils/libdw/libdw.a ../dep/vg/deps/elfutils/libdwfl/libdwfl.a \
                ../dep/vg//lib/libhts.a ../dep/vg//lib/libdeflate.a -ldl -L../dep/vg/lib -fopenmp -Wl,-Bstatic -lvg -lstructures -ltabixpp \
                -lvcflib -lgssw -lssw -lboost_unit_test_framework -lncurses -lgcsa2 -lgbwtgraph -lgbwt -ldivsufsort -ldivsufsort64 -lraptor2 -lpinchesandcacti \
                -l3edgeconnected -lsonlib -lfml -lbdsg -lxg -lsdsl -lhandlegraph -lcairo -lgobject-2.0 -lffi -lglib-2.0 -lpcre -lpixman-1 -lfontconfig -luuid -lexpat \
                -lfreetype -lpng16 -lxcb-shm -lxcb-render -lXext -lX11 -lxcb -lXau -lXdmcp -ljansson -latomic -lsublinearLS "-Wl,-rpath,"../dep/vg/"lib" -rdynamic \
                -ldwfl -ldwelf -lelf -lebl "-lboost_program_options""" -lzstd -lbz2 -llzma -lprotobuf -pthread -lz -Wl,-Bdynamic -lpthread -lm

 all:	vgan hcfiles

.PHONY:	all clean autocomp libgab static hcfiles

autocomp:
	source ../tools/autocomp.bash

.NOTPARALLEL: hcfiles libgab

${VGDIR}:
	rm -rf ${VGDIR} && git clone --recursive --depth 10 -j ${NPROCS} https://github.com/vgteam/vg.git ${VGDIR}

${VGDIR}include/vg/vg.pb.h vgmade: ${VGDIR}
	(cd ${VGDIR} && git checkout de51eb534394bb93a406baa809ac079ed2223eeb && . ./source_me.sh && ${MAKE} -j ${NPROCS}) && touch vgmade

hcfiles hcfilesmade: libgabmade vgmade
	wget -nc -l1 --recursive --no-parent -P ../share/ ftp://ftp.healthtech.dtu.dk:/public/haplocart/hcfiles/ && \
        mv -fv ../share/ftp.healthtech.dtu.dk/public/haplocart/hcfiles/* ../share/hcfiles && rm -rf ../share/ftp.healthtech.dtu.dk && touch hcfilesmade

../lib/libgab/libgab.a: libgabmade
	echo ""

libgab libgabmade: vgmade
	rm -rf ../lib/libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../lib/libgab/
	(cd ../lib/libgab/ && ${MAKE} -j ${NPROCS} && ${MAKE} -j ${NPROCS} -C gzstream/) && touch libgabmade

static:	libgabmade ${VGAN_PREREQS}
	${CXX} -o ${BIN_DIR}vgan $(FLAGS_STATIC)

vgan: libgabmade hcfilesmade ${VGAN_PREREQS}
	${MAKE} autocomp && ${CXX} -o ${BIN_DIR}vgan $(FLAGS_DYNAMIC)

clean:
	rm -f *.o && rm -f ../bin/*

