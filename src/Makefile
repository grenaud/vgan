SHELL := /bin/bash

# Directory where the binaries are placed
BIN_DIR := ../bin/
VGDIR := ../dep/vg/

VGAN_PREREQS =  compute_init_vec.o MCMC.o gam2prof.o baseshift.o haplocart_functions.o map_giraffe_Euka.o damage.o Euka.o vgan.o Clade.o update_likelihood.o ../lib/libgab/libgab.o HaploCart.o get_posterior.o process_mapping.o fa2fq.o get_p_obs_base.o \
        read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.a \
        ../lib/libgab/gzstream/gzstream.o ${VGDIR}obj/subcommand/subcommand.o ${VGDIR}obj/subcommand/giraffe_main.o  \
        ${VGDIR}obj/subcommand/gamsort_main.o ${VGDIR}obj/subcommand/gbwt_main.o ${VGDIR}obj/subcommand/subcommand.o \
        ${VGDIR}obj/config/allocator_config_jemalloc.o ${VGDIR}obj/subcommand/filter_main.o ${VGDIR}obj/subcommand/view_main.o

TEST_PREREQS =  test.o #${VGAN_PREREQS}

VGAN_PREREQS = compute_init_vec.o MCMC.o haplocart_functions.o map_giraffe_Euka.o damage.o Euka.o gam2prof.o vgan.o Clade.o update_likelihood.o ../lib/libgab/libgab.o HaploCart.o get_posterior.o process_mapping.o fa2fq.o get_p_obs_base.o \
	read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.a \
	../lib/libgab/gzstream/gzstream.o ${VGDIR}obj/subcommand/subcommand.o ${VGDIR}obj/subcommand/giraffe_main.o  \
	${VGDIR}obj/subcommand/gamsort_main.o ${VGDIR}obj/subcommand/gbwt_main.o ${VGDIR}obj/subcommand/subcommand.o \
	${VGDIR}obj/config/allocator_config_jemalloc.o ${VGDIR}obj/subcommand/filter_main.o ${VGDIR}obj/subcommand/view_main.o


# https://stackoverflow.com/questions/2527496/how-can-i-write-a-makefile-to-auto-detect-and-parallelize-the-build-with-gnu-mak
NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)

# Only take half as many processors as available
NPROCS := $(shell echo "$(NPROCS)/2"|bc)

CXX      = g++ -s -std=c++2a -w





FLAGS_STATIC =  haplocart_functions.o MCMC.o compute_init_vec.o gam2prof.o baseshift.o     map_giraffe_Euka.o damage.o Euka.o vgan.o Clade.o update_likelihood.o HaploCart.o get_posterior.o \
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

FLAGS_DYNAMIC = haplocart_functions.o \
                 gam2prof.o map_giraffe_Euka.o damage.o Euka.o vgan.o Clade.o update_likelihood.o HaploCart.o get_posterior.o \
                 process_mapping.o fa2fq.o MCMC.o compute_init_vec.o baseshift.o \
                 get_p_obs_base.o read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.a \
                ../lib/libgab/gzstream/gzstream.o ../dep/vg/obj/subcommand/subcommand.o ../dep/vg/obj/subcommand/giraffe_main.o \
                ../dep/vg/obj/subcommand/gamsort_main.o ../dep/vg/obj/subcommand/gbwt_main.o ../dep/vg/obj/config/allocator_config_jemalloc.o \
                ../dep/vg/obj/subcommand/filter_main.o ../dep/vg/lib/libjemalloc.a ../dep/vg//lib/libtabixpp.a ../dep/vg//lib/libdw.a ../dep/vg//lib/libdwfl.a \
                ../dep/vg/lib/libvg.a ../dep/vg/lib/libxg.a ../dep/vg/lib/libvgio.a ../dep/vg/deps/elfutils/libdw/libdw.a ../dep/vg/deps/elfutils/libdwfl/libdwfl.a \
                ../dep/vg//lib/libhts.a ../dep/vg//lib/libdeflate.a -ldl -L../dep/vg/lib -fopenmp -Wl,-Bstatic -lvg  -lstructures  -ltabixpp \
                -lvcflib -lgssw -lssw -lboost_unit_test_framework -lncurses -lgcsa2 -lgbwtgraph -lgbwt -ldivsufsort -ldivsufsort64 -lraptor2 -lpinchesandcacti \
                -l3edgeconnected -lsonlib -lfml -lbdsg -lxg -lsdsl -lhandlegraph -lcairo -lgobject-2.0 -lffi -lglib-2.0 -lpcre -lpixman-1 -lfontconfig -luuid -lexpat \
                -lfreetype -lpng16 -lxcb-shm -lxcb-render -lXext -lX11 -lxcb -lXau -lXdmcp -ljansson -latomic -lsublinearLS "-Wl,-rpath,"../dep/vg/"lib" -rdynamic \
                -ldwfl -ldwelf -lelf -lebl "-lboost_program_options""" -lzstd -lbz2 -llzma -lprotobuf -pthread -lz -Wl,-Bdynamic -lpthread -lm \

FLAGS_TEST =    test.o $(filter-out vgan.o, ${FLAGS_DYNAMIC}) -lboost_unit_test_framework


CXXFLAGS =  -fopenmp -fconcepts -Wno-unused-variable -Wno-char-subscripts -Wno-reorder -Wno-unknown-pragmas -lm -O3 -lz  -I. -I../lib/libgab/ \
	    -I../lib/libgab/gzstream/ -I${VGDIR}include/ -I${VGDIR}deps/libvgio/include/  -I${VGDIR}src/ -I${VGDIR}deps/progress_bar/ -I${VGDIR}src/io \
	    -I/usr/lib -I${VGDIR}src/algorithms -c



FLAGS_STATIC = compute_init_vec.o MCMC.o haplocart_functions.o map_giraffe_Euka.o damage.o Euka.o gam2prof.o vgan.o Clade.o update_likelihood.o HaploCart.o get_posterior.o \
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


FLAGS_DYNAMIC = -L../dep/boost_1_61_0 compute_init_vec.o MCMC.o haplocart_functions.o gam2prof.o map_giraffe_Euka.o damage.o Euka.o vgan.o Clade.o update_likelihood.o HaploCart.o get_posterior.o process_mapping.o fa2fq.o \
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

 all:	vgan test  #hcfiles eukafiles 

.PHONY:	all clean autocomp libgab static hcfiles eukafiles

autocomp:
	source ../tools/autocomp.bash

#.NOTPARALLEL: hcfiles eukafiles vgmade test

${VGDIR}:
	rm -rf ${VGDIR} && git clone --recurse-submodules -j ${NPROCS} https://github.com/vgteam/vg.git ${VGDIR}

${VGDIR}include/vg/vg.pb.h vgmade: ${VGDIR}
	(cd ${VGDIR} && git checkout 8e14f9e43219ff7389cb197c340a8a8e4613f87e && git submodule update --init --recursive)
	(cd ${VGDIR} &&. ./source_me.sh && ${MAKE} -j ${NPROCS} && touch vgmade)

hcfiles hcfilesmade: vgmade
	wget -nc -l1 --recursive --no-parent -P ../share/vgan/ ftp://ftp.healthtech.dtu.dk:/public/haplocart/hcfiles/ && \
        mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/haplocart/hcfiles/* ../share/vgan/hcfiles && rm -rf ../share/vgan/ftp.healthtech.dtu.dk && touch hcfilesmade

eukafiles eukafilesmade: vgmade
	wget -nc -l1 --recursive --no-parent -P ../share/vgan/ ftp://ftp.healthtech.dtu.dk:/public/euka_files/ && \
	mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/euka_files/* ../share/vgan/euka_dir && rm -rf ../share/vgan/ftp.healthtech.dtu.dk && touch eukafilesmade

../lib/libgab/libgab.a: libgabmade
	echo ""

libgab libgabmade: vgmade
	rm -rf ../lib/libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../lib/libgab/
	(cd ../lib/libgab/ && ${MAKE} -j ${NPROCS} && ${MAKE} -j ${NPROCS} -C gzstream/) && touch libgabmade

static:	version.h vgmade libgabmade ${VGAN_PREREQS}
	${CXX} -o ${BIN_DIR}vgan $(FLAGS_STATIC)

version_chmod: version.h
	chmod +x get_version.sh

version.h: version_chmod 
	./get_version.sh

vgan: version.h vgmade libgabmade ${VGAN_PREREQS}
	${MAKE} autocomp && ${CXX} -o ${BIN_DIR}vgan $(FLAGS_DYNAMIC)

test: ${TEST_PREREQS}
	${MAKE} autocomp && ${CXX} -o ${BIN_DIR}test $(FLAGS_TEST)

clean:
	rm -f *.o && rm -f ../bin/* && rm -f -r ../test/output_files/* &&  mkdir -p ../test/output_files/euka


