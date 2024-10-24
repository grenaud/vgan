SHELL := /bin/bash

# Directories
SPIMAPDIR := ../dep/spimap/
VGDIR := ../dep/vg/
RPVGDIR := ../dep/rpvg/
BIN_DIR := ../bin/

# https://stackoverflow.com/questions/2527496/how-can-i-write-a-makefile-to-auto-detect-and-parallelize-the-build-with-gnu-mak
NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)

# Compiler and flags
CXX := g++ -w -g -std=c++2a -Wextra -Wshadow -Wpedantic

# Include directories
CLUSTER_FLAGS := -I${RPVGDIR}deps/gbwt/include/ \
                 -I${RPVGDIR}deps/sdsl-lite/include/ \
                 -I${RPVGDIR}deps/sparsepp/ \
                 -I${RPVGDIR}src/tests/ \
                 -I${RPVGDIR}deps/libvgio/deps/libhandlegraph/src/include/ \
                 -I${RPVGDIR}deps/eigen/ \
                 -I${RPVGDIR}deps/cxxopts/include/ \
                 -I${VGDIR}src/io \
                 -I${VGDIR}src/algorithms \
                 -I${VGDIR}deps/sdsl-lite/include/sdsl/ \
                 -I${VGDIR}src

CXXFLAGS := -fopenmp -fconcepts -Wno-unused-variable -Wno-char-subscripts \
            -Wno-reorder -Wno-unknown-pragmas -lm -O3 -lz -I. -I../lib/libgab/ \
            -I../lib/libgab/gzstream/ -I${VGDIR}include/ -I${VGDIR}deps/progress_bar/ \
            -I${VGDIR}src/io -I/usr/lib -I${VGDIR}src/algorithms \
            -I${VGDIR}deps/libvgio/include/ -I ${RPVGDIR}src/ $(CLUSTER_FLAGS) -c

# Linker flags
CLUSTER_LINKER_FLAGS := -L${RPVGDIR}lib/ \
                        -L${RPVGDIR}deps/sdsl-lite/lib/ \
                        -L${VGDIR}lib \
                        ${VGDIR}lib/libvg.a \
                        ${RPVGDIR}deps/gbwt/lib/libgbwt.a \
                        -ldl

# Prerequisites
VGAN_PREREQS := MCMC.o run_trailmix.o TrailMix.o soibean.o soibean_functions.o \
                compute_init_vec.o haplocart_functions.o map_giraffe_Euka.o damage.o \
                Euka.o gam2prof.o vgan.o Clade.o update_likelihood.o HaploCart.o \
                get_posterior.o process_mapping.o fa2fq.o get_p_obs_base.o \
                read_fasta.o map_giraffe.o rmdup.o load.o readPathHandleGraph.o \
                NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.o \
                ../lib/libgab/gzstream/gzstream.o \
                ${VGDIR}obj/subcommand/subcommand.o \
                ${VGDIR}obj/subcommand/giraffe_main.o \
                ${VGDIR}obj/subcommand/gamsort_main.o \
                ${VGDIR}obj/subcommand/gbwt_main.o \
                ${VGDIR}obj/config/allocator_config_jemalloc.o \
                ${VGDIR}obj/subcommand/filter_main.o \
                ${VGDIR}obj/subcommand/view_main.o

TEST_PREREQS := test.o $(VGAN_PREREQS)

FLAGS_STATIC := ${SPIMAPDIR}src/treevis.o ${SPIMAPDIR}src/newick.o \
                ${SPIMAPDIR}src/logging.o ${SPIMAPDIR}src/parsing.o \
                soibean.o soibean_functions.o haplocart_functions.o MCMC.o \
                compute_init_vec.o TrailMix.o run_trailmix.o gam2prof.o map_giraffe_Euka.o damage.o Euka.o \
                vgan.o Clade.o update_likelihood.o HaploCart.o  get_posterior.o \
                process_mapping.o fa2fq.o get_p_obs_base.o read_fasta.o \
                map_giraffe.o  rmdup.o load.o readPathHandleGraph.o NodeInfo.o \
                AlignmentInfo.o ../lib/libgab/libgab.a \
                ../lib/libgab/gzstream/gzstream.o \
                ${VGDIR}obj/subcommand/subcommand.o \
                ${VGDIR}obj/subcommand/giraffe_main.o \
                ${VGDIR}obj/subcommand/gamsort_main.o \
                ${VGDIR}obj/subcommand/filter_main.o \
                ${VGDIR}obj/subcommand/gbwt_main.o \
                ${VGDIR}obj/config/allocator_config_jemalloc.o \
                ${VGDIR}lib/libjemalloc.a \
                ${VGDIR}lib/libtabixpp.a \
                ${VGDIR}lib/libdw.a \
                ${VGDIR}lib/libdwfl.a \
                ${VGDIR}lib/libvg.a \
                ${VGDIR}lib/libxg.a \
                ${VGDIR}lib/libvgio.a \
                ${VGDIR}deps/elfutils/libdw/libdw.a \
                ${VGDIR}deps/elfutils/libdwfl/libdwfl.a \
                ${VGDIR}lib/libhts.a \
                ${VGDIR}lib/libdeflate.a \
                -Wl,-Bstatic \
                -lvg -lstructures -lprotobuf -ltabixpp -lvcflib -lgssw -lssw \
                -lboost_unit_test_framework -lncurses -lgcsa2 -lgbwtgraph -lgbwt \
                -ldivsufsort -ldivsufsort64 -lraptor2 -lpinchesandcacti \
                -l3edgeconnected -lsonlib -lfml -lbdsg -lxg -lsdsl -lhandlegraph \
                -lcairo -lgobject-2.0 -lffi -lglib-2.0 -lpcre -lpixman-1 \
                -lfontconfig -luuid -lexpat -lfreetype -lpng16 -lxcb-shm \
                -lxcb-render -lXext -lX11 -lxcb -lXau -lXdmcp -ljansson -latomic \
                -lsublinearLS \
                -Wl,-Bdynamic \
                -ldl -lpthread -lm -fopenmp \
                -lzstd -lbz2 -llzma -lz \
                -Wl,-rpath,${VGDIR}lib \
                -rdynamic \
                -ldwfl -ldwelf -lelf -lebl \
                -lboost_program_options \
                -lrpvg \
                -L${RPVGDIR}lib/ \
                -L${RPVGDIR}deps/sdsl-lite/lib/ \
                -L${RPVGDIR}src/ \
                -L${VGDIR}lib



FLAGS_DYNAMIC := ${SPIMAPDIR}src/treevis.o ${SPIMAPDIR}src/newick.o \
                 ${SPIMAPDIR}src/logging.o ${SPIMAPDIR}src/parsing.o \
                 haplocart_functions.o TrailMix.o run_trailmix.o  \
                 map_giraffe_Euka.o rmdup.o damage.o Euka.o vgan.o Clade.o \
                 update_likelihood.o HaploCart.o get_posterior.o \
                 process_mapping.o fa2fq.o MCMC.o compute_init_vec.o \
                 gam2prof.o soibean.o soibean_functions.o get_p_obs_base.o \
                 read_fasta.o map_giraffe.o load.o readPathHandleGraph.o \
                 NodeInfo.o AlignmentInfo.o ../lib/libgab/libgab.a \
                 ../lib/libgab/gzstream/gzstream.o \
                 ${VGDIR}obj/subcommand/subcommand.o \
                 ${VGDIR}obj/subcommand/giraffe_main.o \
                 ${VGDIR}obj/subcommand/gamsort_main.o \
                 ${VGDIR}obj/subcommand/gbwt_main.o \
                 ${VGDIR}obj/config/allocator_config_jemalloc.o \
                 ${VGDIR}obj/subcommand/filter_main.o \
                 ${VGDIR}lib/libjemalloc.a \
                 ${VGDIR}lib/libtabixpp.a \
                 ${VGDIR}lib/libdw.a \
                 ${VGDIR}lib/libdwfl.a \
                 ${VGDIR}lib/libvg.a \
                 ${VGDIR}lib/libxg.a \
                 ${VGDIR}lib/libvgio.a \
                 ${VGDIR}deps/elfutils/libdw/libdw.a \
                 ${VGDIR}deps/elfutils/libdwfl/libdwfl.a \
                 ${VGDIR}lib/libhts.a \
                 ${VGDIR}lib/libdeflate.a \
                 -fopenmp -Wl,-Bstatic -lvg -lstructures -ltabixpp -lvcflib \
                 -lgssw -lssw -lboost_unit_test_framework -lncurses -lgcsa2 \
                 -lgbwtgraph -lgbwt -ldivsufsort -ldivsufsort64 -lraptor2 \
                 -lpinchesandcacti -l3edgeconnected -lsonlib -lfml -lbdsg -lxg \
                 -lsdsl -lhandlegraph -lcairo -lgobject-2.0 -lffi -lglib-2.0 \
                 -lpcre -lpixman-1 -lfontconfig -luuid -lexpat -lfreetype \
                 -lpng16 -lxcb-shm -lxcb-render -lXext -lX11 -lxcb -lXau \
                 -lXdmcp -ljansson -latomic -lsublinearLS -Wl,-rpath,${VGDIR}lib \
                 -rdynamic -ldwfl -ldwelf -lelf -lebl -lboost_program_options \
                 -lzstd -lbz2 -llzma -lprotobuf -pthread -lz -Wl,-Bdynamic -lc \
                 -lpthread -lm -lrpvg -L${RPVGDIR}lib/ \
                 -L${RPVGDIR}deps/sdsl-lite/lib/ -L${VGDIR}lib \
                 ${VGDIR}lib/libvg.a ${RPVGDIR}deps/gbwt/lib/libgbwt.a -ldl

FLAGS_TEST := test.o $(filter-out vgan.o, ${FLAGS_DYNAMIC}) -lboost_unit_test_framework

# Targets
.PHONY: all clean autocomp static hcfiles eukafiles soibeanfiles

all: vgan hcfiles eukafiles soibeanfiles

autocomp:
	# source ../tools/autocomp.bash

hcfiles: hcfilesmade

hcfilesmade: vgmade
	wget -nc -l1 --recursive --no-parent -P ../share/vgan/ ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/
	mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/haplocart/hcfiles/* ../share/vgan/hcfiles
	rm -rf ../share/vgan/ftp.healthtech.dtu.dk
	touch hcfilesmade

eukafiles: eukafilesmade

eukafilesmade: vgmade
	wget -nc -l0 --recursive --no-parent -P ../share/vgan/ ftp://ftp.healthtech.dtu.dk/public/euka_files/
	mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/euka_files/* ../share/vgan/euka_dir
	rm -rf ../share/vgan/ftp.healthtech.dtu.dk
	touch eukafilesmade

soibeanfiles: soibeanfilesmade

soibeanfilesmade: vgmade
	rm -rf ../share/vgan/soibean_dir/tree_dir
	wget -nc -l0 --recursive --no-parent -P ../share/vgan/ ftp://ftp.healthtech.dtu.dk/public/soibean_files/
	mv -fv ../share/vgan/ftp.healthtech.dtu.dk/public/soibean_files/* ../share/vgan/soibean_dir
	rm -rf ../share/vgan/ftp.healthtech.dtu.dk
	touch soibeanfilesmade

# Build libgab if the stamp file doesn't exist
libgabmade:
	rm -rf ../lib/libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../lib/libgab/
	cd ../lib/libgab/ && ${MAKE} -j $(NPROCS) && ${MAKE} -j $(NPROCS) -C gzstream/
	touch libgabmade

vgmade:
	cd $(VGDIR) && ./source_me.sh
	cd $(VGDIR)/deps/sparsehash && ./configure
	$(MAKE) -C $(VGDIR) -j $(NPROCS) || $(MAKE) -C $(VGDIR) -j $(NPROCS)
	touch vgmade

# Build rpvg if the stamp file doesn't exist
rpvgmade: vgmade
	cp $(RPVGDIR)../rpvg_src/*.cpp $(RPVGDIR)src
	cp $(RPVGDIR)../rpvg_src/*.hpp $(RPVGDIR)src
	mkdir -p $(RPVGDIR)build
	(cd ../dep/rpvg/build && cmake .. && ${MAKE} -j $(NPROCS))
	touch rpvgmade

# Build spimap if the stamp file doesn't exist
spimapmade:
	cd $(SPIMAPDIR) && make -j $(NPROCS)
	touch spimapmade

static: libgabmade rpvgmade spimapmade $(VGAN_PREREQS) vgmade
	$(CXX) -o $(BIN_DIR)vgan $(FLAGS_STATIC)

vgan: vgmade rpvgmade libgabmade $(VGAN_PREREQS)
	$(MAKE) autocomp
	$(CXX) $(CLUSTER_FLAGS) -o $(BIN_DIR)vgan $(FLAGS_DYNAMIC) $(CLUSTER_LINKER_FLAGS)

test: $(TEST_PREREQS) vgmade libgabmade
	$(MAKE) autocomp
	$(CXX) -o $(BIN_DIR)test $(FLAGS_TEST) $(CLUSTER_LINKER_FLAGS)
	rm -f ../test/output_files/trailmix/*.txt
	rm -f ../test/output_files/trailmix/*.mcmc
	rm -f ../test/output_files/haplocart/*.txt
	rm -f ../test/output_files/soibean/*.txt
	rm -f ../test/output_files/soibean/*.mcmc
	rm -f ../test/output_files/euka/*.tsv
	rm -f ../test/output_files/euka/*.prof

clean:
	rm -f *.o
	rm -f $(BIN_DIR)*

