FROM ubuntu:20.04


RUN apt-get -qq update

RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata

RUN apt-get -qq install -y \
    sudo \
    pv \
    pigz \
    bsdmainutils \
    build-essential \
    make \
    git \
    zlib1g-dev \
    rs \
    gdb \
    time \
    cmake \
    pkg-config\
    libncurses-dev \
    libbz2-dev  \
    protobuf-compiler \
    libprotoc-dev \
    libprotobuf-dev \
    libjansson-dev \
    automake \
    gettext \
    autopoint \
    libtool \
    jq \
    bsdmainutils \
    bc \
    rs \
    parallel \
    npm \
    curl \
    unzip \
    redland-utils \
    librdf-dev \
    bison \
    flex \
    libzstd-dev \
    lzma-dev \
    liblzma-dev \
    liblz4-dev \
    libffi-dev \
    libcairo-dev \
    libboost-all-dev \
    wget \
    gawk
#ADD deps/bwa_0.7.15-5_amd64.deb /tmp/bwa.deb
#RUN dpkg -i /tmp/bwa.deb


RUN wget -O vgan-main.zip https://github.com/grenaud/vgan/archive/refs/tags/v1.0.2.zip

RUN unzip vgan-main.zip -d /vgan
WORKDIR /vgan/
RUN pwd
RUN ls
WORKDIR /vgan/vgan-1.0.2/
RUN pwd
RUN ls
WORKDIR /vgan/vgan-1.0.2/src/
RUN pwd
RUN ls
RUN make
