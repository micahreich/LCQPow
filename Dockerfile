FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-lc"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    pkg-config \
    ca-certificates \
    python3 \
    python3-dev \
    python3-pip \
    libeigen3-dev \
    sudo \
    vim \
  && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir --upgrade pip \
 && python3 -m pip install --no-cache-dir numpy

WORKDIR /workspace
COPY . /workspace

# Build configuration (can be overridden with --build-arg BUILD_TYPE=Debug)
ARG BUILD_TYPE=Release
ENV SRC_DIR=/workspace
ENV BUILD_DIR=/workspace/build

RUN set -euo pipefail; \
    if [[ ! -d "$SRC_DIR" ]]; then \
      echo "ERROR: expected repo at $SRC_DIR"; \
      exit 1; \
    fi; \
    mkdir -p "$BUILD_DIR"; \
    echo "Configuring in: $BUILD_DIR"; \
    cmake -S "$SRC_DIR" -B "$BUILD_DIR" \
      -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
      -DEXAMPLES=ON \
      -DPYTHON_INTERFACE=ON \
      -DMATLAB_INTERFACE=OFF \
      -DDOCUMENTATION=OFF \
      -DUNIT_TESTS=ON \
      -DPROFILING=OFF \
      -DQPOASES_SCHUR=OFF; \
    echo "Building..."; \
    cmake --build "$BUILD_DIR" -j"$(nproc)"

