# Root wrapper Makefile for CMake transition

BUILD_DIR ?= build
NPROCS := $(shell nproc 2>/dev/null || nprocs 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

.PHONY: all build clean install enzo inits enzohop ring anyl P-GroupFinder

all: build

build:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --parallel $(NPROCS)

clean:
	@rm -rf $(BUILD_DIR)

enzo:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target enzo --parallel $(NPROCS)

inits:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target inits --parallel $(NPROCS)

enzohop:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target enzohop --parallel $(NPROCS)

ring:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target ring --parallel $(NPROCS)

anyl:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target anyl --parallel $(NPROCS)

P-GroupFinder:
	@cmake -B $(BUILD_DIR) -S .
	@cmake --build $(BUILD_DIR) --target P-GroupFinder --parallel $(NPROCS)
