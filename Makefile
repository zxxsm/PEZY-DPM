PZSDK_PATH ?= /opt/pezy/pzsdk
export PZSDK_PATH

PWD := $(shell pwd)

MODE ?= mpi
PRECISION ?= float
KERNEL_VERSION ?= 8th

vpath %.cpp src
vpath %.hpp src
vpath %.F src
vpath %.f src

BASE_CXXFLAGS := -O2 -Wall -D__LINUX__ -DDEBUG -std=c++11
BASE_FFLAGS   := -O2 -mcmodel=medium -fallow-argument-mismatch

ifeq ($(PRECISION),double)
	PRECISION_FLAGS := -DDPM_USE_DOUBLE
else ifeq ($(PRECISION),float)
	PRECISION_FLAGS :=
else
	$(error PRECISION must be float or double)
endif

ifeq ($(MODE),mpi)
	DEFAULT_MAKE := compile/cpu.mk
	TARGET := dpm.mpi
	CPPSRC := Cfun.cpp
	CCOPT := $(BASE_CXXFLAGS) -DMPI_P $(PRECISION_FLAGS)
	FCSRC := dpm.F
	FFSRC := time.f libmath.f libpenmath.f libeloss.f libphoton.f libgeom.f getnam.f
	FCOPT := $(BASE_FFLAGS) -DMPI_P
	CXX := mpicxx
	FC := mpif90
	LIBS := -lstdc++
	PZCL_KERNEL_DIRS :=
else ifeq ($(MODE),openmp)
	DEFAULT_MAKE := compile/cpu.mk
	TARGET := dpm.omp
	CPPSRC := Cfun.cpp
	CCOPT := $(BASE_CXXFLAGS) -DOMP_P -fopenmp $(PRECISION_FLAGS)
	FCSRC := dpm.F
	FFSRC := time.f libmath.f libpenmath.f libeloss.f libphoton.f libgeom.f getnam.f
	FCOPT := $(BASE_FFLAGS) -DOMP_P -fopenmp
	CXX := g++
	FC := gfortran
	LDOPT := -fopenmp
	LIBS := -lstdc++
	PZCL_KERNEL_DIRS :=
else ifeq ($(MODE),mpi_sc3s)
	DEFAULT_MAKE := compile/default_pzcl_host.mk
	TARGET := dpm.sc_mpi
	CPPSRC := Cfun.cpp SC3s.cpp
	CCOPT := $(BASE_CXXFLAGS) -DMPI_P -DPEZY_SCX $(PRECISION_FLAGS)
	FCSRC := dpm.F
	FFSRC := time.f libmath.f libpenmath.f libeloss.f libphoton.f libgeom.f getnam.f
	FCOPT := $(BASE_FFLAGS) -DMPI_P
	CXX := mpicxx
	FC := mpif90
	PZCL_STATIC_LIB ?= 0
	PZCL_KERNEL_DIRS := kernel
else
	$(error MODE must be mpi, openmp, or mpi_sc3s)
endif

INC_DIR ?=
LIB_DIR ?=

# supported architecture: sc1-64, sc2, sc3, sc3s
PZC_TARGET_ARCH ?= sc3s
export PZC_TARGET_ARCH KERNEL_VERSION PRECISION

.DEFAULT_GOAL := all

.PHONY: mpi openmp mpi_sc3s help

mpi:
	$(MAKE) MODE=mpi PRECISION=$(PRECISION)

openmp:
	$(MAKE) MODE=openmp PRECISION=$(PRECISION)

mpi_sc3s:
	$(MAKE) MODE=mpi_sc3s PRECISION=$(PRECISION) KERNEL_VERSION=$(KERNEL_VERSION)

help:
	@echo "Build examples:"
	@echo "  make mpi PRECISION=float"
	@echo "  make mpi PRECISION=double"
	@echo "  make openmp PRECISION=float"
	@echo "  make openmp PRECISION=double"
	@echo "  make mpi_sc3s PRECISION=float KERNEL_VERSION=8th PZC_TARGET_ARCH=sc3s"
	@echo "  make mpi_sc3s PRECISION=double KERNEL_VERSION=4th PZC_TARGET_ARCH=sc3s"

include $(DEFAULT_MAKE)

run:
	@./$(TARGET) < data/input/dpm.in
