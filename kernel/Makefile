PZSDK_PATH ?= /opt/pezy/pzsdk
DEFAULT_MAKE = ../compile/default_pzcl_kernel.mk

PZC_TARGET_ARCH ?= sc3s
KERNEL_VERSION ?= 8th
PRECISION ?= float

ifeq ($(filter $(KERNEL_VERSION),4th 8th),)
	$(error KERNEL_VERSION must be 4th or 8th)
endif

ifeq ($(PRECISION),double)
	CLANG_OPT += -DDPM_USE_DOUBLE
else ifneq ($(PRECISION),float)
	$(error PRECISION must be float or double)
endif

TARGET = kernel.pz
PZCSRC = kernel_$(KERNEL_VERSION).pzc

vpath %.pzc ../pzc

include $(DEFAULT_MAKE)
