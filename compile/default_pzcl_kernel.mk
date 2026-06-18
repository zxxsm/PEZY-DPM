ifndef PZSDK_PATH
$(warning !!!!!! WARNING !!!!!!  PZSDK_PATH is not set. We use /opt/pezy/pzsdk.)
PZSDK_PATH = /opt/pezy/pzsdk
endif

LLVM_PATH  ?= $(PZSDK_PATH)

PZC_RELEASE     ?= 0
PZC_INC_DIR     ?= 
PZC_TARGET_ARCH ?= sc2

ifeq ($(PZC_TARGET_ARCH), sc3s) 
	PZRUNTIME ?= $(PZSDK_PATH)/bin/pzcrt-sc3.pzo
	PZKF      ?= $(PZSDK_PATH)/bin/pzkf-sc3.pzo
	PZATOMIC  ?= $(PZSDK_PATH)/bin/pzatomic-sc3.pzo
else
ifeq ($(PZC_TARGET_ARCH), sc4s) 
	PZRUNTIME ?= $(PZSDK_PATH)/bin/pzcrt-sc3.pzo
	PZKF      ?= $(PZSDK_PATH)/bin/pzkf-sc3.pzo
	PZATOMIC  ?= $(PZSDK_PATH)/bin/pzatomic-sc3.pzo
else
	PZRUNTIME ?= $(PZSDK_PATH)/bin/pzcrt-${PZC_TARGET_ARCH}.pzo
	PZKF      ?= $(PZSDK_PATH)/bin/pzkf-${PZC_TARGET_ARCH}.pzo
	PZATOMIC  ?= $(PZSDK_PATH)/bin/pzatomic-${PZC_TARGET_ARCH}.pzo
endif
endif

CLANG_OPT ?= -O3
LLC_OPT   ?= -O3
OPT_OPT   ?= -O3
MC_OPT    ?=
LINK_OPT  ?=

OBJ_DIR ?= obj

# FIXME These variables will be removed. Do not use!
PZSDK_INC_DIR ?= $(PZSDK_PATH)/inc
PZLINK_DIR    ?= $(PZSDK_PATH)/bin

CLANG  = $(LLVM_PATH)/bin/clang
IRLINK = $(LLVM_PATH)/bin/llvm-link
NM     = $(LLVM_PATH)/bin/llvm-nm
OPT    = $(LLVM_PATH)/bin/opt
LLC    = $(LLVM_PATH)/bin/llc
MC     = $(LLVM_PATH)/bin/llvm-mc
LINK   = $(PZLINK_DIR)/pzlink

ifneq ($(PZC_RELEASE), 0)
	LLC_OPT += -disable-pz-debug-check
endif

PZLIB += $(PZKF) $(PZATOMIC)

######################################################
#  Preprocess
######################################################
ifeq ($(suffix $(TARGET)),.pzo)
TARGET_MAKE_PZ =
TARGET_MAKE_PZO=yes
TARGET_PZ =$(TARGET:.pzo=.pz)
TARGET_PZO=$(TARGET)
else
ifeq ($(suffix $(TARGET)),.pz)
TARGET_MAKE_PZ  =yes
TARGET_MAKE_PZO?=
TARGET_PZ =$(TARGET)
TARGET_PZO=$(TARGET:.pz=.pzo)
else
$(error The suffix of TARGET should be either .pz or .pzo)
endif
endif

TARGET_ALL=
ifdef TARGET_MAKE_PZO
TARGET_ALL+= $(TARGET_PZO)
endif
ifdef TARGET_MAKE_PZ
TARGET_ALL+= $(TARGET_PZ)
endif

PZCIRS   = ${patsubst %.pzc,${OBJ_DIR}/%.ir,${notdir ${PZCSRC}}} 
PZCIROBJ = ${patsubst %.pz,${OBJ_DIR}/%.iro,${notdir ${TARGET_PZ}}} 
PZCSYM   = $(PZCIROBJ:.iro=.sym)
PZCOBJ   = $(PZCIROBJ:.iro=.o)
PZASM    = $(PZCIROBJ:.iro=.pzs)
PZPPOUT  = $(PZCIROBJ:.iro=.ll)
PZHEX    = $(PZCIROBJ:.iro=.hex)
PZDEPS   = ${PZCIROBJ:.iro=.d} 

######################################################
# Targets
######################################################
all: ${TARGET_ALL}


${OBJ_DIR}:
	@if [ ! -d ${OBJ_DIR} ]; then \
		echo "mkdir ${OBJ_DIR}"; mkdir ${OBJ_DIR}; \
	fi

${PZCIRS} : ${OBJ_DIR}/%.ir : %.pzc | ${OBJ_DIR}
	${CLANG} -target pz64 -mcpu=${PZC_TARGET_ARCH} ${CLANG_OPT} ${PZC_INC_DIR} -x c++ -I$(PZSDK_INC_DIR) -fno-threadsafe-statics -fno-exceptions -emit-llvm -c $< -o $@
	${CLANG} -target pz64 -mcpu=${PZC_TARGET_ARCH} ${CLANG_OPT} ${PZC_INC_DIR} -x c++ -I$(PZSDK_INC_DIR) -MM -MT$@ $< > ${OBJ_DIR}/$*.d

${PZASM} : ${PZCIRS}
	${IRLINK} ${PZLIB} ${PZCIRS} -o ${PZCIROBJ}  # Make Obj
	${NM} --defined-only --extern-only ${PZCIROBJ} | awk '/ T _Z[0-9]*pz[c,1]_/ {print $$(NF)}' > ${PZCSYM}
	${OPT} -internalize -internalize-public-api-file=${PZCSYM} ${OPT_OPT} ${PZCIROBJ} -o=${PZCOBJ}
	${IRLINK} ${PZRUNTIME} ${PZCOBJ} -o ${PZCOBJ}
	${LLC} -march=pz64 -mcpu=${PZC_TARGET_ARCH} --no-integrated-as ${LLC_OPT} $(PZCOBJ) -o $(PZASM)

${TARGET_PZO} : ${PZCIRS}
	${IRLINK} ${PZCIRS} -o $@

${TARGET_PZ} : ${PZASM}
	$(CLANG) -target pz64 -mcpu=${PZC_TARGET_ARCH} -x c++ -I$(PZSDK_INC_DIR) -E -DPZ_LLVMMC_ASM $(PZASM) -o $(PZPPOUT)
	$(MC) -filetype=obj -arch=pz64 -mcpu=$(PZC_TARGET_ARCH) $(MC_OPT) $(PZPPOUT) -o $(PZHEX)
	-rm -f $(TARGET_PZ)
	$(LINK) -v $(PZHEX) -o $(TARGET_PZ) $(LINK_OPT)

clean : 
	-rm -f ${PZCIRS} ${TARGET_PZO} ${TARGET_PZ}
	-rm -rf ${OBJ_DIR}

.PHONY: all clean

-include $(PZDEPS)
