GLIBC_LESS_2_17=$(shell ldd --version | head -1 | awk '{print $$NF<2.17?"1":"0"}')

PZCL_STATIC_LIB ?=1
OBJ_DIR ?= obj

ifndef PZSDK_PATH
$(warning !!!!!! WARNING !!!!!!  PZSDK_PATH is not set. We use /opt/pezy/pzsdk.)
PZSDK_PATH = /opt/pezy/pzsdk
endif

PZSDK_INC_PATH ?= $(PZSDK_PATH)/inc
PZSDK_LIB_PATH ?= $(PZSDK_PATH)/lib

INC_DIR += -I$(PZSDK_INC_PATH)

ifeq ($(PZCL_STATIC_LIB),1)
PZCL_LIB ?= $(PZSDK_LIB_PATH)/libpzcl.a
else
PZCL_LIB ?= -lpzcl
endif

LIB_DIR += -L$(PZSDK_LIB_PATH)
LIBS    += $(PZCL_LIB) -lpthread -ldl -lstdc++ -lrt

PZCL_KERNEL_DIRS ?= 

CXX ?= g++
CC ?= gcc
FC ?= gfortran

CCOPT  ?= -D__LINUX__ -O2 -DNDEBUG
COPT   ?= -std=c99 -D__LINUX__ -O2 -DNDEBUG

###########################################
# Preprocess
###########################################
CPPOBJS= ${patsubst %.cpp, ${OBJ_DIR}/%.o, ${notdir ${filter-out %.cc %.cxx,${CPPSRC}}}}
CCOBJS= ${patsubst %.cc, ${OBJ_DIR}/%.o, ${notdir ${filter %.cc,${CPPSRC}}}}
CXXOBJS= ${patsubst %.cxx, ${OBJ_DIR}/%.o, ${notdir ${filter %.cxx,${CPPSRC}}}}
COBJS= ${patsubst %.c, ${OBJ_DIR}/%.o, ${notdir ${CSRC}}}
FCOBJS= ${patsubst %.F, ${OBJ_DIR}/%.o, ${notdir ${FCSRC}}}
FFOBJS= ${patsubst %.f, ${OBJ_DIR}/%.o, ${notdir ${FFSRC}}}

OBJS   = ${CPPOBJS} ${CCOBJS} ${CXXOBJS} ${COBJS} ${FCOBJS} ${FFOBJS}
DEPS   = ${OBJS:.o=.d} 

all:make_pzcl_kernel ${OBJ_DIR} ${OBJS} ${TARGET} 

###########################################
# embedded objs
###########################################
ifdef PZ_EMBED_TARGETS

PZ_NUM_EMBED=$(shell expr $(words $(PZ_EMBED_TARGETS)) / 2 - 1)
PZ_NUM_EMBED_SEQ=$(shell seq 0 $(PZ_NUM_EMBED))

# 
# PZ_EMBED_RULE_TEMPLATE $(1)=<pz file> $(2)=<x86 symbol>
#
define PZ_EMBED_RULE_TEMPLATE

$(1): $(shell dirname $(1))

${OBJ_DIR}/$(strip $(2).cpp):$(shell dirname $(1)) | ${OBJ_DIR}
	@echo -e "#include <utility>\n#include <stddef.h>"	>$$@ && \
	echo "static const unsigned char data[] = {"	>>$$@ && \
	hexdump -ve '/1 "0x%02x,"' $(1) >>$$@ && \
	echo "};"	>>$$@ && \
	echo "std::pair<const void*,size_t> $(2)() { return std::make_pair((const void*)data, sizeof(data)); }"	>>$$@

.DELETE_ON_ERROR: ${OBJ_DIR}/$(strip $(2).cpp)


${OBJ_DIR}/$(strip $(2).o):${OBJ_DIR}/$(strip $(2).cpp)
	${CXX} -o $$@ -c ${CCOPT} $$<

PZ_EMBED_CXXSRCS+=${OBJ_DIR}/$(strip $(2).cpp)
PZ_EMBED_OBJS+=${OBJ_DIR}/$(strip $(2).o)

endef

$(foreach idx,$(PZ_NUM_EMBED_SEQ), \
	 $(eval $(call PZ_EMBED_RULE_TEMPLATE,\
			$(word $(shell expr $(idx) '*' 2 + 1), $(PZ_EMBED_TARGETS)), \
			$(word $(shell expr $(idx) '*' 2 + 2), $(PZ_EMBED_TARGETS)))))

endif

clean : clean_pzcl_kernel
	-rm -f ${TARGET} ${OBJS} ${PZ_EMBED_OBJS} ${PZ_EMBED_CXXSRCS} $(DEPS) 
	-rm -rf ${OBJ_DIR}

${TARGET}: ${OBJS} ${PZ_EMBED_OBJS}
	${FC} -o ${TARGET} ${OBJS} ${PZ_EMBED_OBJS} ${LIB_DIR} ${LIBS} ${LDOPT}

${CPPOBJS}: ${OBJ_DIR}/%.o : %.cpp | ${OBJ_DIR}
	${CXX} ${CCOPT} ${INC_DIR} -c $< -o $@
	${CXX} ${CCOPT} ${INC_DIR} -MM -MT$@ $< > ${OBJ_DIR}/$*.d

${CXXOBJS}: ${OBJ_DIR}/%.o : %.cxx | ${OBJ_DIR}
	${CXX} ${CCOPT} ${INC_DIR} -c $< -o $@
	${CXX} ${CCOPT} ${INC_DIR} -MM -MT$@ $< > ${OBJ_DIR}/$*.d

${CCOBJS}: ${OBJ_DIR}/%.o : %.cc | ${OBJ_DIR}
	${CXX} ${CCOPT} ${INC_DIR} -c $< -o $@
	${CXX} ${CCOPT} ${INC_DIR} -MM -MT$@ $< > ${OBJ_DIR}/$*.d

${COBJS}: ${OBJ_DIR}/%.o : %.c | ${OBJ_DIR}
	${CC} ${COPT} ${INC_DIR} -c $< -o $@
	${CC} ${COPT} ${INC_DIR} -MM -MT$@ $< > ${OBJ_DIR}/$*.d

${FCOBJS}: ${OBJ_DIR}/%.o : %.F | ${OBJ_DIR}
	${FC} ${FCOPT} -c $< -o $@
#	${FC} ${FCOPT} ${INC_DIR} -MM -MT$@ $< > ${OBJ_DIR}/$*.d

${FFOBJS}: ${OBJ_DIR}/%.o : %.f | ${OBJ_DIR}
	${FC} ${FCOPT} -c $< -o $@
#	${FC} ${FCOPT} ${INC_DIR} -MM -MT$@ $< > ${OBJ_DIR}/$*.d

${OBJ_DIR}:
	@if [ ! -d ${OBJ_DIR} ]; then \
		echo "mkdir ${OBJ_DIR}"; mkdir ${OBJ_DIR}; \
	fi


$(PZCL_KERNEL_DIRS): FORCE
	@echo "" ; \
	echo "*******************************************************" ; \
	echo "* making in ./$@" ; \
	echo "*******************************************************" ; \
	$(MAKE) -C $@ ; \
	if [ $$? -ne 0 ] ; then \
		echo "!!!!!! ERROR !!!!!!" ; exit 1; \
	fi

make_pzcl_kernel : $(PZCL_KERNEL_DIRS)

clean_pzcl_kernel:
	@for d in $(PZCL_KERNEL_DIRS); \
	do \
	(cd $$d && \
	echo "" && \
	echo "*******************************************************" && \
	echo "* clean in ./$$d" && \
	echo "*******************************************************" && \
	$(MAKE) clean); \
	if [ $$? -ne 0 ]; then \
		echo "!!!!!! ERROR !!!!!!" ; exit 1; \
	fi \
	done

.PHONY: FORCE all clean make_pzcl_kernel clean_pzcl_kernel

-include $(DEPS)


