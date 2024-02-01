#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#

.PHONY: all clean obj bin lib doc test

all: main library
	@echo "[SKSNSim] Done!"


ifdef SKOFL_ROOT
include $(SKOFL_ROOT)/config.gmk
endif

LOCAL_INC	+= -I./include/

CXX=g++
CXXFLAGS +=$(shell root-config --cflags --libs) -fPIC -lstdc++
CXXFLAGS += -DNO_EXTERN_COMMON_POINTERS #-DDEBUG
CXXFLAGS += $(LOCAL_INC)
ifdef SKOFL_ROOT
CXXFLAGS += -DSKINTERNAL
endif
# if you want to use lates neutrino oscillation parameter, please comment out next line
#CXXFLAGS += -DORIGINAL_NUOSCPARAMETER
FC=gfortran
FCFLAGS += -w -fPIC -lstdc++

ifndef SKOFL_ROOT
LDLIBS=$(shell root-config --libs)
endif

ifdef SKOFL_ROOT
LDFLAGS = $(LOCAL_LIBS) $(LOCAL_INC)
LOCAL_LIBS	= -L$(SKOFL_LIBDIR) -lsnlib_1.0 -lsnevtinfo -lsollib_4.0 -lsklowe_7.0 -lwtlib_5.1 -llibrary 
endif
#INCROOT=-I$(ROOTSYS)/include/
#LIBSROOT=-L$(ROOTSYS)/lib/ -lCint -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui

LN = ln -sf



#
#  Objects
#

SRCS = $(wildcard src/*.cc)
ifdef SKOFL_ROOT
SRCS += $(wildcard src/*.F)
endif
OBJS = $(patsubst src/%.cc, obj/%.o, $(filter %.cc, $(SRCS)))
ifdef SKOFL_ROOT
OBJS += $(patsubst src/%.F, obj/%.o, $(filter %.F, $(SRCS)))
endif

MAINSRCS = $(wildcard *.cc)
MAINSRCS += $(wildcard *.F)
MAINOBJS = $(patsubst %.cc, obj/%.o, $(MAINSRCS))
MAINBINS = $(patsubst %.cc, bin/%, $(MAINSRCS))
SKSNSIMLIBOBJS = $(filter obj/SKSNSim%, $(OBJS))
SKSNSIMLIBOBJS += $(filter obj/elapseday%, $(OBJS))

main: bin obj bin/main_snburst bin/main_dsnb

library: lib lib/libSKSNSim.so

doc: doc/guide_sksnsim.pdf

doc/guide_sksnsim.pdf: doc/guide_sksnsim.texi
	cd doc && pwd && texi2any --pdf guide_sksnsim.texi && cd ../

test:
	@echo "MAINSRCS      "$(MAINSRCS)
	@echo "MAINOBJS      "$(MAINOBJS)
	@echo "MAINBINS      "$(MAINBINS)
	@echo "SRCS:         "$(SRCS)
	@echo "OBJS:         "$(OBJS)
	@echo "INC:          "$(INC)
	@echo "CXXFLAGS:     "$(CXXFLAGS)
	@echo "LDFLAGS:      "$(LDFLAGS)
	@echo "LDLIBS:       "$(LDLIBS)

obj/%.o: %.cc 
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.cc
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.F
	@echo "[SKSNSim] Building FORTRAN code: $*..."
	@$(FC) $(FCFLAGS) -c $< -o $@

bin/main_dsnb: obj/main_dsnb.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

bin/main_snburst: obj/main_snburst.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)


lib/libSKSNSim.so: $(SKSNSIMLIBOBJS)
	@echo "[SKSNSim] Making shared library: $@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(LDFLAGS) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

obj bin lib:
	@mkdir -p $@

clean: 
	$(RM) -r *.o *~ lib/* src/*~ include/*~ core obj/* bin/* obj bin lib $(TARGET)

