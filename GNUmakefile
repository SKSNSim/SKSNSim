#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#

.PHONY: all clean obj

TARGET = main_dsnb main_snburst main_dsnb_new main_snburst_new

all: main
	@echo "[SKSNSim] Done!"


ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif
include $(SKOFL_ROOT)/config.gmk

CXX=g++
CXXFLAGS += -g3 -DNO_EXTERN_COMMON_POINTERS -DDEBUG
FC=gfortran
FCFLAGS += -w -fPIC -lstdc++

LOCAL_INC	+= -I./include/

LOCAL_LIBS	= $(OBJS)\
		-lsnlib_1.0 -lsnevtinfo -lsollib_4.0 -lsklowe_7.0 -llibrary 

LDFLAGS = $(LOCAL_LIBS) $(LOCAL_INC)
#INCROOT=-I$(ROOTSYS)/include/
#LIBSROOT=-L$(ROOTSYS)/lib/ -lCint -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui

#
#  Objects
#

OBJS = $(patsubst src/%.cc, obj/%.o, $(wildcard src/*.cc))
OBJS += $(patsubst src/%.F, obj/%.o, $(wildcard src/*.F))


#SRCS = $(wildcard src/*.cc)
#OBJS = $(patsubst src/%.cc, obj/%.o, $(SRCS))
#FORTRANSRCS = $(wildcard src/*.F)
#FORTRANOBJS = $(patsubst src/%.F, obj/%.o, $(FORTRANSRCS))
#SRCS = $(sort $(shell find src -name '*.cc'))
#OBJS = $(patsubst src/%, obj/%.o, $(basename $(SRCS)))
#FORTRANSRCS = $(sort $(shell find src -name '*.F'))
#FORTRANOBJS = $(patsubst src/%, obj/%.o, $(basename $(FORTRANSRCS)))
#INC := $(addprefix -I , $(sort $(dir $(shell find include -name '*.hh'))))

MAINSRCS = $(wildcard *.cc)
MAINOBJS = $(patsubst %.cc, obj/%.o, $(MAINSRCS))
MAINBINS = $(patsubst %.cc, bin/%, $(MAINSRCS))

main: obj bin bin/main_snburst bin/main_dsnb bin/main_dsnb_new bin/main_snburst_new

test:
	@echo "MAINSRCS      "$(MAINSRCS)
	@echo "MAINOBJS      "$(MAINOBJS)
	@echo "MAINBINS      "$(MAINBINS)
	@echo "SRCS:         "$(SRCS)
	@echo "OBJS:         "$(OBJS)
	@echo "FORTRAN SRCS: "$(FORTRANSRCS)
	@echo "FORTRAN OBJS: "$(FORTRANOBJS)
	@echo "INC:          "$(INC)
	@echo "LDFLAGS:      "$(LDFLAGS)

obj/main_snburst.o: main_snburst.cc
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/main_dsnb.o: main_dsnb.cc
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/main_dsnb_new.o: main_dsnb_new.cc
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@
#main_dsnb_new.d: main_dsnb_new.cc
#	$(CXX) $(CXXFLAGS) -MM -o $@ $<

obj/main_snburst_new.o: main_snburst_new.cc
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.cc
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.F
	@echo "[SKSNSim] Building FORTRAN code: $*..."
	@$(FC) $(FCFLAGS) -c $< -o $@

bin/main_snburst: obj/main_snburst.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

bin/main_dsnb: obj/main_dsnb.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

bin/main_dsnb_new: obj/main_dsnb_new.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

bin/main_snburst_new: obj/main_snburst_new.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

obj bin:
	@mkdir -p $@

clean: 
	$(RM) -r *.o *~ src/*~ include/*~ core obj/* bin/* obj bin $(TARGET)

-include main_dsnb_new.d
