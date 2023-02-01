#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#

.PHONY: all clean obj bin lib doc

all: main library
	@echo "[SKSNSim] Done!"


ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif
include $(SKOFL_ROOT)/config.gmk

CXX=g++
CXXFLAGS += -DNO_EXTERN_COMMON_POINTERS #-DDEBUG
# if you want to use lates neutrino oscillation parameter, please comment out next line
CXXFLAGS += -DORIGINAL_NUOSCPARAMETER
FC=gfortran
FCFLAGS += -w -fPIC -lstdc++

LOCAL_INC	+= -I./include/

LOCAL_LIBS	= -L$(SKOFL_LIBDIR) -lsnlib_1.0 -lsnevtinfo -lsollib_4.0 -lsklowe_7.0 -llibrary 

LDFLAGS = $(LOCAL_LIBS) $(LOCAL_INC)
#INCROOT=-I$(ROOTSYS)/include/
#LIBSROOT=-L$(ROOTSYS)/lib/ -lCint -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui

LN = ln -sf

#
#  Objects
#

SRCS = $(wildcard src/*.cc)
SRCS += $(wildcard src/*.F)
OBJS = $(patsubst src/%.cc, obj/%.o, $(filter %.cc, $(SRCS)))
OBJS += $(patsubst src/%.F, obj/%.o, $(filter %.F, $(SRCS)))

MAINSRCS = $(wildcard *.cc)
MAINSRCS += $(wildcard *.F)
MAINOBJS = $(patsubst %.cc, obj/%.o, $(MAINSRCS))
MAINBINS = $(patsubst %.cc, bin/%, $(MAINSRCS))
SKSNSIMLIBOBJS = $(filter obj/SKSNSim%, $(OBJS))

main: bin obj bin/main_snburst bin/main_dsnb bin/main_snburst_prev bin/main_dsnb_prev bin/main_dsnb_new bin/main_snburst_new

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
	@echo "LDFLAGS:      "$(LDFLAGS)
	@echo "LDLIBS:       "$(LDLIBS)

obj/%.o: %.cc 
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.cc
	@echo "[SKSNSim] Building $* ..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.F
	@echo "[SKSNSim] Building FORTRAN code: $*..."
	@$(FC) $(FCFLAGS) -c $< -o $@

bin/main_snburst_prev: obj/main_snburst.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

bin/main_dsnb_prev: obj/main_dsnb.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

bin/main_dsnb_new: obj/main_dsnb_new.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

bin/main_snburst_new: obj/main_snburst_new.o $(OBJS)
	@echo "[SKSNSim] Building executable:	$@..."
	@LD_RUN_PATH=$(SKOFL) $(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

bin/main_snburst: bin/main_snburst_prev
	@${LN} main_snburst_prev $@
	#@${LN} main_snburst_new $@

bin/main_dsnb: bin/main_dsnb_prev
	@${LN} main_dsnb_prev $@
	#@${LN} main_dsnb_new $@

lib/libSKSNSim.so: $(SKSNSIMLIBOBJS)
	@echo "[SKSNSim] Making shared library: $@..."
	LD_RUN_PATH=$(SKOFL) $(CXX) $(LDFLAGS) $(CXXFLAGS) -shared -o $@ $(SKSNSIMLIBOBJS)  $(LDLIBS) #$(SKOFL_LIBDIR)/libsollib_4.0.a $(SKOFL_LIBDIR)/libsnlib_1.0.a $(SKOFL_LIBDIR)/libiolib.a

obj bin lib:
	@mkdir -p $@

clean: 
	$(RM) -r *.o *~ lib/* src/*~ include/*~ core obj/* bin/* obj bin lib $(TARGET)

-include main_dsnb_new.d
