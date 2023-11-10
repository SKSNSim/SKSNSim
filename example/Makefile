#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#

.PHONY: all clean obj bin lib doc

TARGET = reweight
all: $(TARGET)
	@echo "[SKSNSim] Done!"


ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif
include $(SKOFL_ROOT)/config.gmk


CXX=g++
CXXFLAGS += -DNO_EXTERN_COMMON_POINTERS #-DDEBUG
# if you want to use lates neutrino oscillation parameter, please comment out next line
#CXXFLAGS += -DORIGINAL_NUOSCPARAMETER
FC=gfortran
FCFLAGS += -w -fPIC -lstdc++

LOCAL_INC	+= -I$(SKSNSIM_ROOTDIR)/include

LOCAL_LIBS	= -L$(SKOFL_LIBDIR) -lsnlib_1.0 -lsnevtinfo -lsollib_4.0 -lsklowe_7.0 -lwtlib_5.1 -llibrary -L$(SKSNSIM_ROOTDIR)/lib -lSKSNSim

LDFLAGS = $(LOCAL_LIBS) $(LOCAL_INC)
#INCROOT=-I$(ROOTSYS)/include/
#LIBSROOT=-L$(ROOTSYS)/lib/ -lCint -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui

LN = ln -sf



#
#  Objects
#

test:
	@echo "MAINSRCS      "$(MAINSRCS)
	@echo "MAINOBJS      "$(MAINOBJS)
	@echo "MAINBINS      "$(MAINBINS)
	@echo "SRCS:         "$(SRCS)
	@echo "OBJS:         "$(OBJS)
	@echo "INC:          "$(INC)
	@echo "LDFLAGS:      "$(LDFLAGS)
	@echo "LDLIBS:       "$(LDLIBS)

%.o: %.cc 
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@


reweight: reweight.o 
	@echo "[SKSNSim] Building executable:	$@..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

clean: 
	$(RM) -r *.o *~ $(TARGET)

