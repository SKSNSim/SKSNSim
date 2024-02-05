#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#

.PHONY: all clean obj bin lib doc

TARGET = reweight validate_nuelastic_xsec
all: $(TARGET)
	@echo "[SKSNSim] Done!"

ifdef SKOFL_ROOT
include $(SKOFL_ROOT)/config.gmk
endif

LOCAL_INC	+=  -I$(SKSNSIMINSTALLDIR)/include

CXX=g++
FC=gfortran
LN = ln -sf

CXXFLAGS +=$(shell root-config --cflags --libs) -fPIC -lstdc++ -lgsl -lgslcblas -lm -lboost_filesystem
CXXFLAGS += -DNO_EXTERN_COMMON_POINTERS #-DDEBUG
CXXFLAGS += $(LOCAL_INC)
ifdef SKOFL_ROOT
CXXFLAGS += -DSKINTERNAL
endif
# if you want to use lates neutrino oscillation parameter, please comment out next line
#CXXFLAGS += -DORIGINAL_NUOSCPARAMETER

FCFLAGS += -w -fPIC -lstdc++

ifndef SKOFL_ROOT
LDLIBS=$(shell root-config --libs)
endif

ifdef SKOFL_ROOT
LDFLAGS = $(LOCAL_LIBS) $(LOCAL_INC)
LOCAL_LIBS	= -L$(SKOFL_LIBDIR) -lsnlib_1.0 -lsnevtinfo -lsollib_4.0 -lsklowe_7.0 -lwtlib_5.1 -llibrary 
endif

LOCAL_LIBS += -lSKSNSim
LDLIBS += -L$(SKSNSIMINSTALLDIR)/lib $(LOCAL_LIBS)

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
	@echo "CXXFLAGS:     "$(CXXFLAGS)

%.o: %.cc 
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@


%: %.o 
	@echo "[SKSNSim] Building executable:	$@..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

clean: 
	$(RM) -r *.o *~ $(TARGET)

