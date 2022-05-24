#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#
CXX=g++

ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif

include $(SKOFL_ROOT)/config.gmk


LOCAL_INC	+= -I./include/

LOCAL_LIBS	= $(OBJS) \
		-lsnlib_1.0 -lsnevtinfo -lsollib_4.0 -lsklowe_7.0 -llibrary 

INCROOT=-I$(ROOTSYS)/include/

#
#  Objects
#
SRCS = $(wildcard src/*.cc src/*.F )
OBJS = $(patsubst src/%.cc,  obj/%.o, $(SRCS))
#OBJS = ${SRCS:.cc=.o:.F=.o} 

CFLAGS += -DNO_EXTERN_COMMON_POINTERS

TARGET = main_dsnb main_snburst

all: $(TARGET)

main_snburst: main_snburst.o $(OBJS)
	@LD_RUN_PATH=$(SKOFL_LIBDIR) $(CXX) $(CXXFLAGS) -o $@ main_snburst.o $(OBJS) $(LDLIBS)

main_dsnb: main_dsnb.o $(OBJS)
	@LD_RUN_PATH=$(SKOFL_LIBDIR) $(CXX) $(CXXFLAGS) $(LDLIBS) -o $@ $<

%.o: %.cc
	@$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.F
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.cc obj
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj:
	@mkdir -p $@

clean: 
	$(RM) -r *.o *~ src/*~ include/*~ core obj/* bin/* obj bin $(TARGET)
