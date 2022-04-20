#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#


ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif

include $(SKOFL_ROOT)/config.gmk

LOCAL_INC	+= -I./include/

LOCAL_LIBS	= $(OBJS) \
		-lsnlib_1.0 -lsnevtinfo -lsollib_4.0 -llibrary

#
#  Objects
#
SRCS = $(wildcard src/*.cc)
OBJS = $(patsubst src/%.cc, obj/%.o, $(SRCS))
#OBJS = ${SRCS:.cc=.o:.F=.o} 

CFLAGS += -DNO_EXTERN_COMMON_POINTERS

TARGET = main_srn main_snburst

all: $(TARGET)

main: main.o $(OBJS)
	@LD_RUN_PATH=$(LIBDIR) $(CXX) $(CXXFLAGS) -o $@ main.o $(OBJS) $(LDLIBS)

main_srn: main_srn.o $(OBJS)
	@LD_RUN_PATH=$(LIBDIR) $(CXX) $(CXXFLAGS) $(LDLIBS) -o $@ $<

%.o: %.cc
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.cc obj
	@echo "[SKSNSim] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj:
	@mkdir -p $@

clean: 
	$(RM) -r *.o *~ src/*~ include/*~ core obj/* bin/* obj bin $(TARGET)
