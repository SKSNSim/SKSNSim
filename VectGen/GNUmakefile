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

all: bin/CalcExpect

bin/CalcExpect: obj/main.o $(OBJS) bin data
	@LD_RUN_PATH=$(LIBDIR) $(CXX) $(CXXFLAGS) -o $@ obj/main.o $(OBJS) $(LDLIBS)

obj/main.o: main.cc obj
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: src/%.cc obj
	@echo "[mccomb] Building $* ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

obj bin data:
	@mkdir $@

clean: 
	$(RM) -r *.o *~ src/*~ include/*~ core obj/* bin/* obj bin

GenEvent: obj/GenEvent.o $(OBJS) bin
	LD_RUN_PATH=$(LIBDIR) $(CXX) $(CXXFLAGS) -o $@ obj/GenEvent.o $(OBJS) $(LDLIBS)

