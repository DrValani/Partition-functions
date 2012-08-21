CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-g 
LDFLAGS=-g 
LDLIBS=-lgmpxx -lgmp

SRCS=Energy.cc
OBJS=$(subst .cc,.o,$(SRCS))

all: Energy
	
Energy: $(OBJS)
	g++ $(LDFLAGS) -o Energy $(OBJS) $(LDLIBS) 


clean:
	$(RM) $(OBJS)
