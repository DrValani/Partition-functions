CC=gcc
CXX=g++
RM=rm -fr
CPPFLAGS=-g 
LDFLAGS=-g -fprofile-arcs -ftest-coverage
LDLIBS=-lgmpxx -lgmp

SRC_DIR=lib
FILES=Energy.cc

OBJ_FILES=$(subst .cc,.o,$(FILES))

OUTFILE=main.out



all: $(OBJ_FILES)
	$(CXX) $(LDFLAGS) Debug/$(OBJ_FILES) -o Debug/$(OUTFILE) $(LDLIBS) 

$(OBJ_FILES) : 
	$(CXX) $(LDFLAGS) -c src/$(FILES) -o Debug/$@

clean:
	$(RM) Debug/$(OBJ_FILES) Debug/$(OUTFILE)
