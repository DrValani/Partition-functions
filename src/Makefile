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
	$(CXX) $(LDFLAGS) $(OBJ_FILES) -o $(OUTFILE) $(LDLIBS) 

$(OBJ_FILES) : 
	$(CXX) $(LDFLAGS) -c $(FILES) -o $@

clean:
	$(RM) $(OBJ_FILES) $(OUTFILE)
