CC=gcc
CXX=g++
RM=rm -fr
CPPFLAGS=-g 
LDFLAGS=-g -fprofile-arcs -ftest-coverage
LDLIBS=-lgmpxx -lgmp

SRC_DIR=src
FILES=Polynomial.cc Energy.cc
SRC_FILES=$(addprefix $(FILES)/,$(FILES)
OBJ_DIR=Debug
OBJ_FILES=$(subst .cc,.o,$(FILES))
OBJ_OUT_FILES=$(addprefix $(OBJ_DIR)/,$(OBJ_FILES))

OUTFILE=$(OBJ_DIR)/main.out



all: $(OBJ_FILES)
	$(CXX) $(LDFLAGS) $(OBJ_OUT_FILES) -o $(OUTFILE) $(LDLIBS) 

Polynomial.o : 
	$(CXX) $(LDFLAGS) -c src/Polynomial.cc -o Debug/$@

Energy.o : 
	$(CXX) $(LDFLAGS) -c src/Energy.cc -o Debug/$@


clean:
	$(RM) $(OBJ_OUT_FILES) $(OUTFILE)
