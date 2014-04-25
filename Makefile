# /********************************************************************************************
# * File:		Makefile
# * Author:		$LastChangedBy: matthew $
# * Revision:	$Revision: 233 $
# * Last Updated:	$LastChangedDate: 2006-11-10 15:03:28 -0500 (Fri, 10 Nov 2006) $
# ********************************************************************************************/

PR_TARGET=PoissonRecon
ST_TARGET=SurfaceTrimmer
PR_SOURCE=CmdLineParser.cpp Factor.cpp Geometry.cpp MarchingCubes.cpp PlyFile.cpp Time.cpp PoissonRecon.cpp
ST_SOURCE=CmdLineParser.cpp Factor.cpp Geometry.cpp MarchingCubes.cpp PlyFile.cpp Time.cpp SurfaceTrimmer.cpp

CFLAGS += -fopenmp -std=c++11 -Wall -Wextra -Werror
LFLAGS += -lgomp

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math
LFLAGS_RELEASE = -O3

SRC = Src/
BIN = Bin/

CXX=g++

PR_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PR_SOURCE))))
ST_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(ST_SOURCE))))

PR_DEPENDS=$(addprefix $(BIN), $(addsuffix .d, $(basename $(PR_SOURCE))))
ST_DEPENDS=$(addprefix $(BIN), $(addsuffix .d, $(basename $(ST_SOURCE))))

all: CFLAGS += $(CFLAGS_DEBUG)
all: LFLAGS += $(LFLAGS_DEBUG)
all: $(BIN)$(PR_TARGET)
all: $(BIN)$(ST_TARGET)

release: CFLAGS += $(CFLAGS_RELEASE)
release: LFLAGS += $(LFLAGS_RELEASE)
release: $(BIN)$(PR_TARGET)
release: $(BIN)$(ST_TARGET)

clean:
	rm -f $(BIN)$(PR_TARGET)
	rm -f $(BIN)$(ST_TARGET)
	rm -f $(PR_OBJECTS) $(ST_OBJECTS)
	rm -f $(PR_DEPENDS) $(ST_DEPENDS)

$(BIN)$(PR_TARGET): $(PR_OBJECTS)
	$(CXX) -o $@ $(PR_OBJECTS) $(LFLAGS)

$(BIN)$(ST_TARGET): $(ST_OBJECTS)
	$(CXX) -o $@ $(ST_OBJECTS) $(LFLAGS)

$(BIN)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) $<

$(BIN)%.d: $(SRC)%.cpp
	$(CXX) $(CFLAGS) -MF"$@" -MG -MM -MP -MT "$(addprefix $(BIN), $(addsuffix .o, $(notdir $(basename $<))))" "$<"

test:
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig"
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig"

test-release:
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig-release"
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig-release"

include $(PR_DEPENDS)
include $(ST_DEPENDS)
