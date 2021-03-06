# /********************************************************************************************
# * File:		Makefile
# * Author:		$LastChangedBy: matthew $
# * Revision:	$Revision: 233 $
# * Last Updated:	$LastChangedDate: 2006-11-10 15:03:28 -0500 (Fri, 10 Nov 2006) $
# ********************************************************************************************/

PR_TARGET=PoissonRecon
ST_TARGET=SurfaceTrimmer
PR_SOURCE=CmdLineParser.cpp DumpOutput.cpp Factor.cpp Geometry.cpp MarchingCubes.cpp PlyFile.cpp Time.cpp PoissonRecon.cpp
ST_SOURCE=CmdLineParser.cpp DumpOutput.cpp Factor.cpp Geometry.cpp MarchingCubes.cpp PlyFile.cpp Time.cpp SurfaceTrimmer.cpp

CFLAGS += -fopenmp -std=c++03 -Wall -Wextra -Werror
LFLAGS += -lgomp

CFLAGS_DEBUG = -DDEBUG -g3 -O0
LFLAGS_DEBUG = -O0

CFLAGS_ADDRESS_SANITIZER = $(CFLAGS_DEBUG) -fsanitize=address -fno-omit-frame-pointer
LFLAGS_ADDRESS_SANITIZER = $(LFLAGS_DEBUG) -fsanitize=address

CFLAGS_THREAD_SANITIZER = $(CFLAGS_DEBUG) -fsanitize=thread -fPIC
LFLAGS_THREAD_SANITIZER = $(LFLAGS_DEBUG) -fsanitize=thread -pie

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

address-sanitizer: CFLAGS += $(CFLAGS_ADDRESS_SANITIZER)
address-sanitizer: LFLAGS += $(LFLAGS_ADDRESS_SANITIZER)
address-sanitizer: $(BIN)$(PR_TARGET)
address-sanitizer: $(BIN)$(ST_TARGET)

thread-sanitizer: CFLAGS += $(CFLAGS_THREAD_SANITIZER)
thread-sanitizer: LFLAGS += $(LFLAGS_THREAD_SANITIZER)
thread-sanitizer: $(BIN)$(PR_TARGET)
thread-sanitizer: $(BIN)$(ST_TARGET)

clang-noomp: CFLAGS = -std=c++11 -Wall -Wextra -Werror -DNO_OMP -DCPP11 $(CFLAGS_DEBUG)
clang-noomp: LFLAGS = $(LFLAGS_DEBUG)
clang-noomp: CXX = clang++
clang-noomp: $(BIN)$(PR_TARGET)
clang-noomp: $(BIN)$(ST_TARGET)

release: CFLAGS += $(CFLAGS_RELEASE)
release: LFLAGS += $(LFLAGS_RELEASE)
release: $(BIN)$(PR_TARGET)
release: $(BIN)$(ST_TARGET)

nogradient: CFLAGS += $(CFLAGS_DEBUG) -DNO_GRADIENT_DOMAIN_SOLUTION
nogradient: LFLAGS += $(LFLAGS_DEBUG)
nogradient: $(BIN)$(PR_TARGET)
nogradient: $(BIN)$(ST_TARGET)

noneumann: CFLAGS += $(CFLAGS_DEBUG) -DNO_FORCE_NEUMANN_FIELD
noneumann: LFLAGS += $(LFLAGS_DEBUG)
noneumann: $(BIN)$(PR_TARGET)
noneumann: $(BIN)$(ST_TARGET)

splat1: CFLAGS += $(CFLAGS_DEBUG) -DSPLAT_ORDER_1
splat1: LFLAGS += $(LFLAGS_DEBUG)
splat1: $(BIN)$(PR_TARGET)
splat1: $(BIN)$(ST_TARGET)

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

test: all
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig"
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig"

test-address-sanitizer: address-sanitizer
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig"
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig"

test-thread-sanitizer: thread-sanitizer
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig"
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig"

test-release: release
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig"
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig"

# GRADIENT_DOMAIN_SOLUTION 0
test-nogradient: nogradient
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig"
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig"

# FORCE_NEUMANN_FIELD 0
test-noneumann: noneumann
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig"
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig"

# SPLAT_ORDER 1
test-splat1: splat1
	Test/run-for-dataset.sh "Examples/horse.npts" "horse" "-orig"
	Test/run-for-dataset.sh "Examples/bunny.points.ply" "bunny" "-orig"

test-binary: all
	Test/run-for-dataset.sh "Examples/horse.bnpts" "horse" "-orig"

test-parallel: all
	Test/run-parallel-test.sh "Examples/cube.npts" "cube"
	Test/run-parallel-test.sh "Examples/horse.npts" "horse"

include $(PR_DEPENDS)
include $(ST_DEPENDS)
