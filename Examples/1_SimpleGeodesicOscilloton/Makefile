appname := main


SRC_DIRS ?= ../../Source .
INC_DIRS := $(shell find $(SRC_DIRS) -type d)
OBJ_DIRS := $(filter-out $(SRC_DIRS),$(INC_DIRS))
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CXX := g++
CXXFLAGS := $(INC_FLAGS) -O3 -std=c++14
LDFLAGS := -lgsl -lgslcblas

srcfiles := $(shell find . -name "*.cpp") $(shell find $(OBJ_DIRS) -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)

dist-clean: clean
	rm -f *~ .depend

include .depend
