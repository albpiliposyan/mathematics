CXX = clang++
CXXFLAGS = -g -Wall -std=c++2a
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:src/%.cpp=build/%.o)
EXECUTABLE = main
export INCDIR = -I./include/
export LIBDIR =
export LIBS =

all: clean $(EXECUTABLE)
	@#./$(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) | build
	@$(CXX) $(CXXFLAGS) $(INCDIR) $(LIBDIR) $^ $(LIBS) -o $@

build/%.o: src/%.cpp | build
	@$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@


.PHONY: clean
clean:
	@if test -d build; then rm -rf build; fi
	@if test $(EXECUTABLE); then rm -rf $(EXECUTABLE); fi

build:
	@mkdir -p $@
