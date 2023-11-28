CXX = clang++
CXXFLAGS = -g -Wall -std=c++2a
SOURCES = $(wildcard src/*.cpp)
TEST_FILE = test.cpp
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


test: clean
	@$(CXX) $(CXXFLAGS) $(TEST_FILE) -o test
	@./test

.PHONY: clean
clean:
	@rm -rf build
	@rm -f $(EXECUTABLE)
	@rm -f test

build:
	@mkdir -p $@
