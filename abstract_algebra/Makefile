CXX = g++ 
CXXFLAGS = -g -Wall -std=c++2a
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:src/%.cpp=build/%.o)
EXECUTABLE = main
export INCDIR = -I./include/
export LIBDIR =                                                                               
export LIBS =
      
all: $(EXECUTABLE)
	@#./$(EXECUTABLE)
      
$(EXECUTABLE): $(OBJECTS) | build
	@$(CXX) $(CXXFLAGS) $(INCDIR) $(LIBDIR) $^ $(LIBS) -o $@
      
build/%.o: src/%.cpp | build
	@$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@
      
.PHONY: clean
clean:
	@rm -rf build
	@rm -f $(EXECUTABLE)
      
build:
	@mkdir -p $@
