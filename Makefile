CXX=g++

CFlags=-c -Wall -O3 -pedantic -MMD -std=c++11 -Werror
CFlags+=-ffast-math

Sources=$(wildcard src/*.cpp)
IncludeDir=-I./include -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/
AllObjects=$(addprefix obj/,$(notdir $(Sources:.cpp=.o)))
Executables=main
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p data obj fig
	$(CXX) $(Objects) $(addprefix obj/,$@.o) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFlags) $(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

test: $(Executables)
	$(foreach exe,$(Executables),./$(exe);)

clean:
	rm -rf obj/*.o obj/*.d $(Executables)

