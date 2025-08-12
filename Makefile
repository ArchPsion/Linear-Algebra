CXX=g++ -std=c++23
CXXFLAGS=-W -Wall -Wextra -pedantic -O3
EXEC=foo
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)
INC= $(wildcard *.hpp)

all: $(EXEC)

foo: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

foo.o: $(INC)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

.PHONY: clean mrproper

clean:
	@rm -rf *.o

mrproper: clean
	@rm -rf $(EXEC)
