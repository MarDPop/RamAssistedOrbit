# *****************************************************
# Variables to control Makefile operation
 
CC = g++ -std=c++17
CFLAGS = -Wall -g3

SOURCES = $(wildcard src/*.cpp)
OBJECTS=$(patsubst src/%.cpp, obj/%.o, $(SOURCES))

EXECUTABLE = bin/sim.exe

all:    build $(EXECUTABLE)

$(EXECUTABLE):  $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) 

$(OBJECTS): obj/%.o : src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

build:
	@mkdir -p bin

clean:
	rm -rf $(EXECUTABLE) $(OBJECTS) bin