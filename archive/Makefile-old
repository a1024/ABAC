# -*- Makefile -*-
PROGRAM = abac
SOURCES = $(wildcard src/*.cpp)
CXXFLAGS =

debug:	CXXFLAGS+=-g -D_DEBUG
debug:	$(PROGRAM)

release:	CXXFLAGS+=-O
release:	$(PROGRAM)

$(PROGRAM):	$(SOURCES)
	$(CXX) -no-pie $(CXXFLAGS) -march=core2 $(SOURCES) -o $(PROGRAM)

run:
	./$(PROGRAM) "example.png"

clean:
	rm $(PROGRAM)
