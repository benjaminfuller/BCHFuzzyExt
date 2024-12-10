
CPLUS_INCLUDE_PATH = /usr/include/x86_64-linux-gnu/
CXXFLAGS = -Wno-deprecated 

all: sketch

bch.o: bch.cpp bchsketch.h
	g++ $(CXXFLAGS) -c bch.cpp 

io.o: io.cpp bchsketch.h
	g++ $(CXXFLAGS) -c io.cpp 

sketch.o: sketch.cpp bchsketch.h
	g++ $(CXXFLAGS) -c sketch.cpp 


sketch: sketch.o bch.o io.o bchsketch.h
	g++  $(CXXFLAGS) sketch.o io.o bch.o -lntl -lgmp -o sketch


clean:
	rm  sketch bch.o io.o sketch.o 
