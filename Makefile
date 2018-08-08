# Filename: Makefile
# Author: Scrubbs
# Date: 2018-8-7
# Description: molopt Makefile

#Project Name (executable)
PROJECT = molopt

#Compiler
CXX = g++

#Compiler Directions
CXXFLAGS=

#Source Files
SOURCE = $(wildcard *.c *.cpp)

#Object Files
OBJECTS = $(SOURCE:.cpp=.o)

#Libraries
LIBDIR=lib
LIBRARIES = $(wildcard .lib .a $(LIBDIR)/*.a $(LIBDIR)/*.lib)

#Build Executable
$(PROJECT): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBRARIES) 

#Clean up additional files
.PHONY: clean
clean:
	rm -f *.o *.exe