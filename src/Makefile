CC=g++
CFLAGS=-I.

debug: tp3.cpp
	$(CXX) -std=c++11 -g -Wall -Wextra -Wuninitialized -pedantic -o bin/tp3-debug tp3.cpp $(CFLAGS)

tp3: tp3.cpp
	$(CXX) -std=c++11 -O3 -o bin/tp3 tp3.cpp $(CFLAGS)
