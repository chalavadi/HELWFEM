CC=g++
CFLAGS=-c -Wall
ARMA=-O3 -larmadillo
OBJS=testElement.o element.o
LIBGT=libgtest.a
INCL=-isystem gtest/include

all:
	make testElement
	
testElement: $(OBJS)
	$(CC) $(INCL) -pthread $(OBJS) $(LIBGT) -o testElement $(ARMA)

testElement.o: testElement.cpp ../src/element.cpp ../src/element.h test.h
	$(CC) $(CFLAGS) $(INCL) testElement.cpp -o testElement.o $(ARMA)

element.o: ../src/element.cpp
	$(CC) $(CFLAGS) ../src/element.cpp -o element.o $(ARMA)

clean:
	rm -rf $(OBJS)

