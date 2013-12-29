CC=g++
CFLAGS=-c -Wall
ARMA=-O3 -larmadillo
OBJS=main.o element.o

all: $(OBJS)
	$(CC) $(OBJS) -o fem1 $(ARMA)

main.o: main.cpp element.cpp element.h
	$(CC) $(CFLAGS) main.cpp -o main.o $(ARMA)

element.o: element.cpp
	$(CC) $(CFLAGS) element.cpp -o element.o $(ARMA)

clean:
	rm -rf *.o fem1
