CC = g++

CFLAGS = -std=c++0x -Wall

.PHONY: clean

all: countKmer computeMeasure computeMeasure_onlyd2star

countKmer:
	$(CC) $@.cpp -o $@.out $(CFLAGS)

computeMeasure:
	$(CC) $@.cpp -o $@.out $(CFLAGS)

computeMeasure_onlyd2star:
	$(CC) $@.cpp -o $@.out $(CFLAGS)

clean:
	rm -f *.o countKmer.out computeMeasure.out computeMeasure_onlyd2star.out
