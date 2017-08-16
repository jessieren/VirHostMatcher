CC = g++

LSOURCE =  countKmer.cpp computeMeasure.cpp computeMeasure_onlyd2star.cpp

VirHostMatcher: $(LSOURCE)
	$(CC) countKmer.cpp -o countKmer.out -std=c++11
	$(CC) computeMeasure.cpp -o computeMeasure.out -std=c++11
	$(CC) computeMeasure_onlyd2star.cpp -o computeMeasure_onlyd2star.out -std=c++11

clean:
	rm -f *.o countKmer.out computeMeasure.out computeMeasure_onlyd2star.out
