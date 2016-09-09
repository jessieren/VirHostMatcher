CC = g++

LSOURCE =  countKmer.cpp computeMeasure.cpp computeMeasure_onlyd2star.cpp

VirHostMatcher: $(LSOURCE)
	$(CC) countKmer.cpp -o countKmer.out
	$(CC) computeMeasure.cpp -o computeMeasure.out
	$(CC) computeMeasure_onlyd2star.cpp -o computeMeasure_onlyd2star.out

clean:
	rm -f *.o countKmer.out computeMeasure.out computeMeasure_onlyd2star.out
