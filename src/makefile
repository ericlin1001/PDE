all:realRun

test:pde.cpp
	mpic++ -g -DALGORITHM=1 -o test pde.cpp
	./test

realRun:pde.cpp
	cd .. && make
