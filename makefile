binPath=./bin
cpp=./src/pde.cpp
all:compile run

compile:$(cpp)
	mpic++ -O2 -o $(binPath)/DE1 -DALGORITHM=0 $(cpp)
	mpic++ -O2 -o $(binPath)/DE2 -DALGORITHM=1 $(cpp)
	mpic++ -O2 -o $(binPath)/DE3 -DALGORITHM=2 $(cpp)
	mpic++ -O2 -o $(binPath)/DE4 -DALGORITHM=3 $(cpp)
	mpic++ -O2 -o $(binPath)/PDE -DALGORITHM=4 $(cpp)

run:$(binPath)/LPSO
	cd $(binPath) && make
	@echo ****Result is saved in ./output/all.txt******



