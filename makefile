binPath=./bin
cpp=./src/pde.cpp
all:compile run

alltest:$(cpp)
	mpic++ -O2 -g -o $(binPath)/PDE -DALGORITHM=4 $(cpp)
	mpic++ -O2 -g -o $(binPath)/DE1 -DALGORITHM=0 $(cpp)
#	mpirun -n 5 xterm -e gdb $(binPath)/PDE 
	mpirun -n 5 $(binPath)/PDE 
	$(binPath)/DE1

	

compile:$(cpp)
	mpic++ -O2 -o $(binPath)/DE1 -DALGORITHM=0 $(cpp)
	mpic++ -O2 -o $(binPath)/DE2 -DALGORITHM=1 $(cpp)
	mpic++ -O2 -o $(binPath)/DE3 -DALGORITHM=2 $(cpp)
	mpic++ -O2 -o $(binPath)/DE4 -DALGORITHM=3 $(cpp)
	mpic++ -O2 -o $(binPath)/PDE -DALGORITHM=4 $(cpp)

run:$(binPath)/PDE
	cd $(binPath) && make
	@echo ****Result is saved in ./output/all.txt******



