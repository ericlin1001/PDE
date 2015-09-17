outputPath=../output
all:runDEs output
runDEs:PDE
	echo 'Run Start:'`date` > time.txt
	mpirun -np 5 ./PDE > pde.txt
	./DE1 > de1.txt
	./DE2 > de2.txt
	./DE3 > de3.txt
	./DE4 > de4.txt
	echo 'Run End:'`date` >> time.txt

output:pde.txt
	cat time.txt de1.txt de2.txt de3.txt de4.txt pde.txt > $(outputPath)/all.txt

clean:
	rm -f *.txt
