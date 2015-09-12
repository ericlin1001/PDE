outputPath=../output
all:runPPSO runLPSO  runGPSO runBPSO runCLPSO  output
runPPSO:PPSO
	mpirun -np 5 ./PPSO > ppso.txt

runLPSO:LPSO
	./LPSO > lpso.txt

runGPSO:GPSO
	./GPSO > gpso.txt

runBPSO:BPSO
	./BPSO > bpso.txt

runCLPSO:CLPSO
	./CLPSO > clpso.txt

output:lpso.txt
	cat lpso.txt gpso.txt bpso.txt clpso.txt ppso.txt > $(outputPath)/all.txt

clean:
	rm *.txt
