//To fix:PDE costs more time than sum of the time the four DEs cost.
//#define DEBUG
/*Problem:
  1.what's the strategy when new generated X is out of range.
	Current:is reinitialized among range of X.

  2.within a generation, should a new generated X be updated immediately,
	  so that every following x is generated basing on partial updated population.
	Current:every x is generated basing on last generation population, not on current partial updated population.

  3.What to do when some population do not have enough to be migrated to other population?
	Current:I set MinPopSize=5, and when the size of some population is <= MinPopSize,
		the population will not be migrated to other population.

  4.In PDE, what's value of MaxFEs for every single DE exectuted by PDE?
		Current:MaxFEs(DE)=MaxFEs(PDE)/numberOf(DE)

	--by Eric.  Thu Sep 17 20:14:26 CST 2015
	*/
	
#include<stdio.h>
#include<vector>
#include<algorithm>
#include<unistd.h>
#include<stdlib.h>
#include<iostream>
#include<sys/time.h>
#include<time.h>
#include<string.h>
#include<math.h>
void cec14_test_func(double *, double *,int,int,int);
using namespace std;
//set when compile the code.
//#define ALGORITHM 4

#include "include/template.h"

#if ALGORITHM==4
#define OMPI_IMPORTS
#include "mpi.h"
#endif

#define MAX_BUFFER 100
#define MATH_PI M_PI
#define MATH_EXP M_E
#define Trace(m) {cout<<#m"="<<(m)<<endl;}
#define ASSERT(cond) if(!(cond)){cerr<<"Error: condition("#cond") fails!!"<<endl;};
#define Test(m) cout<<#m"={"; m; cout<<"}"<<endl;


class DE{
	private:
		//about function:f
		double xmin;
		double xmax;
		Function *f;
		int numDim;
		bool isFindMin;
		//
		int numP;
		vector<vector<double> >x;//x,trail x.
		vector<vector<double> >tmpX;
		vector<double>fx;
		vector<double>tmpFx;
		vector<double>tx;
		//
		int bestI;
		//
		int algorithm;
		//
		double F,CR;//algorithm related parameters.
		//
		char algorithmName[100];
	public:
		const char *getAlgorithmName()const{
			/*
			switch(algorithm){
				case 0:return "DE1/best/0.1";
				case 1:return "DE2/best/0.9";
				case 2:return "DE3/rand/0.1";
				case 3:return "DE4/rand/0.9";
				case 4:return "PDE/DE1,DE2,DE3,DE4/AMS";
			}
			*/
			return algorithmName;
		}
		inline bool isFBetter(double fx1,double fx2){
			return f->isFBetter(fx1,fx2);
		}
		void init(int algorithm,int numP){
			//Notice:use numP FEs
			this->algorithm=algorithm;
			this->numP=numP;
			F=0.5;
			switch(algorithm){
				case 0:
				case 2:
					CR=0.1;
					break;
				case 1:
				case 3:
					CR=0.9;
					break;
			}
			if(algorithm<4 &&algorithm>=0){
				sprintf(algorithmName,"DE%d/%s/%f/%f",algorithm,(algorithm==0||algorithm==1)?"best":"rand",CR,F);
			}else if(algorithm==4){
				sprintf(algorithmName,"PDE/DE1/DE2/DE3/DE4");
			}else{
				sprintf(algorithmName,"Unknow DE");
			}
		}
		void begin(Function* f){
			this->f=f;
			numDim=f->getNumDim();
			isFindMin=f->getIsFindMin();
		xmin=f->getRange(0);
		xmax=f->getRange(1);
			popInit();
		}
		void update(int maxGeneration){
			for(int g=1;g<=maxGeneration;g++){
				updateX();
			}
		}
		void solve(Function* f,int maxGeneration,vector<double>&bestX,double &bestF){
			begin(f);
			update(maxGeneration);
			getOutput(bestX,bestF);
		}

		void popInit(){
			tx.resize(numDim);
			//allocate space.
			x.resize(numP);
			//
			fx.resize(numP);
			for(int i=0;i<numP;i++){
				x[i].resize(numDim);
			}
			bestI=0;
			for(int i=0;i<numP;i++){
				for(int d=0;d<numDim;d++){
					x[i][d]=drand(xmin,xmax);
				} 
				fx[i]=(*f)(&x[i][0],x[i].size());
				if(isFBetter(fx[i],fx[bestI])){
					bestI=i;
				}
			}
		}
		//not used.
		void calBestI(){
			bestI=0;
			for(int i=0;i<numP;i++){
				if(isFBetter(fx[i],fx[bestI])){
					bestI=i;
				}
			}
		}
		void updateX(){
//			cin.get();
			ASSERT(numP>=3);
			RandomPermutation perm(numP);
			//for next generation:
//			tmpX=x;
//			tmpFx=fx;
//			int tmpBestI=bestI;
//			cout<<endl;
//			cout<<"next generation:"<<endl;
			for(int i=0;i<numP;i++){
				perm.generate();
				int a; int b=perm.next(); int c=perm.next();
				if(algorithm==0||algorithm==1){
					//DE0,DE1
					a=bestI;
				}else{
					a=perm.next();
				}
				int randDim=rand()%numDim;
				for(int j=0;j<numDim;j++){
					if(j==randDim||drand()<CR){
						tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
						if(tx[j]<xmin || tx[j]>xmax){
							tx[j]=drand(xmin,xmax);
						}
					}else{
						tx[j]=x[i][j];
					}
				}
				double ftx=(*f)(&tx[0],tx.size());
				if(isFBetter(ftx,fx[i])){
//	tmpX[i]=tx;
//	tmpFx[i]=ftx;
					x[i]=tx;
					fx[i]=ftx;
//					cout<<"Update i:"<<i<<endl;
//					cout<<"pop["<<i<<"]:";printVec(x[i]);cout<<endl;
//					if(ftx>tmpFx[tmpBestI]){
					if(isFBetter(ftx,fx[bestI])){
//						tmpBestI=i;
						bestI=i;
//						Trace(bestI);
					}
				}
			}
//			ASSERT(fx[bestI]>=tmpFx[tmpBestI]);
//			x=tmpX;
//			fx=tmpFx;
//			cout<<"fx:";printVec(fx);
//			cout<<endl;
//			bestI=tmpBestI;
//			Trace(fx[bestI]);
		}
		vector<double>del(int i){
			ASSERT(i>=0 && i<numP);
			vector<double>XF;
			XF.resize(numDim+1);
			copy(x[i].begin(),x[i].end(),XF.begin());
			XF[numDim]=fx[i];
			x.erase(x.begin()+i);
			fx.erase(fx.begin()+i);
			//
			numP--;
			if(bestI==i){
				calBestI();
			}
			ASSERT(numP>1);
			ASSERT(x.size()==numP&&fx.size()==numP);
			return XF;
		}
		void add(const vector<double>&XF){
			ASSERT(XF.size()==numDim+1);
			//x.resize(numP+1);
			//fx.resize(numP+1);
			x.push_back(XF);
			fx.push_back(XF.back());
			x[numP].pop_back();
			//
			if(isFBetter(fx[numP],fx[bestI])){
				bestI=numP;
			}
			numP++;
			ASSERT(x.size()==numP&&fx.size()==numP);
		}
		int getNumP()const{return numP;}
		void getOutput(vector<double>&bestX,double &bestF){
			bestX=x[bestI];
			bestF=fx[bestI];
		}
		//double meanF;
		double getMeanF()const{
			double avg=0;
			for(int i=0;i<numP;i++){
				avg+=fx[i];
			}
			ASSERT(numP>0);
			avg/=(double)numP;
			return avg;
		}
};

/*
*/

#if ALGORITHM==4
int PDE(int processId,int numProcess,Function*f,vector<double>&bestX,double &bestF){
	srand(time(NULL));//srand in every process.
	//MPI:
	const int TAG=99;
	const int TAG2=98;
	MPI_Status status;
	//DE:
	const int MaxFE=300000;
	const int NumAlgorithm=min(numProcess-1,4);
	int numP=50;
	int numDim=f->getNumDim();
	//
	DE de;
	if(processId!=0){
		de.init(processId-1,numP);
		de.begin(f);
	}
	ASSERT(NumAlgorithm>0);
	const int MaxGeneration=MaxFE/(NumAlgorithm*numP);
	//const int MaxGeneration=5000;
	ASSERT(MaxGeneration>0);
	//master
	vector<double>meanFs;
	vector<vector<int> > migrateMap;
	//slave:
	vector<int>migrateVec;
	vector<double>XF;
	vector<double>bestXF;
	XF.resize(numDim+1);

	vector<int>popSizes;
	if(processId==0){
		//for master process
		meanFs.resize(NumAlgorithm);
		migrateMap.resize(NumAlgorithm);
		popSizes.resize(NumAlgorithm);
		for(int i=0;i<NumAlgorithm;i++){
			migrateMap[i].resize(NumAlgorithm);
			popSizes[i]=numP;
		}
	}else{
		migrateVec.resize(NumAlgorithm);
	}
	const int MinPopSize=5;
	for(int g=1;g<=MaxGeneration;g++){
		if(processId==0){
			//master process:
#ifdef DEBUG
			if(g==1)cout<<"Master:Generation("<<g<<") begin"<<endl;
#endif
			double Pmigrate=0.01+0.99*(exp((double)10.0*g/(double)MaxGeneration)-1.0)/(exp(10.0)-1.0);
			for(int i=1;i<=NumAlgorithm;i++){
				MPI_Recv(&meanFs[i-1],1,MPI_DOUBLE,i,TAG,MPI_COMM_WORLD,&status);
			}
			//cout<<"Master:end recv meanFs"<<endl;
			//build migrateMap.
			for(int i=0;i<NumAlgorithm;i++){
				migrateMap[i][i]=0;
				for(int j=i+1;j<NumAlgorithm;j++){
					if(drand()<Pmigrate){
						if(f->isFBetter(meanFs[i],meanFs[j])){
							//i is better than j.
							if(popSizes[j]<=MinPopSize){
								migrateMap[i][j]=0;//send to j.
							}else{
								migrateMap[i][j]=-1;//recve from j;
								popSizes[i]++;
								popSizes[j]--;
							}
						}else{
							//j is better than i
							if(popSizes[i]<=MinPopSize){
								migrateMap[i][j]=0;//send to j.
							}else{
								migrateMap[i][j]=1;//send to j.
								popSizes[i]--;
								popSizes[j]++;
							}
						}
					}else{
							migrateMap[i][j]=0;//do nothing.
					}
					migrateMap[j][i]=-migrateMap[i][j];
				}
			}
#ifdef DEBUG
			cout<<"MigrateMap:"<<endl;
			for(int i=0;i<NumAlgorithm;i++){
				for(int j=0;j<NumAlgorithm;j++){
					cout<<migrateMap[i][j]<<",";
				}
				cout<<endl;
			}
			cout<<endl;
			if(g!=MaxGeneration)cout<<"Master:Generation("<<g+1<<") begin"<<endl;
#endif
			for(int i=1;i<=NumAlgorithm;i++){
				MPI_Send(&migrateMap[i-1][0],migrateMap[0].size(),MPI_INT,i,TAG,MPI_COMM_WORLD);
			}
		}else{
			//slave processes:
			de.update(1);
			//cout<<"Process("<<processId<<"):"<<"end de.update()"<<endl;
			double meanF=de.getMeanF();
#ifdef DEBUG
			cout<<"Process("<<processId<<"):"<<"meanF:"<<meanF<<endl;
#endif
			MPI_Send(&meanF,1,MPI_DOUBLE,0,TAG,MPI_COMM_WORLD);
		
			MPI_Recv(&migrateVec[0],migrateVec.size(),MPI_INT,0,TAG,MPI_COMM_WORLD,&status);
			//cout<<"Process("<<processId<<"):"<<"end recv migrateMap"<<endl;
			//Attention: without explicit synchronization, there might be a deadlock or error.
			for(int i=0;i<NumAlgorithm;i++){
				if(migrateVec[i]==1){
					//1:send
					ASSERT(de.getNumP()>1);
					int randI=rand()%de.getNumP();
					XF=de.del(randI);
#ifdef DEBUG
			cout<<"Process("<<processId<<"):"<<"sending to P("<<i+1<<")"<<endl;
#endif
					MPI_Send(&XF[0],XF.size(),MPI_DOUBLE,i+1,TAG,MPI_COMM_WORLD);
				}else if(migrateVec[i]==-1){
					//-1:recv
					MPI_Recv(&XF[0],XF.size(),MPI_DOUBLE,i+1,TAG,MPI_COMM_WORLD,&status);
					de.add(XF);
				}
			}
		}
	}
	//return only in  process0.
	if(processId==0){
		for(int i=1;i<=NumAlgorithm;i++){
			MPI_Recv(&XF[0],XF.size(),MPI_DOUBLE,i,TAG2,MPI_COMM_WORLD,&status);
			if(i==1){
				bestXF=XF;
			}else{
				if(f->isFBetter(XF.back(),bestXF.back())){
					bestXF=XF;
				}
			}
		}
		bestF=bestXF.back();
		bestX=bestXF;
		bestX.pop_back();
		return 0;
	}else{
		de.getOutput(XF,bestF);
		XF.push_back(bestF);
		MPI_Send(&XF[0],XF.size(),MPI_DOUBLE,0,TAG2,MPI_COMM_WORLD);
	}
	return -1;
}
vector<double> runPDE(int id,int idSize,Function*f,int maxRun){
	vector<double>results;
	results.resize(maxRun);
	for(int run=0;run<maxRun;run++){
		vector<double>bestX;
		double bestF;
		PDE(id,idSize,f,bestX,bestF);
		results[run]=fabs(bestF-(f->getBest()));
	}
	return results;
}
#endif
vector<double> runSerialDE(DE &de,Function*f,int maxRun){
	vector<double>results;
	const int MaxFE=300000;
	int numDim=f->getNumDim();
	vector<double>bestX;
	double bestF;
	results.resize(maxRun);
	for(int run=0;run<maxRun;run++){
		de.solve(f,MaxFE/de.getNumP(),bestX,bestF);
		results[run]=fabs(bestF-(f->getBest()));
	}
	return results;
}


//int old_main(int argc,char *argv[]){
int main(int argc,char *argv[]){
	srand(time(NULL));
	const int maxRun=30;
	const int numDim=30;
	const int numP=50;
	cec14_test_func();
	FunctionFactory &funGenerator=FunctionFactory::Instance(numDim);
	const int numTestFunction=funGenerator.getNumFunction();
#if ALGORITHM==4 
	int id,idSize;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	MPI_Comm_size(MPI_COMM_WORLD,&idSize);
	if(id==0){
		cout<<"Runing PDE "<<maxRun<<"times."<<endl;
		printf("F\tmean\tstd\n");
		Tic::tic("begin");
	}
	for(int i=0;i<numTestFunction;i++){
		Function*f=funGenerator.getFunction(i);
		vector<double>results=runPDE(id,idSize,f,maxRun);
		if(id==0){
			double min,max,mean,std;
			calStatistics(results,min,max,mean,std);
			printf("%s\t%g\t%g\n",f->getShortName(),mean,std);
		}
	}
	if(id==0){
		Tic::tic("end");
		cout<<"end."<<endl;
		cout<<endl;
	}
	MPI_Finalize();
#else
	DE de;
	de.init(ALGORITHM,numP);
	cout<<"Runing DE("<<de.getAlgorithmName()<<") "<<maxRun<<"times."<<endl;
	printf("F\tmean\tstd\n");
	Tic::tic("begin");
	for(int i=0;i<numTestFunction;i++){
		//if(i==3){
		Function*f=funGenerator.getFunction(i);
		vector<double>results=runSerialDE(de,f,maxRun);
		double min,max,mean,std;
		calStatistics(results,min,max,mean,std);
		printf("%s\t%g\t%g\n",f->getShortName(),mean,std);
//		printf("%s\t%g\t%g\n",f->getName(),mean,std);
		//}
	}
	Tic::tic("end");
	cout<<"end."<<endl;
#endif
	return 0;
}
void test(){}
int unused_main(int argc,char *argv[]){
//int main(int argc,char *argv[]){
	srand(time(NULL));
	const int MaxFE=300000;
	int numDim=30;
	const int numP=50;
	int algorithm=0;
	FunctionFactory &funGenerator=FunctionFactory::Instance(numDim);
	const int numTestFunction=funGenerator.getNumFunction();
	DE de;
	de.init(algorithm,numP);
	Tic::tic("begin");
	cout<<"Runing DE("<<de.getAlgorithmName()<<") MaxGeneration:"<<MaxFE/numP<<"in serial way."<<endl;
	for(int i=0;i<numTestFunction;i++){
		if(i==0){
		vector<double>bestX;
		double bestF;
		Function*f=funGenerator.getFunction(i);
		de.solve(f,MaxFE/numP,bestX,bestF);
//		de.solve(f,1000,bestX,bestF);
		cout<<f->getName()<<" bestF:"<<bestF;
		cout<<" bestX:";printVec(bestX);
		cout<<endl;
		cout<<endl;
			}
	}
	Tic::tic("end");
	test();
	return 0;
}
//0:LPSO 1:GPSO 2:BPSO 3:CLPSO
//GPSO seems useless.
