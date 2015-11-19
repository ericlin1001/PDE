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
using namespace std;
//set when compile the code.
//#define ALGORITHM 4

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
void Tagg(const char *str){
	static int count=0;
	count++;
	printf("******Tag%d:%s\n",count,str);
}
double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;
	if ( phase == 0 ) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2.0 * U1 - 1.0;
			V2 = 2.0 * U2 - 1.0;
			S = V1 * V1 + V2 * V2;
		} while(S >= 1 || S == 0);
		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);
	phase = 1 - phase;
	return X;
}
inline double NormD(double u,double t){
	return gaussrand()*t+u;
}
double drand(){
	//[0,1);
	double r=rand();
	r/=((double)RAND_MAX+1);
	return r;
}
double drand(double min,double max){
	return drand()*(max-min)+min;
}

template<class T>
void printVec(const vector<T>&arr){
	cout<<"(";
	for(int i=0;i<arr.size();i++){
		if(i!=0)cout<<',';
		cout<<arr[i];
	}
	cout<<")";
}
class Tic{
	//accuration in milliseconds
	private:
		static long lastTime;
		Tic(){}
		inline static long getTimeMs(){
			timeval timeStart;
			gettimeofday(&timeStart,NULL);
			long res=((long)timeStart.tv_sec)*1000+(long)timeStart.tv_usec/1000;
			return res;
		}
	public:
		static long mtic(){
			//in milliseconds.
			long currentTime=getTimeMs();
			long dur=currentTime-lastTime;
			lastTime=currentTime;
			return dur;
		}
		static void tic(const char *tag="begin"){
			if(strcmp(tag,"begin")==0){
				cout<<"Tic::"<<tag<<endl;
				dtic();
			}else{
				cout<<"Tic::"<<tag<<" used:"<<dtic()<<"(seconds)."<<endl;
			}
		}
		inline static double dtic(){
			//in seconds.
			return (double)mtic()/1000.0;
		}
		static void test(){
			Tic::mtic();
			usleep(1234000);//sleep for 1234 milliseconds.(1.234seconds)
			cout<<Tic::dtic()<<"seconds"<<endl;
		}

};
long Tic::lastTime=0;

void printArr(int *arr,int size){
	cout<<"(";
	for(int i=0;i<size;i++){
		if(i!=0)cout<<',';
		cout<<arr[i];
	}
	cout<<")";
}
class Function{
#define MAX_FUNCTION_NAME 150
	private:
		char shortName[50];
		char funName[MAX_FUNCTION_NAME];
		double xlow,xup;
		double fbest;
		bool isFindMin;
		int numDim;
		//
		int feCounter;
	private:
	public:
		static double u(double x,double a,double k,double m){
			if(x>a)return k*pow(x-a,m);
			if(x<-a)return k*pow(-x-a,m);
			return 0;
		}
		virtual double operator()(const double *xs,int size){
			feCounter++;
			return 0;
		}
		inline double f(const vector<double>&xs){
			return operator()(&xs[0],xs.size());
		}
	public:
		Function(const char *s,double xlow,double xup,double fbest,bool isFindMin,int numDim){
			this->xlow=xlow;
			this->xup=xup;
			this->fbest=fbest;
			this->isFindMin=isFindMin;
			this->numDim=numDim;
			strcpy(shortName,s);
			if(isFindMin){
				sprintf(funName,"{%s(%f,%f)^%d fmin:%f}",s,xlow,xup,numDim,fbest);
			}else{
				sprintf(funName,"{%s(%f,%f)^%d fmax:%f}",s,xlow,xup,numDim,fbest);
			}
			feCounter=0;
		}
		int getfeCounter()const{return feCounter;}
		double getBest()const{return fbest;}
		bool getIsFindMin()const{return isFindMin;}
		inline bool isFBetter(double a,double b){
			return isFindMin^(a>=b);
		}
		int getNumDim()const{return numDim;}
		double getRange(int botOrUp){
			if(botOrUp==0)return xlow;
			return xup;
		}
		const char *getShortName()const{return shortName;}
		const char *getName()const{return funName;}
};

#define DefFunction(name,xlow,xup,fbest,isFindMin) class name : public Function{\
	public: name(int numDim):Function(#name,xlow,xup,fbest,isFindMin,numDim){}\
			virtual double operator()(const double *xs,int size){\
				Function::operator()(xs,size);
#define EndDef }};	
DefFunction(PDEF3,-10,10,0,true)
	double res=0.0;
	for(int i=0;i<size-1;i++){
		double t=xs[i+1]-xs[i];
		double t1=xs[i]-1;
		res+=t*t*100.0+t1*t1;
	}
	return res;
EndDef

DefFunction(PDEF4,-500,500,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		res+=-xs[i]*sin(sqrt(fabs(xs[i])));
	}
	res+=(double)418.9829*(double)size;
	return res;
EndDef

DefFunction(F1,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double x=xs[i];
		res+=x*x;
	}
return res;
	EndDef

DefFunction(F2,-10,10,0,true)
	double res=0.0;
	double sum=0.0;
	double mul=1.0;
	for(int i=0;i<size;i++){
		double fabsx=fabs(xs[i]);
		sum+=fabsx;
		mul*=fabsx;
	}
res=sum+mul;
return res;
EndDef

DefFunction(F3,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double insum=0.0;
		for(int j=0;j<=i;j++){
			insum+=xs[j];
		}
		res+=insum*insum;
	}
return res;
EndDef

DefFunction(F4,-100,100,0,true)
	double res=fabs(xs[0]);
	for(int i=1;i<size;i++){
		double tmp=fabs(xs[i]);
		if(tmp<res)res=tmp;
	}
return res;
EndDef

//untest:
DefFunction(F5,-30,30,0,true)
	double res=0.0;
	for(int i=0;i<size-1;i++){
		double tmp=pow(xs[i+1]-xs[i]*xs[i],2)*100.0+pow(xs[i]-1.0,2);
		res+=tmp;
	}
return res;
EndDef

DefFunction(F6,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		int tmp=floor(xs[i]+0.5);
		res+=tmp*tmp;
	}
return res;
EndDef

DefFunction(F7,-1.28,1.28,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=pow(xs[i],4)*(double)(i+1);
		res+=tmp;
	}
	res+=drand();
return res;
EndDef

DefFunction(F8,-500,500,-12569.5,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=-xs[i]*sin(sqrt(fabs(xs[i])));
		res+=tmp;
	}
return res;
EndDef

DefFunction(F9,-5.12,5.12,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=pow(xs[i],2)-(double)10.0*cos(xs[i]*2.0*MATH_PI)+10.0;
		res+=tmp;
	}
return res;
EndDef

DefFunction(F10,-32,32,0,true)
	double res=0.0;
	double sumx2=0.0;
	double sumcosx=0.0;
	for(int i=0;i<size;i++){
		sumx2+=pow(xs[i],2);
		sumcosx+=cos(xs[i]*MATH_PI*2.0);
	}
res=-20.0*exp(-0.2*sqrt(sumx2/(double)size))-exp(sumcosx/(double)size)+20.0+MATH_EXP;
return res;
EndDef

DefFunction(F11,-600.0,600.0,0,true)
	double res=0.0;
	double sumx2=0.0;
	double mulcos=1.0;
	for(int i=0;i<size;i++){
		sumx2+=pow(xs[i],2);
		mulcos*=cos(xs[i]/sqrt((double)i+1));
	}
res=sumx2/4000.0-mulcos+1.0;
return res;
EndDef

DefFunction(F12,-50,50,0,true)
	double res=0.0;
	double y1=1.0+(xs[0]+1.0)/4.0;
	double yd=1.0+(xs[size-1]+1.0)/4.0;
	double sumy=0.0;
	double sumu=0.0;
	//
	double yi,yi1;
	yi=y1;
	for(int i=0;i<size-1;i++){
		yi1=1.0+(xs[i+1]+1.0)/4.0;
		sumy+=pow(yi-1.0,2)*(1.0+10.0*pow(sin(MATH_PI*yi1),2));
		yi=yi1;
	}
for(int i=0;i<size;i++){
	sumu+=Function::u(xs[i],10,100,4);
}
res=MATH_PI/(double)size*(10.0*pow(sin(MATH_PI*y1),2)+sumy+pow(yd-1,2))
	+sumu;
	return res;
	EndDef

DefFunction(F13,-50,50,0,true)
	double res=0.0;
	double sumx=0.0;
	double sumu=0.0;
	for(int i=0;i<size-1;i++){
		sumx+=pow(xs[i]-1,2)*(1+pow(sin(3.0*MATH_PI*xs[i+1]),2));
	}
for(int i=0;i<size;i++){
	sumu+=Function::u(xs[i],5,100,4);
}
double xd=xs[size-1];
res=0.1*(pow(sin(3.0*MATH_PI*xs[0]),2)+sumx+
		pow(xd-1.0,2)*(1+pow(sin(2.0*MATH_PI*xd),2)))+sumu;
return res;
EndDef

class RandomPermutation{
	vector<int>p;
	int i;
	int n;
	public:
		RandomPermutation(int tn):n(tn),i(0){
			p.resize(n);
			for(int i=0;i<n;i++){
				p[i]=i;
			}
		}
		void generate(){
			i=0;
		}
		int next(){
			//at most invoked n times before next generate.
			ASSERT(i<n&&i>=0);
			int r=rand()%(n-i);
			i++;
			swap(p[n-i],p[r]);
			return p[n-i];
		}
};
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

class FunctionFactory{
	private:
		vector<Function*>fs;
		FunctionFactory(int numDim){
			fs.resize(4);
			fs[0]=new F1(numDim);
			fs[1]=new F3(numDim);
			fs[2]=new PDEF3(numDim);
			fs[3]=new PDEF4(numDim);
		}
		static FunctionFactory*instance;
	public:
		static FunctionFactory &Instance(int numDim){
			if(instance==0)instance=new FunctionFactory(numDim);
			return *instance;
		}
		/*
		   void setNumDim(int numDim){
		   }
		 */
		Function*getFunction(int index)const{
			return fs[index];
		}
		int getNumFunction()const{
			return fs.size();
		}
		~FunctionFactory(){
			for(int i=0;i<getNumFunction();i++){
				delete fs[i];
			}
		}
};
FunctionFactory*FunctionFactory::instance=0;

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
	int numP=50;
	const int MaxFE=300000;
	const int NumAlgorithm=min(numProcess-1,4);
	int numDim=f->getNumDim();
	//
	DE de;
	if(processId!=0){
		de.init(processId-1,numP);
		de.begin(f);
	}
	ASSERT(NumAlgorithm>0);
	const int MaxGeneration=MaxFE/(numP*NumAlgorithm);
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
void calStatistics(const vector<double>&arr,double &min,double &max,double &mean,double &std){
	min=arr[0];
	max=arr[0];
	mean=0.0;
	for(int i=0;i<arr.size();i++){
		double x=arr[i];
		if(x<min)min=x;
		if(x>max)max=x;
		mean+=x;
	}
	mean/=(double)arr.size();
	std=0.0;
	for(int i=0;i<arr.size();i++){
		double x=arr[i];
		std+=pow(x-mean,2);
	}
	std/=(double)arr.size();
	std=sqrt(std);
}

//int old_main(int argc,char *argv[]){
int main(int argc,char *argv[]){
	srand(time(NULL));
	const int maxRun=30;
	const int numDim=30;
	const int numP=50;
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
