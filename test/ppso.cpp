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
#define ALGORITHM 0

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
/*
double sampleNormal() {
	double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
	double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
	double r = u * u + v * v;
	if (r == 0 || r > 1) return sampleNormal();
	double c = sqrt(-2 * log(r) / r);
	return u * c;
}
*/
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
			ASSERT((2*U1)==(2.0*U1));
			Trace(U1);
			Trace(2*U1);
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
	//return sampleNormal()*t+u;
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
		static inline double u(double x,double a,double k,double m){
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

class PSO{
	private:
		//about function:f
		Function *f;
		int numDim;
		bool isFindMin;
		//
		int numP;
		vector<vector<double> >v,x,pBest,gBest;
		vector<double>fx,fpBest,fgBest;
		int algorithm;
	public:
		const char *getAlgorithmName(){
			switch(algorithm){
				case 0:return "LPSO";
				case 1:return "GPSO";
				case 2:return "BPSO";
				case 3:return "CLPSO";
				case 4:return "PPSO";
			}
			return 0;
		}
		bool isFBetter(double fx1,double fx2){
			if(fx1<fx2){
				return isFindMin;
			}else{
				return !isFindMin;
			}
		}
		/*
		   PSO(int algorithm,int numP){
		   this->algorithm=algorithm;
		   this->numP=numP;
		   }
		 */
		void init(int algorithm,int numP){
			//Notice:use numP FEs
			this->algorithm=algorithm;
			this->numP=numP;
		}
		void begin(Function* f){
			this->f=f;
			numDim=f->getNumDim();
			isFindMin=f->getIsFindMin();
			popInit();
		}
		void replaceParticle(int i,vector<double>&bestX,double &bestF){
			x[i]=bestX;
			fx[i]=bestF;
			//
			if(isFBetter(fx[i],fpBest[i])){
				pBest[i]=x[i];
				fpBest[i]=fx[i];
			}
			updateGBest();
		}
		void update(int maxGeneration){
			for(int g=1;g<=maxGeneration;g++){
				updateX(g,maxGeneration);
				updateGBest();
			}
		}
		void solve(Function* f,int maxGeneration,vector<double>&bestX,double &bestF){
			begin(f);
			update(maxGeneration);
			getOutput(bestX,bestF);
		}

		void popInit(){
			//allocate space.
			v.resize(numP);
			x.resize(numP);
			pBest.resize(numP);
			gBest.resize(numP);
			//
			fx.resize(numP);
			fpBest.resize(numP);
			fgBest.resize(numP);
			for(int i=0;i<numP;i++){
				v[i].resize(numDim);
				x[i].resize(numDim);
				pBest[i].resize(numDim);
				gBest[i].resize(numDim);
			}
			//init v,x,pBest,gBest.
			const double xmin=f->getRange(0);
			const double xmax=f->getRange(1);
			const double vmin=xmin/5.0;
			const double vmax=xmax/5.0;
			for(int i=0;i<numP;i++){
				//init x,v
				for(int d=0;d<numDim;d++){
					x[i][d]=drand(xmin,xmax);
					v[i][d]=drand(vmin,vmax);
				}
				fx[i]=(*f)(&x[i][0],x[i].size());
			}
			pBest=x;
			fpBest=fx;
			updateGBest();
		}
		void updateX(int g,int maxGeneration){
			const double xmin=f->getRange(0);
			const double xmax=f->getRange(1);
			const double vmin=xmin/5.0;
			const double vmax=xmax/5.0;

			const double wmax=0.9,
				  wmin=0.4;
			double w;
			switch(algorithm){
				case 0:
				case 1:
					//LPSO,GPSO
					{
						const double  c1=2.0,
							  c2=2.0;
						if(algorithm==1){
							//GPSO
							w=wmax-(wmax-wmin)*(double)g/(double)maxGeneration;
							double t=(double)g/(double)maxGeneration;
							t-=1;
							w=wmin+(wmax-wmin)*t*t;
						}else{
							//LPSO
							w=wmax-(wmax-wmin)*g/maxGeneration;
						}
						for(int i=0;i<numP;i++){
							for(int d=0;d<numDim;d++){
								//update V,X
								v[i][d]=w*v[i][d]+c1*drand()*(pBest[i][d]-x[i][d])+
									c2*drand()*(gBest[i][d]-x[i][d]);
								x[i][d]=x[i][d]+v[i][d];

								if(x[i][d]<xmin){
									x[i][d]=xmin;
									v[i][d]=drand(vmin,vmax);
								}else if(xmax<x[i][d]){
									x[i][d]=xmax;
									v[i][d]=drand(vmin,vmax);
								}
								if(v[i][d]<vmin){
									v[i][d]=vmin;
								}else if(vmax<v[i][d]){
									v[i][d]=vmax;
								}
							}
							//update Fx
							fx[i]=(*f)(&x[i][0],x[i].size());
							//update PBest.
							if(isFBetter(fx[i],fpBest[i])){
								pBest[i]=x[i];
								fpBest[i]=fx[i];
							}
						}
					}
					break;
				case 2:
					//BPSO
					for(int i=0;i<numP;i++){
						for(int d=0;d<numDim;d++){
							double p_u=0.5*(pBest[i][d]+gBest[i][d]);
							double p_t=fabs(pBest[i][d]-gBest[i][d]);

							x[i][d]=NormD(p_u,p_t);

							//
							//BPSO-redistribution of X
							while(x[i][d]<xmin||x[i][d]>xmax){
								x[i][d]=NormD(p_u,p_t);
							}


						}
						//update Fx
						fx[i]=(*f)(&x[i][0],x[i].size());
						//update PBest.
						if(isFBetter(fx[i],fpBest[i])){
							pBest[i]=x[i];
							fpBest[i]=fx[i];
						}
					}
					break;
				case 3:
					//CLPSO
					{
						const double c1=1.49445;
						w=wmax-(wmax-wmin)*g/maxGeneration;
						for(int i=0;i<numP;i++){
							vector<int>fi=CLSPO_createFi(i);
							bool isUpdatePbest=true;
							for(int d=0;d<numDim;d++){
								v[i][d]=w*v[i][d]+c1*drand()*(pBest[fi[d]][d]-x[i][d]);
								x[i][d]=x[i][d]+v[i][d];

								if(isUpdatePbest){
									if(x[i][d]<xmin||xmax<x[i][d]){
										isUpdatePbest=false;
									}
								}
								if(v[i][d]<vmin){
									v[i][d]=vmin;
								}else if(vmax<v[i][d]){
									v[i][d]=vmax;
								}

								//checkVX(x[i][d],v[i][d]);
							}
							//update Fx
							fx[i]=(*f)(&x[i][0],x[i].size());
							//update PBest.
							if(isUpdatePbest){
								if(isFBetter(fx[i],fpBest[i])){
									pBest[i]=x[i];
									fpBest[i]=fx[i];
								}
							}

						}
					}
					break;

				default:
					cerr<<"Error: PSO.update()"<<endl;
					break;
			}
			//checkXVInRange();
		}
		vector<int>CLSPO_createFi(int i){
			vector<int>fi;
			fi.resize(numDim);
			double Pci=0.05+0.45*(exp(10.0*i/(numP-1))-1)/(exp(10.0)-1);
			bool isAllItSelf=true;
			for(int d=0;d<numDim;d++){
				if(drand()<Pci){
					int i1=i+(int)(drand()*(numP-i));
					int i2=i+(int)(drand()*(numP-i));
					if(isFBetter(fpBest[i1],fpBest[i2])){
						fi[d]=i1;
					}else{
						fi[d]=i2;
					}
					if(isAllItSelf&&fi[d]!=i){
						isAllItSelf=false;
					}
				}else{
					fi[d]=i;
				}
			}
			if(isAllItSelf){
				//rule3
				int randD=drand()*numDim;
				int i1;
				do{
					i1=drand()*numP;
				}while(i==i1);
				fi[randD]=i1;
			}
			return fi;
		}
		//update gBest base on pBest.
		void updateGBest(){
			//This's LPSO
			//TODO:be careful the strategy of update gBest.
			//update gBest,fgBest
			switch(algorithm){
				case 0:
					//LPSO
					for(int i=0;i<numP;i++){
						int l=(i+numP-1)%numP;
						int r=(i+1)%numP;
						int bi=i;
						if(isFBetter(fpBest[l],fpBest[bi])){
							bi=l;
						}
						if(isFBetter(fpBest[r],fpBest[bi])){
							bi=r;
						}
						gBest[i]=pBest[bi];
						fgBest[i]=fpBest[bi];
					}
					break;
				case 1:
				case 2:
					{
						//GPSO
						//BPSO
						int bestI=0;
						for(int i=1;i<numP;i++){
							if(isFBetter(fpBest[i],fpBest[bestI])){
								bestI=i;
							}
						}
						//
						vector<double>&tmpGBest=pBest[bestI];
						double tmpGBestF=fpBest[bestI];
						for(int i=0;i<numP;i++){
							gBest[i]=tmpGBest;
							fgBest[i]=tmpGBestF;
						}
					}
					break;
				case 3:
					//CLPSO
					//GBest is not used.
					break;
				default:
					cerr<<"Error:updateGBest() type is not supported!"<<endl;
					break;
			}
		}
		int getNumP()const{return numP;}
		void getOutput(vector<double>&bestX,double &bestF){
			int bestI;
			switch(algorithm){
				case 0:
					if(isFindMin){
						bestI=distance(fgBest.begin(),min_element(fgBest.begin(),fgBest.end()));
					}else{
						bestI=distance(fgBest.begin(),max_element(fgBest.begin(),fgBest.end()));
					}
					bestX=gBest[bestI];
					bestF=fgBest[bestI];
					break;
				case 1:
				case 2:
					bestX=gBest[0];
					bestF=fgBest[0];
					break;
				case 3:
					if(isFindMin){
						bestI=distance(fpBest.begin(),min_element(fpBest.begin(),fpBest.end()));
					}else{
						bestI=distance(fpBest.begin(),max_element(fpBest.begin(),fpBest.end()));
					}
					bestX=pBest[bestI];
					bestF=fpBest[bestI];
					break;
				default:
					break;
			}
		}
};

class FunctionFactory{
	private:
		//Function*fs[14];
		vector<Function*>fs;
		//int numFuns;
		FunctionFactory(int numDim){
			//numFuns=13;
			fs.resize(13);
			//fs[0]=new F1(numDim);
			fs[0]=new F1(numDim);
			fs[1]=new F2(numDim);
			fs[2]=new F3(numDim);
			fs[3]=new F4(numDim);
			fs[4]=new F5(numDim);
			fs[5]=new F6(numDim);
			fs[6]=new F7(numDim);
			fs[7]=new F8(numDim);
			fs[8]=new F9(numDim);
			fs[9]=new F10(numDim);
			fs[10]=new F11(numDim);
			fs[11]=new F12(numDim);
			fs[12]=new F13(numDim);
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
			//return numFuns;
		}
		~FunctionFactory(){
			for(int i=0;i<getNumFunction();i++){
				delete fs[i];
			}
		}
};
FunctionFactory*FunctionFactory::instance=0;

#if ALGORITHM==4
int PPSO(int processId,int numProcess,Function*f,vector<double>&bestX,double &bestF){
	//MPI:
	const int TAG=99;
	MPI_Status status;
	//PSO:
	const int MaxStages=4;
	const int MaxFE=300000/MaxStages;
	const int NumAlgorithm=min(numProcess-1,4);
	int numDim=f->getNumDim();
	//
	PSO pso;
	int numP=20;

	vector<double>tmpbestXF;
	vector<double>bestXF;
	//
	tmpbestXF.resize(numDim+1);
	bestXF.resize(numDim+1);

	if(processId>=1 && processId<=4){
		pso.init(processId-1,numP);
		pso.begin(f);
	}
	for(int stage=1;stage<=MaxStages;stage++){
		switch(processId){
			case 0:	//Master:
				for(int i=1;i<=NumAlgorithm;i++){
					MPI_Recv(&tmpbestXF[0],tmpbestXF.size(),MPI_DOUBLE,i,TAG,MPI_COMM_WORLD,&status);
					if(i==1){
						bestXF=tmpbestXF;
					}else{
						if((tmpbestXF.back()>bestXF.back())^(f->getIsFindMin())){
							bestXF=tmpbestXF;
						}
					}
				}
				for(int i=1;i<=NumAlgorithm;i++){
					MPI_Send(&bestXF[0],bestXF.size(),MPI_DOUBLE,i,TAG,MPI_COMM_WORLD);
				}
				break;
			case 1:
			case 2:
			case 3:
			case 4:
				if(stage>1){
					pso.replaceParticle(drand()*numP,bestX,bestF);
				}
				pso.update(MaxFE/numP);
				pso.getOutput(bestX,bestF);
				bestX.push_back(bestF);
				MPI_Send(&bestX[0],bestX.size(),MPI_DOUBLE,0,TAG,MPI_COMM_WORLD);
				MPI_Recv(&bestX[0],bestX.size(),MPI_DOUBLE,0,TAG,MPI_COMM_WORLD,&status);
				bestF=bestX.back();
				bestX.pop_back();
				break;
			default:
				cerr<<"Error: numProcess must <=4"<<endl;
				break;
		}
	}
	//return only in master process.
	if(processId==0){
		bestF=bestXF.back();
		bestX=bestXF;
		bestX.pop_back();
		return 0;
	}
	return -1;
}
vector<double> runPPSO(int id,int idSize,Function*f,int maxRun){
	vector<double>results;
	results.resize(maxRun);
	for(int run=0;run<maxRun;run++){
		vector<double>bestX;
		double bestF;
		PPSO(id,idSize,f,bestX,bestF);
		results[run]=fabs(bestF-(f->getBest()));
	}
	return results;
}
#endif
vector<double> runSerialPSO(PSO &pso,Function*f,int maxRun){
	vector<double>results;
	const int MaxFE=300000;
	int numDim=f->getNumDim();
	vector<double>bestX;
	double bestF;
	results.resize(maxRun);
	for(int run=0;run<maxRun;run++){
		pso.solve(f,MaxFE/pso.getNumP(),bestX,bestF);
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
const char *getAlgorithmName(int algorithm){
	switch(algorithm){
		case 0:return "LPSO";
		case 1:return "GPSO";
		case 2:return "BPSO";
		case 3:return "CLPSO";
		case 4:return "PPSO";
	}
	return 0;
}
int old_main(int argc,char *argv[]){
	//int main(int argc,char *argv[]){
	const int maxRun=25;
	const int numDim=30;
	const int numP=20;
	FunctionFactory &funGenerator=FunctionFactory::Instance(numDim);
	const int numTestFunction=funGenerator.getNumFunction();
	srand(time(NULL));
#if ALGORITHM==4 
	int id,idSize;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	MPI_Comm_size(MPI_COMM_WORLD,&idSize);
	if(id==0){
		cout<<"Runing PPSO "<<maxRun<<"times."<<endl;
		printf("F\tmin\tmean\tstd\n");
		Tic::tic("begin");
	}
	for(int i=0;i<numTestFunction;i++){
		Function*f=funGenerator.getFunction(i);
		vector<double>results=runPPSO(id,idSize,f,maxRun);
		if(id==0){
			double min,max,mean,std;
			calStatistics(results,min,max,mean,std);
			printf("%s\t%g\t%g\t%g\n",f->getShortName(),min,mean,std);
		}
	}
	if(id==0){
		Tic::tic("end");
		cout<<"end."<<endl;
		cout<<endl;
	}
	MPI_Finalize();
#else
	PSO pso;
	pso.init(ALGORITHM,numP);
	cout<<"Runing PSO("<<pso.getAlgorithmName()<<") "<<maxRun<<"times."<<endl;
	printf("F\tmin\tmean\tstd\n");
	Tic::tic("begin");
	for(int i=0;i<numTestFunction;i++){
		//if(i==3){
		Function*f=funGenerator.getFunction(i);
		vector<double>results=runSerialPSO(pso,f,maxRun);
		double min,max,mean,std;
		calStatistics(results,min,max,mean,std);
		printf("%s\t%g\t%g\t%g\n",f->getShortName(),min,mean,std);
		//}
	}
	Tic::tic("end");
	cout<<"end."<<endl;
#endif
	return 0;
}
//int unused_main1(int argc,char *argv[]){
int main(int argc,char *argv[]){
	srand(time(NULL));
	const int MaxFE=300000;
	int numDim=30;
	int algorithm=2;
	const int numP=20;
	FunctionFactory &funGenerator=FunctionFactory::Instance(numDim);
	const int numTestFunction=funGenerator.getNumFunction();
	PSO pso;
	pso.init(algorithm,numP);
	Tic::tic("begin");
	cout<<"Runing PSO("<<pso.getAlgorithmName()<<") in serial way."<<endl;
	for(int i=0;i<numTestFunction;i++){
		//			if(i==0){
		vector<double>bestX;
		double bestF;
		Function*f=funGenerator.getFunction(i);
		pso.solve(f,MaxFE/numP,bestX,bestF);
		cout<<f->getName()<<" bestF:"<<bestF;
		cout<<" bestX:";printVec(bestX);
		cout<<endl;
		cout<<endl;
		//			}
	}
	Tic::tic("end");
	return 0;
}
int unused_main2(int argc,char *argv[]){
	//int main(){
	/*
	   Tic::tic("begin");
	   vector<double>a,b;
	   cout<<"NormD:";
	   const int times=1000000;
	   const int jtimes=10;
	   for(int i=0;i<times;i++){
	   for(int j=0;j<jtimes;j++){
	   a.push_back(NormD(0,1));
	//cout<<NormD(0,1)<<',';
	//		cout<<'.';
	}
	}
	cout<<endl;
	Tic::tic("NormD");
	for(int i=0;i<times;i++){
	for(int j=0;j<jtimes;j++){
	a.push_back(NormD1(0,1));
	//		cout<<'.';
	//cout<<NormD1(0,1)<<',';
	}
	}
	Tic::tic("NormD1");
	 */
	return 0;
}
//0:LPSO 1:GPSO 2:BPSO 3:CLPSO
//GPSO seems useless.
