//Author:Ericlin 
//Date:Tue Oct 20 21:43:43 CST 2015
//#define DEBUG
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

#include "include/functions.h"
#include "include/utils.h"
#include "inlcude/settingParser.h"

#if ALGORITHM==4
#define OMPI_IMPORTS
#include "mpi.h"
#endif

#define Trace(m) {cout<<#m"="<<(m)<<endl;}

class DMPDE{
	double F,CR;
	int NumDim;
	Function *f;
	double xmin,xmax;
	SettingParser setting;
	public:
	SettingParser &getSetting(){
		return setting;
	}
	Function*getF()const{return f;}
	void setSetting(SettingParser &p){
		setting=p;
		initSetting(p.getDouble("F"),p.getDouble("CR"),p.getInt("NumDim"));
	}
	void initSetting(double F,double CR,int NumDim){
		this->F=F;
		this->CR=CR;
		this->NumDim=NumDim;
	}
	void setFunction(Function*f){
		this->f=f;
		assert(NumDim==f->getNumDim());
		xmin=f->getRange(0);
		xmax=f->getRange(1);
	}
	//the last value of X is FX.
	//the size of X is NumDim+1.
	//x is the array of X.
	void initPopulation(vector<vector<double> >&x,int begin,int size,int &bestI){
		assert(begin>=0&&begin+size<=x.size());
		int end=begin+size;
		for(int i=begin;i<end;i++){
			for(int j=0;j<NumDim;j++){
				x[i][j]=drand(xmin,xmax);
			}
			double fx=(*f)(&x[i][0],NumDim);
			x[i][NumDim]=fx;
			if(f->isFBetter(fx,x[bestI][NumDim])){
				bestI=i;
			}
		}
	}
	void update(vector<vector<double> >&x,int begin,int size,int &bestI){
		//update locally
		assert(size>=3);
		assert(begin>=0&&begin+size<=x.size());
		assert(begin<=bestI&&bestI<begin+size);
		if(x[0].size()!=(NumDim+1)){
			Trace(begin);
			Trace(x.size());
			Trace(x[1].size());
			Trace(x[0].size());
		}
		assert(x[0].size()==(NumDim+1));
		vector<double>tx;
		RandomPermutation perm(size);
		tx.resize(NumDim+1);
		int end=begin+size;
		for(int i=begin;i<end;i++){
			perm.generate();
			int a=begin+perm.next();
			int b=begin+perm.next();
			int c=begin+perm.next();
			int randDim=rand()%NumDim;
			for(int j=0;j<NumDim;j++){
				if(j==randDim||drand()<CR){
					tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
					if(tx[j]<xmin || tx[j]>xmax){
						tx[j]=drand(xmin,xmax);
					}
				}else{
					tx[j]=x[i][j];
				}
			}
			double ftx=(*f)(&tx[0],NumDim);
			tx[NumDim]=ftx;
			if(f->isFBetter(ftx,x[i][NumDim])){
				x[i]=tx;
				if(f->isFBetter(ftx,x[bestI][NumDim])){
					bestI=i;
				}
			}
		}
	}
	//update base on the infomation of all.
	void updateG(vector<vector<double> >&x,int begin,int size,int&bestI){
		//update globally
		assert(size>=3);
		assert(begin>=0&&begin+size<=x.size());
		assert(begin<=bestI&&bestI<begin+size);
		assert(x[0].size()==(NumDim+1));
		vector<double>tx;
		RandomPermutation perm(x.size());
		int end=begin+size;
		tx.resize(NumDim+1);
		for(int i=begin;i<end;i++){
			perm.generate();
			int a=perm.next();
			int b=perm.next();
			int c=perm.next();
			int randDim=rand()%NumDim;
			for(int j=0;j<NumDim;j++){
				if(j==randDim||drand()<CR){
					tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
					if(tx[j]<xmin || tx[j]>xmax){
						tx[j]=drand(xmin,xmax);
					}
				}else{
					tx[j]=x[i][j];
				}
			}
			double ftx=(*f)(&tx[0],NumDim);
			tx[NumDim]=ftx;
			if(f->isFBetter(ftx,x[i][NumDim])){
				x[i]=tx;
				if(f->isFBetter(ftx,x[bestI][NumDim])){
					bestI=i;
				}
			}
		}
	}
};


void serialOut(const vector<vector<double> >&x,int begin,int size,double *xbuf){
	const int end=begin+size;
	const int len=x[begin].size();
	int counter=0;
	for(int i=begin;i<end;i++){
			for(int j=0;j<len;j++){
				//xbuf[i*len+j]=x[i][j];
				xbuf[counter++]=x[i][j];
			}
	}
}
void serialIn(vector<vector<double> >&x,const double *xbuf){
	const int end=x.size();
	const int len=x[0].size();
	for(int i=0;i<end;i++){
		for(int j=0;j<len;j++){
			x[i][j]=xbuf[i*len+j];
		}
	}
}
int doDMPDE(int id,int idSize,DMPDE*de,vector<double>&bestX,double &bestF){
	srand(time(NULL));//srand in every process.
	//MPI:
	const int TAG=99;
	const int TAG2=98;
	MPI_Status status;
	//DE:
	SettingParser s=de->getSetting();
	const int MaxFEs=s.getInt("MaxFEs");
	const int PopSize=s.getInt("PopSize");
	const int NumSubPop=s.getInt("NumSubPop");
	const int NumGenPerComm=s.getInt("NumGenPerComm");
	const int NumDim=s.getInt("NumDim");
	//
	const int size=PopSize/NumSubPop;
	const int begin=id*size;
	const int MaxGen=MaxFEs/PopSize-1;//since the init use one generation.

	int bestI=begin;//must be initialized.

	assert(idSize==NumSubPop);
	assert(PopSize%NumSubPop==0);
	//
	vector<vector<double> > x;
	x.resize(PopSize);
	for(int i=0;i<PopSize;i++){
		x[i].resize(NumDim+1);
	}
	de->initPopulation(x,begin,size,bestI);
	double *xbuf=new double[PopSize*(NumDim+1)];
	double *xsendbuf=new double[size*(NumDim+1)];
	assert(xbuf!=NULL);

	const int sendCount=size*(NumDim+1);
	const int sendBegin=sendCount*id;
	int g=0;
	while(true){
		//locally update
		for(int l=0;l<NumGenPerComm;l++){
			g++;
			de->update(x,begin,size,bestI);
			if(g>=MaxGen)break;
		}
		//exchange info globally.
		serialOut(x,begin,size,xsendbuf);
		MPI_Allgather(xsendbuf,sendCount,MPI_DOUBLE,
				xbuf,sendCount,MPI_DOUBLE,MPI_COMM_WORLD);
		serialIn(x,xbuf);

		//globally update
		g++;
		de->updateG(x,begin,size,bestI);
		if(g>=MaxGen)break;
	}
		serialOut(x,begin,size,xsendbuf);
		MPI_Allgather(xsendbuf,sendCount,MPI_DOUBLE,
				xbuf,sendCount,MPI_DOUBLE,MPI_COMM_WORLD);
		serialIn(x,xbuf);
	delete []xsendbuf;
	delete []xbuf;

	int tmpBestI=0;
	for(int i=0;i<PopSize;i++){
		if(de->getF()->isFBetter(x[i][NumDim],x[tmpBestI][NumDim])){
			tmpBestI=i;
		}
	}
	bestX=x[tmpBestI];
	bestF=bestX.back();
	bestX.pop_back();
}
vector<double> runDMPDE(int id,int idSize,DMPDE*de,int MaxRun){
	vector<double>results;
	results.resize(MaxRun);
	for(int run=0;run<MaxRun;run++){
		vector<double>bestX;
		double bestF;
		doDMPDE(id,idSize,de,bestX,bestF);
		results[run]=fabs(bestF-(de->getF()->getBest()));
	}
	return results;
}

int main(int argc,char *argv[]){
	srand(time(NULL));
	//parse the setting.
	SettingParser setting("setting.json");
	char author[100];
	char name[100];
	setting.getAuthorInfo("Author",author);
	setting.getString("Name",name);
	const int NumDim=setting.getInt("NumDim");
	const int MaxRun=setting.getInt("MaxRun");

	//cout<<"Reading Setting [OK]"<<endl;
	//set up the test functions;
	FunctionFactory &funGenerator=FunctionFactory::Instance(NumDim);
	const int numTestFunction=funGenerator.getNumFunction();

	//set up the MPI
	int id,idSize;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	MPI_Comm_size(MPI_COMM_WORLD,&idSize);

	//start run
	if(id==0){
		cout<<"Runing "<<name<<"(Author:"<<author<<") "<<MaxRun<<" times."<<endl;
		printf("F\tmean\tstd\n");
		Tic::tic("begin");
	}
	for(int i=0;i<numTestFunction;i++){
		Function*f=funGenerator.getFunction(i);
		DMPDE de;
		de.setSetting(setting);
		de.setFunction(f);
		vector<double>results=runDMPDE(id,idSize,&de,MaxRun);
		if(id==0){
			double min,max,mean,std;
			calStatistics(results,min,max,mean,std);
			printf("%s\t%g\t%g\n",f->getShortName(),mean,std);
		}
	}
	//end run.

	//do some end up works.
	if(id==0){
		Tic::tic("end");
		cout<<endl;
	}
	MPI_Finalize();
	return 0;
}
