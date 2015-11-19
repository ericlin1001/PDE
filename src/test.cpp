#include "include/settingParser.h"
#include<iostream>
#include<vector>
using namespace std;

#define T(fun,input) {fun(input,buffer);cout<<#fun"("#input")="<<buffer<<endl;}
#define Trace(m) cout<<#m"="<<(m)<<endl;
template<class T>
void printV(const vector<T>&arr){
	cout<<"[";
	for(int i=0;i<arr.size();i++){
		if(i!=0)cout<<",";
		cout<<arr[i];
	}
	cout<<"]";
}
template<class T>
void printVV(const vector<vector<T> >&arr){
	cout<<"[";
	for(int i=0;i<arr.size();i++){
		if(i!=0)cout<<",";
		const vector<T>&brr=arr[i];
		cout<<"[";
		for(int j=0;j<brr.size();j++){
			if(j!=0)cout<<",";
			cout<<brr[i];
		}
		cout<<"]";
	}
	cout<<"]";
}
int main(){
	SettingParser setting("setting.json");
	char buffer[100];
	T(setting.getAuthorInfo,"Email");
	T(setting.getAuthorInfo,"Author");

	T(setting.getString,"Name");

	Trace(setting.getDouble("F"));
	Trace(setting.getDouble("CR"));
	Trace(setting.getInt("NumDim"));
	Trace(setting.getInt("MaxFEs"));
	Trace(setting.getInt("PopSize"));
	Trace(setting.getInt("NumSubPop"));
	Trace(setting.getInt("NumGenPerComm"));
	vector<vector<double> > range;
	setting.getBiVector("Range",range);
	printVV(range);
	cout<<endl;
	return 0;
}
