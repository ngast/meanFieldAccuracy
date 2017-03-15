#include<iostream>
#include<cstdlib>

double r_unif(){
  return (rand()/(1.+RAND_MAX));
}
void simulate(double lambda,int T=100000000){
  int X1=0,X2=0;
  double rate = 2*lambda+2;
  for(int t=0;t<T;t++){
    double u = rate*r_unif();
    if (u < 2*lambda){
      if (X1<=X2) X1++; else X2++;
    }
    else if (u < 2*lambda + 1)
      {if (X1>0) X1--;}
    else
      {if (X2>0) X2--;}
    if (t%1000==0)
      std::cout << X1 << " " << X2 <<"\n";
  }
}
int main(){
  simulate(0.80);
}
