#include<iostream>
#include<math.h>

using namespace std;


#define Scalar double
#define pi 3.141592653589793238L


Scalar hu=0.5;
Scalar d=0.6;


Scalar ff(Scalar t){
  return  hu/d-(1-cos(t/2))/2 -(t-sin(t))/(16*sin(t/2));
}

Scalar dff(Scalar t){
  return  (-1./2*sin(t/2)*1./2
                    -( (1-cos(t))*16*sin(t/2) 
                    -(t-sin(t))*16*cos(t/2)*1./2 )  
                    /( pow(16*sin(t/2),2) ));
}


Scalar newton ( Scalar (&f)(Scalar), Scalar (&df)(Scalar) ,Scalar x,int n, Scalar eps){
  int it=0;
  while ( (fabs(f(x)) > eps) && (it<n)){
    x=x-f(x)/df(x);
    it=it+1;
    cout<<x<<endl;
  }
  return x;
}



int main(){
  
//   cout<<ff(0.001)<<endl;
  Scalar x=2.*pi*0.9999;
  newton(ff,dff,x,20,1e-12);
//   cout<<newton(ff,dff,0.001,20,1e-6)<<endl;

  
}
