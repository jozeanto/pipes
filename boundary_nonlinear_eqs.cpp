#include<iostream>
#include<math.h>

using namespace std;


#define Scalar double
#define pi 3.141592653589793238L


Scalar hu=0.5;
Scalar d=0.6;
Scalar Q=0.1;
Scalar g=9.81;


/*
 * Flujo entrada supercritico en lamina libre
 * -------------------------------------------
 */
Scalar ff(Scalar t){
  return  hu/d-(1-cos(t/2))/2 -(t-sin(t))/(16*sin(t/2));
}

Scalar dff(Scalar t){
  return  (-1./2*sin(t/2)*1./2
                    -( (1-cos(t))*16*sin(t/2) 
                    -(t-sin(t))*16*cos(t/2)*1./2 )  
                    /( pow(16*sin(t/2),2) ));
}

/*
 * Flujo de salida subcritico en lamina libre
 * --------------------------------------------
 */
Scalar fg(Scalar t){
  return  512.*Q/(g*pow(d,5))-pow(t-sin(t),3)/sin(t/2.);
}

Scalar dfg(Scalar t){
  return  -(3.*pow(t-sin(t),2)*(1-cos(t))*sin(t/2)-pow(t-sin(t),3)*cos(t/2)*1./2 )/pow(sin(t/2),2);
}


/*
 * Flujo de entrada subcritico en lamina libre
 * --------------------------------------------
 */
Scalar fh(Scalar t){
  return  512.*Q/(g*pow(d,5))-pow(t-sin(t),3)/sin(t/2.);
}

Scalar dfh(Scalar t){
  return  -(3.*pow(t-sin(t),2)*(1-cos(t))*sin(t/2)-pow(t-sin(t),3)*cos(t/2)*1./2 )/pow(sin(t/2),2);
}

Scalar newton ( Scalar (&f)(Scalar), Scalar (&df)(Scalar) ,Scalar x,int n, Scalar eps){
  int it=0;
  while ( (fabs(f(x)) > eps) && (it<n)){
    x=x-f(x)/df(x);
    it=it+1;
//     cout<<x<<endl;
  }
  return x;
}



int main(){
  
//   cout<<ff(0.001)<<endl;
  Scalar x=2.*pi*0.9999;
  newton(ff,dff,x,20,1e-12);
  
  x=2.*pi*0.9999;
  cout<<"newton value "<<newton(fg,dfg,x,40,1e-12)<<endl;
//   cout<<newton(ff,dff,0.001,20,1e-6)<<endl;

  
}
