#include<iostream>
#include<cmath>

using namespace std;

const double g=0.08;

//Se definen las dos funciones
double f1(double t, double x1, double x2, double b)
{
  return -b*x1*x2;
}

double f2(double t, double x1, double x2, double b)
{
  return (b*x1*x2)-(g*x2);
}

void UnPasoDeRKs(double & t0, double & x10, double & x20, double dt, double & b)
{
  double dx11,dx21,dx31,dx41;
  double dx12,dx22,dx32,dx42;
  
  dx11=dt*f1(t0,x10,x20,b);
  dx12=dt*f2(t0,x10,x20,b);
  
  dx21=dt*f1(t0+dt/2,x10+dx11/2,x20+dx12/2,b);
  dx22=dt*f2(t0+dt/2,x10+dx11/2,x20+dx12/2,b);
  
  dx31=dt*f1(t0+dt/2,x10+dx21/2,x20+dx22/2,b);
  dx32=dt*f2(t0+dt/2,x10+dx21/2,x20+dx22/2,b);
  
  dx41=dt*f1(t0+dt,x10+dx31,x20+dx32,b);
  dx42=dt*f2(t0+dt,x10+dx31,x20+dx32,b);
  
  x10+=(dx11+2*(dx21+dx31)+dx41)/6;
  x20+=(dx12+2*(dx22+dx32)+dx42)/6;
  t0+=dt;
}

int main()
{
  double t,x1,x2,b; 
  double dt=0.02; //Paso de integracion

  for(b=0.09;b<1;) //Ciclo para variar beta
  {
    for(t=0,x1=0.999,x2=0.001; t<4000;)
      {
        UnPasoDeRKs(t,x1,x2,dt,b);
      }
    cout<< b/g << " " << x1 << endl; //Imprime ultimo valor de s
    b+=0.01; //añade el paso para beta
  }
  return 0;
}
