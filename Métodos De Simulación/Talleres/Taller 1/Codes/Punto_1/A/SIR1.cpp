#include<iostream>
#include<cmath>

using namespace std;

const double b=0.35;
const double g=0.08;

//Se definen las dos funciones
double f1(double t, double x1, double x2)
{
  return -b*x1*x2;
}

double f2(double t, double x1, double x2)
{
  return (b*x1*x2)-(g*x2);
}

double f3(double t, double x2, double x3)
{
  return (g*x2);
}


void UnPasoDeRKs(double & t0, double & x10, double & x20, double dt)
{
  double dx11,dx21,dx31,dx41;
  double dx12,dx22,dx32,dx42;
  
  dx11=dt*f1(t0,x10,x20);
  dx12=dt*f2(t0,x10,x20);
  
  dx21=dt*f1(t0+dt/2,x10+dx11/2,x20+dx12/2);
  dx22=dt*f2(t0+dt/2,x10+dx11/2,x20+dx12/2);
  
  dx31=dt*f1(t0+dt/2,x10+dx21/2,x20+dx22/2);
  dx32=dt*f2(t0+dt/2,x10+dx21/2,x20+dx22/2);
  
  dx41=dt*f1(t0+dt,x10+dx31,x20+dx32);
  dx42=dt*f2(t0+dt,x10+dx31,x20+dx32);
  
  x10+=(dx11+2*(dx21+dx31)+dx41)/6;
  x20+=(dx12+2*(dx22+dx32)+dx42)/6;
  t0+=dt;
}

void UnPasoDeRK4(double & t0, double & x0 ,double & x2, double dt)
{
  double dx1,dx2,dx3,dx4;
  dx1=dt*f3(t0,x2,x0);
  dx2=dt*f3(t0+dt/2,x2,x0+dx1/2);
  dx3=dt*f3(t0+dt/2,x2,x0+dx2/2);
  dx4=dt*f3(t0+dt,x2,x0+dx3);  
  x0+=(dx1+2*(dx2+dx3)+dx4)/6;
  t0+=dt;
}

int main()
{
  double t,x1,x2,x3;
  double dt=0.01;
  
  for(t=0,x1=0.999,x2=0.001,x3=0; t<300; )
    {
      cout<<t<< " "<<x1<< " " << x2 << " " << x3 <<endl;
      UnPasoDeRKs(t,x1,x2,dt);
      UnPasoDeRK4(t,x3,x2,dt);
    }
  return 0;
}
