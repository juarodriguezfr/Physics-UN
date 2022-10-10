#include<iostream>
#include<cmath>

using namespace std;

double f1(double t, double x1, double x2, double l)
{
  return -(x1/t)-(l*l*x2);
}

double f2(double t, double x1, double x2, double l)
{
  return x1;
}


void UnPasoDeRKs(double & t0, double & x10, double & x20, double dt, double & l)
{
  double dx11,dx21,dx31,dx41;
  double dx12,dx22,dx32,dx42;
  
  dx11=dt*f1(t0,x10,x20,l);
  dx12=dt*f2(t0,x10,x20,l);
  
  dx21=dt*f1(t0+dt/2,x10+dx11/2,x20+dx12/2,l);
  dx22=dt*f2(t0+dt/2,x10+dx11/2,x20+dx12/2,l);
  
  dx31=dt*f1(t0+dt/2,x10+dx21/2,x20+dx22/2,l);
  dx32=dt*f2(t0+dt/2,x10+dx21/2,x20+dx22/2,l);
  
  dx41=dt*f1(t0+dt,x10+dx31,x20+dx32,l);
  dx42=dt*f2(t0+dt,x10+dx31,x20+dx32,l);
  
  x10+=(dx11+2*(dx21+dx31)+dx41)/6;
  x20+=(dx12+2*(dx22+dx32)+dx42)/6;
  t0+=dt;
}

double f(double l)
{
  double t,x1,x2;
  double dt=0.001;

  for(t=0.01,x1=0,x2=1; t<=1; )
      {
          UnPasoDeRKs(t,x1,x2,dt,l);
      }
   return x2;
}

//Error tolerado
const double ErrMax=1e-9;

//Método de bisección
double CerosPorBiseccion(double a, double b)
{
  double x, m,fa,fm;
  fa=f(a);
  while(b-a>=ErrMax)
    {
      m=(b+a)/2; 
      fm=f(m);

      if(fa*fm>0)
      {
        a=m;
        fa=fm;
      }
      else
      b=m;
    }
  return (a+b)/2;
}


int main()
{
  double a, b;
  for(a=0,b=3;a<15;)
  {
    //cout <<"El cero es " <<CerosPorBiseccion(a,b) <<endl;
    cout << CerosPorBiseccion(a,b) << " " << f(CerosPorBiseccion(a,b)) << endl;
    a+=3;
    b+=3;
  }
  return 0;
}
