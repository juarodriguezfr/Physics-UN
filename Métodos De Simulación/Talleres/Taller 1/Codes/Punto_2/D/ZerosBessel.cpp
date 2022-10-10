#include<iostream>
#include<cmath>
using namespace std;


//_________________________________________________________________________
//Sistema de ecuaciones Diferenciales de primer orden
//_________________________________________________________________________

double f1(double t, double x1, double x2, double l)
{
  return -(x1/t)-(l*l*x2);
}

double f2(double t, double x1, double x2, double l)
{
  return x1;
}


//_________________________________________________________________________
//Solución de sistema de EDO por Runge Kutta 4
//_________________________________________________________________________

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


//_________________________________________________________________________
//función R(1,lambda) con dl=0.001
//_________________________________________________________________________

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


//_________________________________________________________________________
//Método de bisección con precisión 1e-9 para la solución numérica
//_________________________________________________________________________

double CerosPorBiseccion(double a, double b)
{
  const double ErrMax=1e-9;
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


//_________________________________________________________________________
//Argumento de la integral de Bessel
//_________________________________________________________________________

double ff(double alpha, double x, double t)
{
  return cos(alpha*t-x*sin(t));
}


//_________________________________________________________________________
//Método de integración de Simpson
//_________________________________________________________________________

double IntegralSimpson(double alpha, double x, double a, double b, int n)
{
  double t, h,suma;
  int i;
  n*=2;
  h=(b-a)/n;
  
  for(suma=0,i=0;i<=n;i++)
    {
      t=a+i*h;
      if (i==0 || i==n) //si es el primero o el ultimo
	suma+=ff(alpha,x,t);
      else if(i%2==0) // si es par
	suma+=2*ff(alpha,x,t);
      else            //si es impar
	suma+=4*ff(alpha,x,t);
    }
  return suma*h/3;
}


//_________________________________________________________________________
//Cálculo funciones de Bessel
//_________________________________________________________________________

double Bessel(double alpha, double x)
{
  double a=0,b=M_PI;
  int n=100;
  return IntegralSimpson(alpha,x,a,b,n)/M_PI;
}


//_________________________________________________________________________
//Método de bisección con precisión 1e-9 para las función de Bessel
//_________________________________________________________________________

double CerosPorBiseccionBessel(double a, double b)
{
  const double ErrMax=1e-9;
  double x, m,fa,fm;
  fa=Bessel(0,a);
  while(b-a>=ErrMax)
    {
      m=(b+a)/2; 
      fm=Bessel(0,m);

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

//_________________________________________________________________________
//Programa que calcula los ceros por ambos métodos
//_________________________________________________________________________

int main()
{
  double a, b;

  cout << "Ceros Solución RK4: " << "," << "Ceros Solución Bessel: " << ","<< "Valor de la función numérica " << "," << "Valor de la función de Bessel " <<endl;
  for(a=0,b=3;a<15;)
  {
    cout << CerosPorBiseccion(a,b) << "," << CerosPorBiseccionBessel(a,b) << "," << f(CerosPorBiseccion(a,b)) << "," <<Bessel(0,CerosPorBiseccionBessel(a,b))<< endl;
    a+=3;
    b+=3;
  }
  return 0;
}
