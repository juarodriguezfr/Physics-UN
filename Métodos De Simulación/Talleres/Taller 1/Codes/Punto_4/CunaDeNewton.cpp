#include <iostream>
#include <cmath>
#include <fstream>
//#include "vector.h"
using namespace std;

//Constantes globales

const int N=3; // poner un numero mas grande
const double g=980; // unidades cgs
const double K=1e9;

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  //vector3D r,V,F;
  double theta, omega, tau;
  double m,R,l,I,x;
public:
  void Inicie(double theta0,double omega0,double m0,double R0, double l0, double x0);
  void BorreTorque(void){tau=0;};
  void SumeTorque(double tau0){tau+=tau0;};
  void Mueva_theta(double dt,double coeficiente);
  void Mueva_omega(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return x+l*sin(theta);}; //Inline
  double Gety(void){return -l*cos(theta);}; //Inline
  double Gettau(void){return tau;}; //Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double theta0,double omega0,double m0,double R0, double l0, double x0){
  theta=theta0; omega=omega0; l=l0; m=m0; R=R0; x=x0; I=m*l*l;
}
void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta+=omega*(dt*coeficiente);
}
void Cuerpo::Mueva_omega(double dt,double coeficiente){
  omega+=tau*(dt*coeficiente/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<x<<"+"<<l/7<<"*t*sin("<<theta<<"),-"<<l/7<<"*t*cos("<<theta<<")";
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeTorques(Cuerpo * Pendulo);
  void CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);    
};
void Colisionador::CalculeTorques(Cuerpo * Pendulo){
  int i,j; double tau0;
  //Borrar torques
  for(i=0;i<N;i++){
    Pendulo[i].BorreTorque();
    tau0=-Pendulo[i].l*Pendulo[i].m*g*sin(Pendulo[i].theta);
    Pendulo[i].SumeTorque(tau0);
  }
  //Calcular los torques entre todas las parejas de pendulos
  for(i=N-1;i>0;i--)
    CalculeTorqueEntre(Pendulo[i],Pendulo[i-1]);
}
void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  double s=(Pendulo2.Getx()+Pendulo2.R)-(Pendulo1.Getx()-Pendulo1.R); double F=0;
  if(s>0) F=K*pow(s,1.5); //Fuerza de Hertz
  Pendulo1.SumeTorque(F*Pendulo1.l); Pendulo2.SumeTorque(-F*Pendulo2.l);
  }

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'DosPendulos.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-5:12]"<<endl;
  cout<<"set yrange[-15:2]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}

void ExporteDatos(double t, double tau){
  ofstream archivo;
  archivo.open("torque.txt",ios::out);

  if(archivo.fail()){
    cout<<"No se pudo abrir el archivo"<<endl;
  }

  archivo<<t<<'\t'<<tau<<endl;
  archivo.close();
}

int main(){
  Cuerpo Pendulo[N];
  Colisionador Newton;
  double m0=100, l0=12, R=1.5;
  double T=2*M_PI*sqrt(l0/g);
  double t,tmax=1.5*T,dt=0.001;
  double tdibujo, tcuadro=T/200; // tcuadro > dt
  int i,k;
  //double Kvec[] = {10,20,30};
  
  //---------------(theta0,omega0,m0,R0,l0,x0)
  Pendulo[0].Inicie(-0.2618, 0, m0, R, l0, 0); // theta negativo
  for(i=1;i<N;i++)
     Pendulo[i].Inicie(0, 0, m0, R, l0, 2*R*i);

  /*ofstream archivo; // archivo con los datos de tau vs t
    archivo.open("prueba.txt",ios::out);*/

  /*for(k=0;k<2;k++){
    double K=Kvec[k];*/
  
  InicieAnimacion();
  
    for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
      //cout<<tdibujo<<'\t'<<tcuadro<<endl;
      //Dibujar
      if(tdibujo>tcuadro){
	InicieCuadro();
	for(i=0;i<N;i++) Pendulo[i].Dibujese();
	TermineCuadro();
	tdibujo=0;
	}       
    
      //cout<<Pendulo[1].Getx()<<" "<<Pendulo[1].Gety()<<endl;
      // Mover por PEFRL
      for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);
      Newton.CalculeTorques(Pendulo);
      for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
      for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
      Newton.CalculeTorques(Pendulo);
      for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
      for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Coeficiente2);
      Newton.CalculeTorques(Pendulo);
      for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
      for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
      Newton.CalculeTorques(Pendulo);
      for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
      for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);

      // función torque
      //cout<<t<<'\t'<<Pendulo[1].Gettau()<<endl;
      //archivo<<t<<'\t'<<Pendulo[1].Gettau()<<endl;
    }
    //}
    //archivo.close();
  
  return 0;
}
