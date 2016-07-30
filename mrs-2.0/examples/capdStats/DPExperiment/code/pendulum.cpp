/*   File: pendulum.cpp
     
     Step 1: Given an initial conditions x_0, a parameter p, and times 
     {\t_i}_0^N_, we generate data points {\y_i}_0^N. Note that the times 
     may be interval-valued.
     
     Step 2: Given the data {\t_i; \y_i}_0^N and a parameter q, we check 
     for consistency.
      
     All of the above should be looped so we can check several different
     trajectories (one for each new initial condition).

     Ex: ./pendulum 127 238 10 9.81
     
     Author: Warwick Tucker <warwick@math.uu.se>
     Latest edit: Sat Jul 30 19:37:38 CEST 2016 (raaz testing in git)
*/

#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// Our vector field for the double pendulum. Should be tweaked into something else...
void pendulumVectorField(Node t, Node in[], int dimIn, Node out[], int dimOut, Node params[], int noParams)
{
	Node g  = params[0];
	Node l1 = params[1];
	Node l2 = params[2];
	Node m1 = params[3];
  Node m2 = params[4]; 
  Node F1 = params[5];
  Node F2 = params[6];
  
  Node I1 = m1*l1^2/12.0;
  Node I2 = m2*l2^2/12.0;
  Node h1 = 0.5*l1;                    
  Node h2 = 0.5*l2;                   
  Node mu1 = 1.0/(I1 + m2*l1^2);
  Node mu2 = 1.0/(I2 + m2*h2^2);
  Node w12 = g*(m1*h1 + m2*l1)/(I1 + m2*l1^2);
  Node w22a = g*m2^2;
  Node w22 = w22a/(I2 + m2*h2^2);
	Node x1 = in[0];
	Node x2 = in[1];	
	Node x3 = in[2];
	Node x4 = in[3];
		
	Node e1a = w12*sin(x1) + F1*mu1*x3;        // Note: the nodes can not
	Node e1b = m2*mu1*l1*h2*sin(x1-x2)*x4^2;   // consist of heavily
	Node e1 = e1a + e1b;                       // compound expressions.
	Node e2 = m2*mu1*l1*h2*cos(x1-x2);

  Node f1a = w22*sin(x2) + F2*mu2*x4;
  Node f1b = m2*mu2*l2*h2*sin(x1-x2)*x3^2;
	Node f1 = f1a - f1b;
	Node f2 = m2*mu2*l2*h2*cos(x1-x2);

	out[0] = x3;
	out[1] = x4;
	out[2] = (e1 - f1*e2)/(e2*f2 - 1.0);
	out[3] = (f1 - e1*f2)/(e2*f2 - 1.0);
}

void printSample(const double &t, const IVector &x)
{
	cout.setf(ios::scientific);
  cout.setf(ios::showpos);
  cout.precision(10);
  cout << t << " ";
  for(int i=0; i<x.dimension(); ++i)
    cout << x[i].leftBound() << " " << x[i].rightBound() << " ";
  cout << endl; 
}

int main(int argc, const char*argv[]){

 if (argc != 5) {
    cerr << "Syntax: pendulum theta1 theta2 Nt g" << endl;
    exit(0); 
  }
  
  double theta1(atof(argv[1])); // Initial angle 1 in degrees.
  double theta2(atof(argv[2])); // Initial angle 2 in degrees.
  int    Nt(atoi(argv[3]));     // The number of time samples.
  double g(atof(argv[4]));      // Guess for gravity.
  
  // Set remaining parameter values. Can be double on interval type.
  double l1(1.0);
  double l2(1.0);
  double m1(1.0);
  double m2(1.0);
  double F1(0.2);
  double F2(0.2);  
  
  // Specify the vector field, and set the parameters.
  int dimIn=4, dimOut=4, noParam=7;
	IMap pendulum(pendulumVectorField,dimIn,dimOut,noParam);
  pendulum.setParameter(0, g);
  pendulum.setParameter(1,l1);
  pendulum.setParameter(2,l2);
  pendulum.setParameter(3,m1);
  pendulum.setParameter(4,m2);
  pendulum.setParameter(5,F1);
  pendulum.setParameter(6,F2);
  
  // Define an instance of the ODE solver.
  ITaylor solver(pendulum,20);
  solver.setAbsoluteTolerance(1e-12);
  solver.setRelativeTolerance(1e-12);
  ITimeMap timeMap(solver);

  // Specify the initial conditions.
  IVector x(4);
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = theta1/360*2*DInterval::pi();
  x[3] = theta2/360*2*DInterval::pi();

	DInterval dx(-0.001,+0.001); // Just checking.
  x[3] += dx;  

  double T0 = 0.0; // Initial time.
  double Tf = 1.0; // Final time.
  
  // Use a doubleton representation of initial condition.
  C0HORect2Set x0(x,T0);
  // Define a functional object...
  ITimeMap::SolutionCurve solution(T0);  
  // ...and integrate.
  timeMap(Tf,x0,solution);
  
  // Get intermediate points on the trajectory.
  printSample(T0, solution(T0));
  //cout << "t = " << T0 << "; x = " << solution(T0) << endl;    
  for (int i = 1; i < Nt; ++i) { 
    double ti = 1.0/Nt*i;
    printSample(ti, solution(ti));
  }
  printSample(Tf, solution(Tf));
  
  return 0;
}
