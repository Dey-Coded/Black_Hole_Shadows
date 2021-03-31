#include "math.h"
#include <iostream>
#include <vector>
using namespace std;

double k1,k2,k3,k4,l1,l2,l3,l4;

double f(double u){
	return -u*(1-3*u/2);
}

//Function for increment in u, v using RK4 scheme
void rkstep(double& u, double& v, double h){
	k1 = h*v;
	l1 = h*f(u);
	k2 = h*(v+0.5*l1);
	l2 = h*f(u+0.5*k1);
	k3 = h*(v+0.5*l2);
	l3 = h*f(u+0.5*k2); 
	k4 = h*(v+l3);
	l4 = h*f(u+k3);
	u += (k1+2*k2+2*k3+k4)/6;
	v += (l1+2*l2+2*l3+l4)/6;
}

int main(){
	vector<double> X;
	vector<double> Y;

	//The initial conditions for light emitted at infinity parallel to (theta=0,phi=0) direction towards the black hole
	double u=0, phi=0, r;
	//The initial value of v=du/d(phi) = beta(in this case)
	double v = sqrt(0.25);
	double x,y;
	do{
		rkstep(u,v,0.02);
		phi += 0.02; //Since phi increases along the geodesics considered
		r = 1/u;
		x = r*cos(phi);
		y = r*sin(phi);
		//The following condition is put for imaging the trajectory properly in gnuplot so that only a small but relevant frame is shown
		if(fabs(x)<6){
			X.push_back(x);
			Y.push_back(y);
		}
	}while(x>-4 && y>-4 && r>1);

	//The values will be printed in the output file
	for (int i = 0; i < X.size(); ++i)
	{
		cout<<X[i]<<" "<<Y[i]<<endl;
	}

	return 0;
}
