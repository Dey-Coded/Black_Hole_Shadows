#include "math.h"
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
using namespace Eigen;
using namespace cv;
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

//Function for calculating the 3d rotation matrix. The axis of rotation is (cos(phi),sin(phi),0) and theta is the angle of rotation
void calc_Rmat(MatrixXd& R, double theta, double phi){
	double n1 = cos(phi);
	double n2 = sin(phi);
	R(0,0) = cos(theta) + n1*n1*(1-cos(theta));
	R(0,1) = n1*n2*(1-cos(theta));
	R(0,2) = n2*sin(theta);
	R(1,0) = n1*n2*(1-cos(theta));
	R(1,1) = cos(theta) + n2*n2*(1-cos(theta));
	R(1,2) = -n1*sin(theta);
	R(2,0) = -n2*sin(theta);
	R(2,1) = n1*sin(theta);
	R(2,2) = cos(theta);
}

int main(){
	//Setting the parameters for the camera:
	//Radial coordinate = rc
	//Azimutal coordinate = phic
	//Angular resolution = max_res
	//Pixel resolution = 100
	double rc=10, phic=M_PI/8, max_res=M_PI/4, np=200;
	double f, psi, theta, beta2;
	double u, v, phi, r, x, y, z, xp, yp,ra,dir=1;

	//Initially an cv::image object initiall coloured all black
	Mat M((int)np,(int)np,CV_8UC3,Vec3b(0,0,0));

	//R=Rotation matrix, Xp=geodesic coordinate in (x,y) plane where it is computed, X=rotated coordinates corresponding to 
	// the pixel's direction
	MatrixXd R = MatrixXd::Zero(3,3);
	MatrixXd Xp = MatrixXd::Zero(3,1);
	MatrixXd X = MatrixXd::Zero(3,1);

	//Assigning the mid-point pixel coordinate and calculating focus 
	int im, jm;
	if(((int)np)%2==0){
		f = (np/2-1)/tan(max_res/2);
		im = np/2 - 1; jm = np/2 - 1;
	}
	else{
		f = 0.5*(np-1)/tan(max_res/2);
		im = (np-1)/2; jm = (np-1)/2;
	}
	cout<<"focus: "<<f<<endl;
	for (int i = 0; i < np; ++i)
	{
		for (int j = 0; j < np; ++j)
		{
			if(i==im && j==jm){
				continue;
			}
			//The direction of emmision of ray for pixel (i,j)
			psi = atan(sqrt(pow(i-im,2)+pow(j-jm,2))/f);
			
			//initial parameters for the null geodesic
			r=rc; u=1/rc; phi=phic;
			beta2 = (r-1)/(pow(sin(psi),2)*pow(r,3));
			v = sqrt(beta2 - u*u*(1-u));

			//Calculating the angle by which the geodesic should be rotated from the xy plane 
			if(i==im){
				if(j<jm)
					theta = M_PI/2;
				else
					theta = -M_PI/2;
			}
			else if(j==jm){
				if(i<im){
					theta = 0;
				}
				else{
					theta = 0;
					dir = -1;
				}
			}
			else{
				if(i>im && j<jm){
					theta = -atan(((double)(jm-j))/((double)(i-im)));
					dir = -1;
				}
				else if(i>im && j>jm){
					theta = atan(((double)(j-jm))/((double)(i-im)));
					dir = -1;
				}
				else{
					theta = atan(((double)(jm-j))/((double)(im-i)));
				}
			}
			//Computing the rotation matrix
			calc_Rmat(R,theta,phic);

			do{
				rkstep(u,v,0.02);
				phi += dir*0.02; //Based on the initial direction of beam, phi is either incremented or decremented 
				r = 1/u;
				xp = r*cos(phi);
				yp = r*sin(phi);				
				Xp(0,0)=xp; Xp(1,0)=yp; Xp(2,0)=0;//Coordinates in xy plan
				X = R*Xp;//Coordinates in rotated plan
				if(r>1 && r<2 && fabs(X(1,0))<0.05){//The conditions for ray to way in the annular region of accretion disk 
					M.at<Vec3b>(i,j)[0] = 19;
					M.at<Vec3b>(i,j)[1] = 90;
					M.at<Vec3b>(i,j)[2] = 242;
					break;
				}
				if(X(0,0)<-2){
					break;
				}
				
			}while(r>1);

		}
	}
	namedWindow("image",CV_WINDOW_FREERATIO);
	imshow("image",M);
	waitKey(0);
	imwrite("phi22_5.png",M);
	return 0;
}

