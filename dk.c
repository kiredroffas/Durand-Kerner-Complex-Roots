/*  Erik Safford
 *  Programming Assignment 1 : Durand-Kerner Method
 *  CS 330
 *  September 2018  */

/* - This program implements the Durand-Kerner method for finding all n (complex) roots of a nth degree polynomial 
 * - n is equal to the highest degree of the polynomial
 * - Program reads in polynomials until something that is not %lf %lf is entered (i.e. a letter), or can pipe a .in file  */

#include <stdio.h>
#include <math.h>
#include <complex.h>

#define M_PI 3.14159265358979323846

void durandKerner(); //Prototypes
double findThetaj();

int main() {

	double complex cList[20];  //cList is a list of coefficients
	double complex z; //z used to combine real/imaginary parts into one
	double x,y; //x = real, y = imaginary
	int n=0; //n is number of complex roots/highest degree of polynomial

//------Read in Coefficients------------------------------------------------
	printf("Enter coefficiens and enter 'd' when done, or pipe a .in file\n");	
	while(scanf("%lf %lf",&x,&y) == 2)  { //Read coefficients from stdin
		cList[n] = (x + y*I);
		n++;
	}
	x = 1;  //Cn = 1, and is not specified
	y = 0;
	z = (x + y*I);
	cList[n] = z; //Store in cList[]
	
	durandKerner(cList,n);
}
//===========================================================================
double R=0; //R holds radius of circle in complex plane that contains f(z) roots
double complex z[20]; //z[] holds initial guesses for the n roots
double complex deltaZ[20]; //deltaZ[] holds values computed w/Horner's method and Equation 2
double deltaZMax;
double epsilon = 1e-6; //E = 1x10^-6,targets 5 or 6 digits of precision
double complex QsubJ,fz; //fz holds the function evaluated at some z
int kMAX = 50;

void durandKerner(double complex cList[20],int n) {
	
	for(int j=0;j < n;j++) { //Start Equation 5
		if(cabs(cList[j]) > R) { //Finds largest coefficient
			R = cabs(cList[j]);
		}
	} //Adds 1 to largest coefficient to find R
	R = 1 + R;  //End Equation 5

	for(int j=0;j < n;j++) { //Start Equation 6
		z[j] = ( cos(findThetaj(n,j)) + (I*sin(findThetaj(n,j))) )*R;
	} //End Equation 6


	for(int k=1;k <= kMAX;k++) { //Start Durand-Kerner algorithm
		
		printf("iter %d \n",k);
		for(int i=0;i < n;i++) {  //Print z[]
                	printf("z[%d] = %0.10f + %0.10f*I\n",i,creal(z[i]),cimag(z[i]));
                fflush(stdout);
        	}

		deltaZMax = 0; //Set change in ZMax to 0

		for(int j=0;j < n;j++) { //Start Equation 2
			
			QsubJ = 1; //Initialize to 1 so can be multiplied
			for(int i=0;i < n;i++) { //Start Equation 3
				if(i != j) { //Stops multiplying QsubJ if i=j
					QsubJ = (z[j]-z[i])*QsubJ;
				}
			} //End Equation 3	
			fz = 1; //Evaluate f(zsubj) by using Horner's Method
			for(int k = n-1;k >= 0;k--) {
				fz = fz*z[j] + cList[k];
			}

			deltaZ[j] = (-fz/QsubJ);

			if(cabs(deltaZ[j]) > deltaZMax) {
				deltaZMax = cabs(deltaZ[j]);
			}
		} //End Equation 2
		
		for(int j=0;j < n;j++) { //Update zsubj
			z[j] = z[j] + deltaZ[j];
		}

		if(deltaZMax <= epsilon) { //If deltaZMax <= epsilon
			break;   //breaks out of loop
		}
		
	}

}
//===========================================================================
double findThetaj(int n,int j) { //Equation 7
	double theta = 0;
	theta = j*((2*M_PI)/n);
	return(theta);
}
