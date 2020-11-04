
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct Result {
	size_t length;
	double* times;
	double* values;
} Result;

double newtonSchritt(double x0, double f, double fAbleitung){
	return x0 - f/fAbleitung;
}

#define MAX(a,b) (a > b ? a : b)
#define errorFunction(yi) (yAbleitungsGleichung(yi,tn) - (yi-yn)/h)

Result implizitesEulerVerfahren(double (*yAbleitungsGleichung)(double y, double t), double y0, double t0, double tEnd, double h){
	
	size_t length = (size_t) MAX(((tEnd-t0)/h+1), 1);
	double* times = (double*) malloc(sizeof(double) * length);
	double* values = (double*) malloc(sizeof(double) * length);
	
	double yn = y0;
	
	double maxError = 1e-3 * h;
	
	times[0] = t0;
	values[0] = y0;
	
	for(size_t i=1;i<length;i++){
		
		double tn = t0 + h * i;
		
		// newtonSchritt
		int maxSteps = 5000, stepCount = 0;
		double yn_p1 = yn, error;
		do {
			
			// f ist die Errorfunktion...
			double f0 = errorFunction(yn_p1);
			double f1 = errorFunction(yn_p1+h);
			double dfdx = (f1-f0)/h;
			double xn_p1 = newtonSchritt(yn_p1, f0, dfdx);
			error = abs(xn_p1-yn_p1);
			yn_p1 = xn_p1;
			stepCount++;
			
		} while(error > maxError && stepCount < maxSteps);
		
		times[i] = tn;
		values[i] = yn_p1;
		
		yn = yn_p1;
		
	}
	
	Result result;
	result.length = length;
	result.times = times;
	result.values = values;
	return result;
	
}

double funcB(double y, double t){
	return 2*t*exp(y)/(1+t*t);
}

void main(){
	
	for(int i=1;i<=8;i++){
		double h = pow(10, -i);
		Result iev = implizitesEulerVerfahren(funcB, 0.0, -2.0, 0.0, h);
		// todo print it or export it
		printf("10^-%d: %.17g, target: %.17g\n", i, iev.values[iev.length-1], -log(1.0 + log(5.0)));
	}
	
}