
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// war unbekannt
#define bool int
#define false 0

// Struktur, um das Ergebnis zurückzugeben
typedef struct Result {
	size_t length;
	double* times;
	double* values;
} Result;

double newtonSchritt(double x0, double f, double fAbleitung){
	return x0 - f/fAbleitung;
}

// yAbleitungsGleichung ist die rechte Seite der Gleichung y' = f(y,t)

// inline-Funktionsdefinitionen gehen in C leider nicht,
// deshalb, um Dinge zu vereinfachen, ist das Newton-Verfahren in dieser Funktion
// und die Fehlerfunktion ist als define definiert
#define MAX(a,b) (a > b ? a : b)
#define errorFunction(yi) (yAbleitungsGleichung(yi,tn) - (yi-yn)/h)

Result implizitesEulerVerfahren(
	double (*yAbleitungsGleichung)(double y, double t), 
	double y0, double t0, // start value and time
	double tEnd, // end time
	double h, // step size
	bool saveAllValues // whether the result shall contain all values
	){
	
	// berechne die benötigte Anzahl an Schritten, und wieviele Speicherstellen für
	// das Ergebnis benötigt werden
	size_t stepCount0 = (size_t) MAX(((tEnd-t0)/h), 0);
	size_t length  = saveAllValues ? stepCount0+1 : 1;
	double* times  = (double*) malloc(sizeof(double) * length);
	double* values = (double*) malloc(sizeof(double) * length);
	
	double maxError = 1e-3 * h;
	
	double yn = y0;
	times[0] = t0;
	values[0] = y0;
	
	// berechne alle Schritte
	for(size_t i=1;i<=stepCount0;i++){
		
		// berechne die n-te Zeit
		double tn = t0 + h * i;
		
		// Newton-Verfahren
		// max. 5000 Schritte, damit sich das Programm nicht
		// aufhängt falls es nicht konvergiert.
		int maxSteps = 5000, stepCount = 0;
		double yn_p1 = yn, error;
		do {
			
			// f ist die Errorfunktion...
			double f0 = errorFunction(yn_p1);
			double f1 = errorFunction(yn_p1+h);
			double dfdx = (f1-f0)/h;
			double nextValue = newtonSchritt(yn_p1, f0, dfdx);
			error = abs(nextValue-yn_p1);
			yn_p1 = nextValue;
			stepCount++;
			
		} while(error > maxError && stepCount < maxSteps);
		
		// speichere die Ergebnisse, falls wir die Zwischenwerte auch brauchen
		if(saveAllValues){
			times[i] = tn;
			values[i] = yn_p1;
		}
		
		yn = yn_p1;
		
	}
	
	// speichere das Ergebnis, wenn noch nicht getan
	if(!saveAllValues){
		times[0] = t0 + h * stepCount0;
		values[0] = yn;
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
	
	double targetValue = -log(1.0 + log(5.0));
	
	for(int i=1;i<=8;i++){
		double h = pow(10, -i);
		Result iev = implizitesEulerVerfahren(funcB, 0.0, -2.0, 0.0, h, false);
		// print it or export it
		double lastValue = iev.values[iev.length-1];
		printf("\n 10^-%d: %.17g,\ntarget: %.17g\n", i, lastValue, targetValue);
	}
	printf("\n");
	
}