/*
 * Convolver.cpp
 *
 *  Created on: Dec 1, 2015
 *      Author: Danny Lewis
 */

#include "Convolver.h"

Convolver::Convolver() {
	// TODO Auto-generated constructor stub

}

WavData * Convolver::Convolve(WavData drySound, WavData impulseResponse)
{
	WavData * convolved = new WavData();


	convolved->setChannels(drySound.getChannels());
	convolved->setNumSamples(drySound.getNumberOfSamples() + impulseResponse.getNumberOfSamples() - 1);
	convolved->setBitsPerSample(drySound.getBitsPerSample());
	convolved->setSampleRate(drySound.getSampleRate());
	convolved->setData((short*) malloc(convolved->getNumberOfSamples() * convolved->getChannels() * (convolved->getBitsPerSample()/8)));

	for(int i = 0; i < convolved->getNumberOfSamples(); i++)
		convolved->getData()[i] = 0;

	double * tempValues = new double[convolved->getNumberOfSamples()];
	double * xn = new double[drySound.getNumberOfSamples()];
	double * hn = new double[impulseResponse.getNumberOfSamples()];


	// Add data to arrays
	for(int i = 0; i < drySound.getNumberOfSamples(); i++)
		xn[i] = drySound.getData()[i];
	for(int i =0; i < impulseResponse.getNumberOfSamples(); i++)
		hn[i] = impulseResponse.getData()[i];

	Convolve(xn, drySound.getNumberOfSamples(), 			 // x[n]
			hn, impulseResponse.getNumberOfSamples(), 		 // h[n]
			tempValues, convolved->getNumberOfSamples()); 	 // y[n]

	// Find min and max
	double min = tempValues[0];
	double max = min;
	double current;
	for(int i =0; i < convolved->getNumberOfSamples(); i++)
	{
		current = tempValues[i];

		if(current < min)
			min = current;
		if(current > max)
			max = current;
	}

	// Normalize the values between -1 and 1
	for(int i = 0; i < convolved->getNumberOfSamples(); i++)
		tempValues[i] = Normalize(tempValues[i], min, max, -1, 1);

	// Scale the data back to short
	for(int i = 0; i < convolved->getNumberOfSamples(); i++)
		convolved->getData()[i] = tempValues[i] * SHRT_MAX;

	delete[] tempValues;
	delete[] xn;
	delete[] hn;

	return convolved;
}


/*
 * Written by Dr. Leonard Manzara. Uses doubles so values don't wrap around.
 */
void Convolver::Convolve(const double x[], int N, const double h[], int M, double y[], int P)
{
	int n, m;

	for (n = 0; n < P; n++)
		y[n] = 0.0;


	for (n = 0; n < N; n++) {
		for (m = 0; m < M; m++) {
			y[n+m] += (double) (x[n] * h[m]);
		}
	}
}

double Convolver::Normalize(double value, double oldMin, double oldMax, double newMin, double newMax)
{
	double a = (newMax - newMin) / (oldMax - oldMin);
	double b = newMax - a * oldMax;

	return a * value + b;
}
