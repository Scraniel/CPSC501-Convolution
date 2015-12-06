/*
 * Convolver.cpp
 *
 *  Created on: Dec 1, 2015
 *      Author: Danny Lewis
 */

#include "Convolver.h"
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr;

Convolver::Convolver() {
	// TODO Auto-generated constructor stub

}

WavData * Convolver::FFTConvolve(WavData drySound, WavData impulseResponse)
{
	WavData * convolved = new WavData();

	convolved->setChannels(drySound.getChannels());
	convolved->setNumSamples(drySound.getNumberOfSamples() + impulseResponse.getNumberOfSamples() - 1);
	convolved->setBitsPerSample(drySound.getBitsPerSample());
	convolved->setSampleRate(drySound.getSampleRate());
	convolved->setData((short*) malloc(convolved->getNumberOfSamples() * convolved->getChannels() * (convolved->getBitsPerSample()/8)));


	int newsize = NextHighestPowerOf2(convolved->getNumberOfSamples());
	newsize = newsize << 1;


	double * drySoundPadded = ZeroPadding(drySound.getData(), drySound.getNumberOfSamples(), newsize);

	double * impulseResponsePadded = ZeroPadding(impulseResponse.getData(), impulseResponse.getNumberOfSamples(), newsize);

	FFTConvolve(drySoundPadded-1, newsize/2, 1);
	FFTConvolve(impulseResponsePadded-1, newsize/2, 1);

	double * multiplied = ComplexMultiplication(drySoundPadded, impulseResponsePadded, newsize);

	FFTConvolve(multiplied-1, newsize/2, -1);

	// Scale and Find min/max
	double min = multiplied[0]/(double)newsize;
	double max = min;
	double current;
	for(int i =0; i < newsize; i++)
	{
		multiplied[i] /= (double)newsize;
		current = multiplied[i];

		if(current < min)
			min = current;
		if(current > max)
			max = current;
	}

	// Normalize the values between -1 and 1
	for(int i = 0; i < newsize; i++)
		multiplied[i] = Normalize(multiplied[i], min, max, -1, 1);

	// Scale the data back to short
	for(int i = 0; i < convolved->getNumberOfSamples()*2; i+=2)
		convolved->getData()[i/2] = rint(multiplied[i] * SHRT_MAX);

	delete[] drySoundPadded;
	delete[] impulseResponsePadded;
	delete[] multiplied;

	return convolved;
}

WavData * Convolver::TimeDomainConvolve(WavData drySound, WavData impulseResponse)
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

	TimeDomainConvolve(xn, drySound.getNumberOfSamples(), 	 // x[n]
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
		convolved->getData()[i] = rint(tempValues[i] * SHRT_MAX);

	delete[] tempValues;
	delete[] xn;
	delete[] hn;

	return convolved;
}


/*
 * Written by Dr. Leonard Manzara. Uses doubles so values don't wrap around.
 */
void Convolver::TimeDomainConvolve(const double x[], const int & N, const double h[], const int & M, double y[], const int & P)
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

/*
 * Found on Page 507 - 508 of Numerical Recipes in C: The art of Scientific Computing.
 *
 * Modified by Danny Lewis
 */
void Convolver::FFTConvolve(double data[], unsigned long nn, int isign)
{

    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

double Convolver::Normalize(const double & value, const double & oldMin, const double & oldMax, const double & newMin, const double & newMax)
{
	double a = (newMax - newMin) / (oldMax - oldMin);
	double b = newMax - a * oldMax;

	return a * value + b;
}

long Convolver::NextHighestPowerOf2(long number)
{
	number--;
	number |= number >> 1;
	number |= number >> 2;
	number |= number >> 4;
	number |= number >> 8;
	number |= number >> 16;
	number++;

	return number;
}

// returns a zero padded (double) array representing the original array as
// series of complex numbers. The length is guaranteed to be a power of 2
double * Convolver::ZeroPadding(const short toPad[], const int & size, const int & newsize)
{
	double * padded = new double [newsize];

	for(int i =0, k = 0 ; k < newsize; i++, k+=2)
	{
		if(i < size)
			padded[k] = toPad[i];
		else
			padded[k] = 0;
		padded[k+1] = 0;
	}


	return padded;
}

double * Convolver::ComplexMultiplication(const double first[], const double second[], const int & length)
{
	double * multiplied = new double[length];

	for(int i = 0; i < length; i+=2)
	{

		multiplied[i] = first[i] * second[i] - first[i+1] * second[i+1];
		multiplied[i+1] = first[i] * second[i+1] + first[i+1] * second[i];
	}

	return multiplied;
}
