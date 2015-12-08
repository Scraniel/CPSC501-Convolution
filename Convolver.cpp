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

void Convolver::FindMinMaxAndScale(double array[], const int & length, double & min, double & max, const double & scale)
{
	min = array[0]/scale;
	max = min;
	double current;
	int i;
	for(i =0; i < length-2; i+=4)
	{
		array[i] /= scale;
		current = array[i];

		if(current < min)
			min = current;
		if(current > max)
			max = current;

		array[i+2] /= scale;
		current = array[i+2];

		if(current < min)
			min = current;
		if(current > max)
			max = current;
	}
	if(i == length - 2)
	{
		array[i] /= scale;
		current = array[i];

		if(current < min)
			min = current;
		if(current > max)
			max = current;
	}

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

	FFTConvolve(drySoundPadded-1, newsize >> 1, 1);
	FFTConvolve(impulseResponsePadded-1, newsize >> 1, 1);

	double * multiplied = ComplexMultiplication(drySoundPadded, impulseResponsePadded, newsize);

	FFTConvolve(multiplied-1, newsize >> 1, -1);

	// Scale and Find min/max
	double min, max;
	FindMinMaxAndScale(multiplied, newsize, min, max, drySound.getNumberOfSamples());


	// Jamming optimization & partial unrolling
	// Normalize the values between -1 and 1
	int i;
	for(i = 0; i < newsize - 2; i+=4)
	{
		multiplied[i] = Normalize(multiplied[i], min, max, -1, 1);
		convolved->getData()[i >> 1] = rint(multiplied[i] * SHRT_MAX);

		multiplied[i+2] = Normalize(multiplied[i+2], min, max, -1, 1);
		convolved->getData()[(i+2) >> 1] = rint(multiplied[i+2] * SHRT_MAX);
	}
	if(i == newsize - 2)
	{
		multiplied[i] = Normalize(multiplied[i], min, max, -1, 1);
		convolved->getData()[i >> 1] = rint(multiplied[i] * SHRT_MAX);
	}

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

    for (i = 1; i < n-2; i += 4) {
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

		if (j > i+2) {
			SWAP(data[j], data[i+2]);
			SWAP(data[j+1], data[i+3]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
    }
    if(i == n-2)
    {
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
	if(oldMax == oldMin)
		return oldMax;
	if(newMin == newMax)
		return newMax;

	double a = (newMax - newMin) / (oldMax - oldMin);
	double b = newMax - a * oldMax;

	return a * value + b;
}

long Convolver::NextHighestPowerOf2(long number)
{
	if(!number)
		return 2;

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
// series of complex numbers. The new length must be at least 2 * the original length.
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

/***********************************************************************************/
/*** TEST DEFINITIONS - These tests can be run by invoking Convolver::RunTests() ***/
/***********************************************************************************/


bool Convolver::RunTests()
{
	bool result = true;

	result &= test_Normalize();
	result &= test_NextHighestPowerOf2();
	result &= test_ZeroPadding();
	result &= test_ComplexMultiplication();
	result &= test_TimeDomainConvolve();
	result &= test_FFTConvolve();

	if(!result)
		std::cerr << "One or more tests failed. Check the console to see details." << std::endl;

	return result;
}

bool Convolver::test_Normalize()
{

	bool passed = true;
	double value, min, max, newMin, newMax, result, expected;

	// Edge case: All Zero
	value = 0;
	min = 0;
	max = 0;
	newMin = 0;
	newMax = 0;
	expected = 0;
	result = Normalize(value, min, max, newMin, newMax);

	if(result != expected)
	{
		std::cerr << "test_Normalize() failed on test 'All Zero'. Result: " << result << std::endl;
		passed = false;
	}

	// Typical case: moving a number between -1 and 1
	value = 75;
	min = 0;
	max = 100;
	newMin = -1;
	newMax = 1;
	expected = 0.5;
	result = Normalize(value, min, max, newMin, newMax);

	if(result != expected)
	{
		std::cerr << "test_Normalize() failed on test 'Typical Case: -1 -> 1'. Result: " << result << std::endl;
		passed = false;
	}

	return passed;
}

bool Convolver::test_NextHighestPowerOf2()
{
	bool passed = true;
	long value, result, expected;

	// Edge case: 0
	value = 0;
	expected = 2;
	result = NextHighestPowerOf2(value);

	if(result != expected)
	{
		std::cerr << "test_NextHighestPowerOf2() failed on test 'Edge case: 0'. Result: " << result << std::endl;
		passed = false;
	}

	// Typical case: high number
	value = 10000000;
	expected = 16777216;
	result = NextHighestPowerOf2(value);

	if(result != expected)
	{
		std::cerr << "test_NextHighestPowerOf2() failed on test 'Typical Case: high number'. Result: " << result << std::endl;
		passed = false;
	}

	return passed;
}

bool Convolver::test_ZeroPadding()
{

	bool passed = true;
	double * result;
	int size, newSize;

	// Edge Case: newLength = 2*length
	size = 5;
	newSize = 10;
	short value1[] = {1,2,3,4,5};
	double expected1[] = {1,0,2,0,3,0,4,0,5,0};
	result = ZeroPadding(value1, size, newSize);

	for(int i = 0; i < newSize; i++)
	{
		if(result[i] != expected1[i])
		{
			std::cerr << "test_ZeroPadding failed on test 'Edge Case: newLength = 2*length'. Result: {";
			for(int i = 0; i < newSize; i++)
			{
				std::cerr << result[i];
				if(i < newSize -1)
					std::cerr << ", ";
			}
			std::cerr << "}\n";

			passed = false;
			break;
		}
	}
	delete[] result;

	// Typical case: newLength is a power of 2, 0's at end
	size = 5;
	newSize = 16;
	short value2[] = {1,2,3,4,5};
	double expected2[] = {1,0,2,0,3,0,4,0,5,0,0,0,0,0,0,0};

	result = ZeroPadding(value2, size, newSize);

	for(int i = 0; i < newSize; i++)
	{
		if(result[i] != expected2[i])
		{
			std::cerr << "test_ZeroPadding failed on test 'Typical case: newLength is a power of 2, 0's at end'. Result: {";
			for(int i = 0; i < newSize; i++)
			{
				std::cerr << result[i];
				if(i < newSize -1)
					std::cerr << ", ";
			}
			std::cerr << "}\n";

			passed = false;
			break;
		}
	}
	delete[] result;

	return passed;
}

bool Convolver::test_ComplexMultiplication()
{
	bool passed = true;

	int length;
	double * result;

	// Edge Case: One number each
	length = 2;
	double value1[] = {3,4};
	double value2[] = {7,1};
	double expected1[] = {17,31};

	result = ComplexMultiplication(value1, value2, length);

	for(int i = 0; i < length; i++)
	{
		if(result[i] != expected1[i])
		{
			std::cerr << "test_ComplexMuliplication failed on test 'Edge Case: One number each'. Result: {";
			for(int i = 0; i < length; i++)
			{
				std::cerr << result[i];
				if(i < length -1)
					std::cerr << ", ";
			}
			std::cerr << "}\n";

			passed = false;
			break;
		}
	}
	delete[] result;

	// Typical Case: multiple numbers, some imaginary 0s and negatives
	length = 6;
	double value3[] = {3,4,5,6,7,8};
	double value4[] = {7,1,8,0,9,-1};
	double expected2[] = {17,31,40,48,71,65};

	result = ComplexMultiplication(value3, value4, length);

	for(int i = 0; i < length; i++)
	{
		if(result[i] != expected2[i])
		{
			std::cerr << "test_ComplexMuliplication failed on test 'Typical Case: multiple numbers, some imaginary 0s and negatives'. Result: {";
			for(int i = 0; i < length; i++)
			{
				std::cerr << result[i];
				if(i < length -1)
					std::cerr << ", ";
			}
			std::cerr << "}\n";

			passed = false;
			break;
		}
	}
	delete [] result;


	return passed;
}

bool Convolver::test_TimeDomainConvolve()
{
	bool passed = true;

	WavData drySound, impulseResponse, expected;
	WavData * result;


	short xn[] = {3,4,5};
	short hn[] = {2,1};
	short yn[] = {6,11,14,5}; // in a perfect world, this would be the answer. However,
							  // since we scale the numbers by MAX_SHRT, the values get
							  // skewed up heavily. This being said, the intermediate
							  // value (just after the convolution) IS this exactly.
	short ynr[] = {-25485, 10922, 32767, -32767};

	drySound.setData(xn);
	drySound.setChannels(1);
	drySound.setNumSamples(3);
	drySound.setSampleRate(44100);
	drySound.setBitsPerSample(16);

	impulseResponse.setData(hn);
	impulseResponse.setChannels(1);
	impulseResponse.setNumSamples(2);
	impulseResponse.setSampleRate(44100);
	impulseResponse.setBitsPerSample(16);

	expected.setData(ynr);
	expected.setChannels(1);
	expected.setNumSamples(4);
	expected.setSampleRate(44100);
	expected.setBitsPerSample(16);

	result = TimeDomainConvolve(drySound, impulseResponse);

	if(!result->equals(expected))
	{
		std::cerr << "test_TimeDomainConvolve() failed on test 'Typical Case'. Result: {";
		for(int i = 0; i < result->getNumberOfSamples(); i++)
		{
			std::cerr << result->getData()[i];
			if(i < result->getNumberOfSamples() -1)
				std::cerr << ", ";
		}
		std::cerr << "}\n";
		passed = false;
	}

	delete result;

	return passed;
}

bool Convolver::test_FFTConvolve()
{
	bool passed = true;

		WavData drySound, impulseResponse, expected;
		WavData * result;


		short xn[] = {3,4,5};
		short hn[] = {2,1};
		short yn[] = {6,11,14,5}; // in a perfect world, this would be the answer. However,
								  // since we scale the numbers by MAX_SHRT, the values get
								  // skewed up heavily. This being said, the intermediate
								  // value (just after the convolution) IS this exactly.
		short ynr[] = {-25485, 10922, 32767, -32767};

		drySound.setData(xn);
		drySound.setChannels(1);
		drySound.setNumSamples(3);
		drySound.setSampleRate(44100);
		drySound.setBitsPerSample(16);

		impulseResponse.setData(hn);
		impulseResponse.setChannels(1);
		impulseResponse.setNumSamples(2);
		impulseResponse.setSampleRate(44100);
		impulseResponse.setBitsPerSample(16);

		expected.setData(ynr);
		expected.setChannels(1);
		expected.setNumSamples(4);
		expected.setSampleRate(44100);
		expected.setBitsPerSample(16);

		result = FFTConvolve(drySound, impulseResponse);

		if(!result->equals(expected))
		{
			std::cerr << "test_FFTConvolve() failed on test 'Typical Case'. Result: {";
			for(int i = 0; i < result->getNumberOfSamples(); i++)
			{
				std::cerr << result->getData()[i];
				if(i < result->getNumberOfSamples() -1)
					std::cerr << ", ";
			}
			std::cerr << "}\n";
			passed = false;
		}

		delete result;

		return passed;
}

