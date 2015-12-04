/*
 * Convolver.h
 *
 *  Created on: Dec 1, 2015
 *      Author: Danny Lewis
 */

#ifndef CONVOLVER_H_
#define CONVOLVER_H_

#include "WavData.h"
#include <climits>

class Convolver {
public:
	Convolver();

	static WavData * TimeDomainConvolve(WavData, WavData);
	static WavData * FFTConvolve(WavData, WavData);

private:
	static void TimeDomainConvolve(const double[], const int &, const double[], const int &, double[], const int &);
	static void FFTConvolve(double [], unsigned long, int);

	static double Normalize(const double &, const double &, const double &, const double &, const double &);
	static int NextHighestPowerOf2(int);
	static double * ZeroPadding(const short [], const int &);
};

#endif /* CONVOLVER_H_ */
