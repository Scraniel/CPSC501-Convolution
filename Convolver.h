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
	static bool RunTests();

private:
	static void TimeDomainConvolve(const double[], const int &, const double[], const int &, double[], const int &);
	static void FFTConvolve(double [], unsigned long, int);

	static double Normalize(const double &, const double &, const double &, const double &, const double &);
	static long NextHighestPowerOf2(long);
	static double * ZeroPadding(const short [], const int &, const int &);
	static double * ComplexMultiplication(const double [], const double [], const int & );
	static void FindMinMaxAndScale(double [], const int &, double &, double &, const double &);

	/*** TESTS ***/
	static bool test_Normalize();
	static bool test_NextHighestPowerOf2();
	static bool test_ZeroPadding();
	static bool test_ComplexMultiplication();
	static bool test_TimeDomainConvolve();
	static bool test_FFTConvolve();
};

#endif /* CONVOLVER_H_ */
