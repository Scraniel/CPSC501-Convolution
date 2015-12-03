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

	static WavData * Convolve(WavData, WavData);

private:
	static void Convolve(const double[], int, const double[], int, double[], int);
	static double Normalize(double, double, double, double, double);
};

#endif /* CONVOLVER_H_ */
