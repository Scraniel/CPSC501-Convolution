/*
 * Convolver.h
 *
 *  Created on: Dec 1, 2015
 *      Author: Danny Lewis
 */

#ifndef CONVOLVER_H_
#define CONVOLVER_H_

#include "WavData.h"

class Convolver {
public:
	Convolver();

	static WavData * Convolve(WavData, WavData);
};

#endif /* CONVOLVER_H_ */
