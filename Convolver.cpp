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


	for(int i =0; i < drySound.getNumberOfSamples(); i++)
	{
		for(int j = 0; j < impulseResponse.getNumberOfSamples(); j++)
		{
			convolved->getData()[i + j] += drySound.getData()[i] * impulseResponse.getData()[j];
		}
	}

	return convolved;
}
