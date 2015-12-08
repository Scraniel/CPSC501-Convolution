/*
 * WavData.h
 *
 *  Created on: Dec 1, 2015
 *      Author: Danny Lewis
 */

#ifndef WAVDATA_H_
#define WAVDATA_H_

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

class WavData {
	public:
        WavData() {
            data = 0;
            samples = 0;
        }

        void loadWaveFile(char*);
        void writeWaveFile(FILE*);

        unsigned long getNumberOfSamples();
        unsigned long getSampleRate();
        short getChannels();
        short getBitsPerSample();

        short* getData();
        void setData(short*);
        void setChannels(short);
        void setNumSamples(unsigned long);
        void setSampleRate(unsigned long);
        void setBitsPerSample(short);

        bool equals(WavData);
	private:
        short* data;
        unsigned long samples;
        short formatTag, channels, blockAlign, bitsPerSample;
        unsigned long formatLength, sampleRate, avgBytesSec, dataSize;

        size_t fwriteIntLSB(int, FILE*);
        size_t fwriteShortLSB(short int, FILE*);


};

#endif /* WAVDATA_H_ */
