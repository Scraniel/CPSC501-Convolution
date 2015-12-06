/*
 * WavData.cpp
 *
 *  Created on: Dec 1, 2015
 *      Author: Danny Lewis
 */
#include "WavData.h"

void WavData::loadWaveFile(char *fileName)
{

	FILE* fp = fopen(fileName,"rb");
	if (fp) {
		char id[5];
		unsigned long fsize;
		fread(id, sizeof(unsigned char), 4, fp);
		id[4] = '\0';

		if (!strcmp(id, "RIFF")) {
			fread(&fsize, sizeof(unsigned char)*4, 1, fp);
			fread(id, sizeof(unsigned char), 4, fp);
			id[4] = '\0';

			if (!strcmp(id,"WAVE")) {
				fread(id, sizeof(unsigned char), 4, fp);
				fread(&formatLength, sizeof(unsigned char)*4,1,fp);
				fread(&formatTag, sizeof(short), 1, fp);
				fread(&channels, sizeof(short),1,fp);
				fread(&sampleRate, sizeof(unsigned char)*4, 1, fp);
				fread(&avgBytesSec, sizeof(unsigned char)*4, 1, fp);
				fread(&blockAlign, sizeof(short), 1, fp);
				fread(&bitsPerSample, sizeof(short), 1, fp);
				fread(id, sizeof(unsigned char), 4, fp);
				if(id[0] != 'd')
					fread(id, sizeof(unsigned char), 2, fp); // sometimes there are 2 extra 0's at the end
				fread(&dataSize, sizeof(unsigned char)*4, 1, fp);

				samples = dataSize/(channels * (bitsPerSample/8));
				data = (short*) malloc(dataSize);
				fread(data, sizeof(short), samples, fp);
				}
			else {
				std::cout << "Error: RIFF file but not a wave file\n";
				}
			}
		else {
			std::cout << "Error: not a RIFF file\n";
			}
		}
	fclose(fp);
}

unsigned long WavData::getNumberOfSamples()
{
	return samples;
}

unsigned long WavData::getSampleRate()
{
	return sampleRate;
}

short WavData::getBitsPerSample()
{
	return bitsPerSample;
}

short WavData::getChannels()
{
	return channels;
}

short* WavData::getData()
{
	return data;
}

void WavData::setData(short* data)
{
	this->data = data;
}
void WavData::setChannels(short channels)
{
	this->channels = channels;
}

void WavData::setNumSamples(unsigned long samples)
{
	this->samples = samples;
}

void WavData::setSampleRate(unsigned long sampleRate)
{
	this->sampleRate = sampleRate;
}

void WavData::setBitsPerSample(short bitsPerSample)
{
	this->bitsPerSample = bitsPerSample;
}

/******************************************************************************
*		Author: 		Dr. Leonard Manzara
*		Modified By: 	Daniel Lewis
*
*       function:       writeWaveFileHeader
*
*       purpose:        Writes the header in WAVE format to the output file.
*
*       arguments:      channels:  the number of sound output channels
*                       numberSamples:  the number of sound samples
*                       outputRate:  the sample rate
*                       outputFile:  the output file stream to write to
*
*       internal
*       functions:      fwriteIntLSB, fwriteShortLSB
*
*       library
*       functions:      ceil, fputs
*
******************************************************************************/

void WavData::writeWaveFile(FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * samples * (bitsPerSample/8);

    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;

    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * (bitsPerSample/8);

    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(sampleRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);

    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);

    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);

    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)sampleRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(bitsPerSample, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);

	for (int i = 0; i < samples; i++) {

		/*  Convert the value to a 16-bit integer, with the
			range -maximumValue to + maximumValue.  The calculated
			value is rounded to the nearest integer  */
		short int sampleValue = data[i];

		/*  Write out the sample as a 16-bit (short) integer
			in little-endian format  */
		fwriteShortLSB(sampleValue, outputFile);
	}
}

/******************************************************************************
*		Author: 		Dr. Leonard Manzara
*
*       function:       fwriteIntLSB
*
*       purpose:        Writes a 4-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t WavData::fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}



/******************************************************************************
*		Author: 		Dr. Leonard Manzara
*
*       function:       fwriteShortLSB
*
*       purpose:        Writes a 2-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t WavData::fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}
