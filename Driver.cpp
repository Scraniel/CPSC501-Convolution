#include "WavData.h"
#include "Convolver.h"

int main()
{
	WavData drySound, impulseResponse;
	drySound.loadWaveFile("DrySounds/DrumsDry.wav");
	impulseResponse.loadWaveFile("ImpulseResponses/BigHall.wav");


	std::cout << "Samples: " << impulseResponse.getNumberOfSamples() << std::endl;
	std::cout << "Channels: " << impulseResponse.getChannels() << std::endl;
	std::cout << "Bits / Sample: " << impulseResponse.getBitsPerSample() << std::endl;
	std::cout << "Sample Rate: " << impulseResponse.getSampleRate() << std::endl;

	char* filename = "test.wav";
	FILE *outputFileStream = fopen(filename, "wb");
	if (outputFileStream == NULL) {
		fprintf(stderr, "File %s cannot be opened for writing\n", filename);
			return -1;
	}

	WavData * convolved = Convolver::TimeDomainConvolve(drySound, impulseResponse);

	convolved->writeWaveFile(outputFileStream);

	fclose(outputFileStream);

}
