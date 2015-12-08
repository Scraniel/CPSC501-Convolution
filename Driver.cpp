#include "WavData.h"
#include "Convolver.h"

int main()
{
	WavData drySound, impulseResponse;
	drySound.loadWaveFile("DrySounds/DrumsDry.wav");
	impulseResponse.loadWaveFile("ImpulseResponses/BigHall.wav");

	std::cout << "Samples: " << drySound.getNumberOfSamples() << std::endl;
	std::cout << "Channels: " << drySound.getChannels() << std::endl;
	std::cout << "Bits / Sample: " << drySound.getBitsPerSample() << std::endl;
	std::cout << "Sample Rate: " << drySound.getSampleRate() << std::endl;

	char* filename = "TDCtest.wav";
	FILE *outputFileStream = fopen(filename, "wb");
	if (outputFileStream == NULL) {
		fprintf(stderr, "File %s cannot be opened for writing\n", filename);
			return -1;
	}


	std::cout << "Running tests...\n";
	if(!Convolver::RunTests())
		return 1;
	std::cout << "Tests Passed! Continuing program.\n";


	//WavData * convolved = Convolver::TimeDomainConvolve(drySound, impulseResponse);
	//WavData * convolved = Convolver::FFTConvolve(drySound, impulseResponse);

	//convolved->writeWaveFile(outputFileStream);

	fclose(outputFileStream);

}
