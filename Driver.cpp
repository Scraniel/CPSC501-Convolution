#include "WavData.h"
#include "Convolver.h"

int main(int argc, char * argv[])
{
	if(argc != 4)
	{
		std::cerr << "Usage: 'convolve <input file> <IR file> <output file>'\n";
		return -1;
	}

	WavData drySound, impulseResponse;
	drySound.loadWaveFile(argv[1]);
	impulseResponse.loadWaveFile(argv[2]);

	char* filename = argv[3];
	FILE *outputFileStream = fopen(filename, "wb");
	if (outputFileStream == NULL) {
		fprintf(stderr, "File %s cannot be opened for writing\n", filename);
			return -1;
	}


	std::cout << "Running tests...\n";
	if(!Convolver::RunTests())
		return -1;
	std::cout << "Tests Passed! Continuing program.\n";


	//WavData * convolved = Convolver::TimeDomainConvolve(drySound, impulseResponse);
	WavData * convolved = Convolver::FFTConvolve(drySound, impulseResponse);

	convolved->writeWaveFile(outputFileStream);

	fclose(outputFileStream);

}
