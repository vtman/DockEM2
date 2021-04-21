#include "xdock_header.h"

InputParameters::InputParameters() {
	numberOfThreads = 1;
	maskPixel = 0.5f;

	writeMask = false;
	writeROI = false;
	writeSearch = false;
	writeTarget = false;
	
	minCCpeak = 0.1f;
	nRelativeInfo = 10;
	nPDBout = 20;
	angularResolution = 6.0f;
	angleForEstimation = 12.0f;
	maskPixel = 0.5f;
}

int InputParameters::readInputParameters(int argc, char** argv) {
	int m;
	int isSearch, isOutput, isPrefix, isTarget, isCentre, isBoxSize;
	int isSignal, isMaskLevel, isRadius, isResolution;

	bool t;

	char p1[1000];
	
	isResolution = 0;
	isTarget = 0;
	isSearch = 0;
	isOutput = 0;
	isPrefix = 0;
	isCentre = 0;
	isBoxSize = 0;
	isSignal = 0;
	isMaskLevel = 0;
	isRadius = 0;

	outputReduction = 1;
	sprintf(systemFolder, "\0");
	sprintf(inputTargetFile, "\0");
	sprintf(inputSearchFile, "\0");
	sprintf(outputFolder, "\0");
	sprintf(outputPrefix, "\0");
		

	printf("Number of input parameters: %i\n\n", argc);

	if (argc < 20) {
		printUsage();
		return -1;
	}

	for (int i = 0; i < argc; i++) {

		printf("# %i: \"%s\"\n", i, argv[i]);
		m = strlen(argv[i]);

		if (i % 2 == 1) {
			t = false;
			for (int j = 0; j < m; j++) {
				p1[j] = tolower(argv[i][j]);
			}
			p1[m] = '\0';

			if (p1[0] != '-') {
				printf("Error: wrong parameter\n");
				return -1;
			}

			if (i + 1 == argc) {
				printf("Error: wrong number of parameters\n");
				return -3;
			}

			if (strcmp(p1, "-target") == 0) {
				isTarget = 1;
				sprintf(inputTargetFile, argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-search") == 0) {
				isSearch = 1;
				sprintf(inputSearchFile, argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-outputfolder") == 0) {
				isOutput = 1;
				sprintf(outputFolder, argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-systemfolder") == 0) {
				sprintf(systemFolder, argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-prefix") == 0) {
				isPrefix = 1;
				sprintf(outputPrefix, argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-writemask") == 0) {
				if (argv[i + 1][0] == 'y' || argv[i + 1][0] == 'Y') writeMask = true;
				t = true;
			}

			if (strcmp(p1, "-writesearch") == 0) {
				if (argv[i + 1][0] == 'y' || argv[i + 1][0] == 'Y') writeSearch = true;
				t = true;
			}

			if (strcmp(p1, "-writeroi") == 0) {
				if (argv[i + 1][0] == 'y' || argv[i + 1][0] == 'Y') writeROI = true;
				t = true;
			}

			if (strcmp(p1, "-writetarget") == 0) {
				if (argv[i + 1][0] == 'y' || argv[i + 1][0] == 'Y') writeTarget = true;
				t = true;
			}

			if (strcmp(p1, "-minccpeak") == 0) {
				minCCpeak = (float)atof(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-nrelativeinfo") == 0) {
				nRelativeInfo = atoi(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-npdbout") == 0) {
				nPDBout = atoi(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-outputreduction") == 0) {
				outputReduction = atoi(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-threads") == 0) {
				numberOfThreads = atoi(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-angularresolution") == 0) {
				angularResolution = (float)atof(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-angleforestimation") == 0) {
				angleForEstimation = (float)atof(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-boxsize") == 0) {
				isBoxSize = 1;
				if (sscanf(argv[i + 1], "%i,%i,%i", &boxSizeX, &boxSizeY, &boxSizeZ) != 3) {
					if (sscanf(argv[i + 1], "%i", &boxSizeX) != 1) {
						printf("Error: cannot read 3 numbers");
						return -7;
					}
					else {
						boxSizeY = boxSizeX;
						boxSizeZ = boxSizeX;
					}
				}
				t = true;
			}

			if (strcmp(p1, "-resolution") == 0) {
				isResolution = 1;
				searchResolution = (float)atof(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-signallevel") == 0) {
				isSignal = 1;
				signalLevel = (float)atof(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-masklevel") == 0) {
				isMaskLevel = 1;
				maskLevel = (float)atof(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-maskradius") == 0) {
				isRadius = 1;
				maskRadius = (float)atof(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-maskpixel") == 0) {
				maskPixel = (float)atof(argv[i + 1]);
				t = true;
			}

			if (strcmp(p1, "-center") == 0) {
				isCentre = 1;
				if (sscanf(argv[i + 1], "%i,%i,%i", &iCenterX, &iCenterY, &iCenterZ) != 3) {
					printf("Error: cannot read 3 numbers");
					return -7;
				}
				t = true;
			}
						
			if (!t) {
				printf("Error: unknown parameter\n");
				return -2;
			}

		}

	}

	printf("\n");

	t = false;
	
	if (isTarget == 0) {
		printf("Error: unknown target file (-target)\n");
		t = true;
	}

	if (isSearch == 0) {
		printf("Error: unknown search file (-search)\n");
		t = true;
	}

	if (isOutput == 0) {
		printf("Error: unknown output folder (-outputFolder)\n");
		t = true;
	}

	if (isPrefix == 0) {
		printf("Error: unknown output prefix (-prefix)\n");
		t = true;
	}

	if (isBoxSize == 0) {
		printf("Error: unknown size of the ROI box (-boxSize)\n");
		t = true;
	}

	if (isCentre == 0) {
		printf("Error: unknown position of ROI box (-center)\n");
		t = true;
	}
	
	if (isSignal == 0) {
		printf("Error: unknown min signal level (-signalLevel, [0.0, 1.0])\n");
		t = true;
	}

	if (isResolution == 0) {
		printf("Error: unknown resolution (-resolution, [0.5, 10.0])\n");
		t = true;
	}

	if (isMaskLevel == 0) {
		printf("Error: unknown mask level (-maskLevel, [0.0, 1.0])\n");
		t = true;
	}
	
	if (isRadius == 0) {
		printf("Error: unknown radius for the mask (-maskRadius)\n");
		t = true;
	}
	
	if (t) return -1;

	return 0;
}

int InputParameters::printUsage() {
	
	printf("Usage (s - string, f - real, i - integer):\n");
	printf("\t -target s               Target file in MRC format\n");
	printf("\t -search s               Search file in PDB\n");
	printf("\t -outputFolder s         Output folder\n");
	printf("\t -systemFolder s         System folder\n");
	printf("\t -prefix s               Prefix for output files\n");
	printf("\t -threads i              Number of threads\n");
	printf("\t -angularResolution f    Angular resolution (in degrees, >= 2)\n");
	printf("\t -angleForEstimation f   Angle to estimate max Std value\n");
	printf("\t -centre i,i,i           Centre of the ROI box\n");
	printf("\t -boxSize i,i,i or i     Size of the ROI box\n");
	printf("\t -signalLevel f          Minimum signal level [0.0, 1.0]\n");
	printf("\t -resolution f           Search map resolution (in A, [0.5, 10.0])\n");
	printf("\t -maskLevel f            Mask level [0.0, 1.0]\n");
	printf("\t -maskRadius f           Radius to create mask (in A)\n");
	printf("\t -maskPixel f            Pixels size for mask/search maps (in A)\n");
	printf("\t -outputRreduction i     Values for each m-th voxel (m = 1, 2, 3)\n");
	printf("\t -nPDBout i              Number of output PDB files (default: 20)\n");
	printf("\t -writeMask s            Write mask map to file (yes, NO)\n");
	printf("\t -writeSearch s          Write search map to file (yes, NO)\n");
	printf("\t -writeTarget s          Write target map [used for FFT] to file (yes, NO)\n");
	printf("\t -writeROI s             Write ROI of target map to file (yes, NO)\n");
	printf("\t -minCCpeak f            Min value of a peak to be reported (default: 0.1)\n");
	printf("\t -nRelativeInfo i        Extra information for the first n peaks (default: 10)\n");
	printf("\n");

	return 0;
}
