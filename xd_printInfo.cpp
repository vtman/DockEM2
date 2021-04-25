#include "xdock_header.h"

int myImin(int a, int b) {
	if (a < b) {
		return a;
	}
	return b;
}

int myImax(int a, int b) {
	if (a < b) {
		return b;
	}
	return a;
}


float myFmin(float a, float b) {
	if (a < b) {
		return a;
	}
	return b;
}

float myFmax(float a, float b) {
	if (a < b) {
		return b;
	}
	return a;
}

double myDmin(double a, double b) {
	if (a < b) {
		return a;
	}
	return b;
}

double myDmax(double a, double b) {
	if (a < b) {
		return b;
	}
	return a;
}



int XData::printParameters() {
	int max_threads;

	fprintf(flog, "Input parameters:\n\n");

	fprintf(flog, "Target file: \"%s\"\n", ip->inputTargetFile);
	fprintf(flog, "Search file: \"%s\"\n", ip->inputSearchFile);
	fprintf(flog, "Output folder: \"%s\"\n", ip->outputFolder);
	if (strlen(ip->systemFolder) > 0) {
		fprintf(flog, "System folder: \"%s\"\n", ip->systemFolder);
	}
	else {
		fprintf(flog, "System folder: \"\"\n");
	}
	fprintf(flog, "Prefix for output files: \"%s\"\n", ip->outputPrefix);

#pragma omp parallel
	{
		max_threads = omp_get_num_threads();
	}

	fprintf(flog, "Number of threads (available): %i\n", max_threads);
	fprintf(flog, "Number of threads: %i\n", ip->numberOfThreads);
	ip->numberOfThreads = myImin(max_threads, myImax(1, ip->numberOfThreads));
	fprintf(flog, "Number of threads (corrected): %i\n", ip->numberOfThreads);
	
	fprintf(flog, "Angular resolution: %f\n", ip->angularResolution);
	ip->angularResolution = myFmin(30.0f, myFmax(2.0, ip->angularResolution));
	fprintf(flog, "Angular resolution (corrected): %f\n", ip->angularResolution);

	fprintf(flog, "Angle to estimate Std: %f\n", ip->angleForEstimation);
	ip->angleForEstimation = myFmin(30.0, myFmax(2.0, ip->angleForEstimation));
	fprintf(flog, "Angle to estimate Std (corrected): %f\n", ip->angleForEstimation);

	fprintf(flog, "ROI centre: %i, %i, %i\n", ip->iCenterX, ip->iCenterY, ip->iCenterZ);
	fprintf(flog, "ROI size: %i x %i x %i\n", ip->boxSizeX, ip->boxSizeY, ip->boxSizeZ);

	ip->boxSizeX = myImin(200, myImax(16, ip->boxSizeX));
	ip->boxSizeY = myImin(200, myImax(16, ip->boxSizeY));
	ip->boxSizeZ = myImin(200, myImax(16, ip->boxSizeZ));

	fprintf(flog, "ROI size (corrected): %i x %i x %i\n", ip->boxSizeX, ip->boxSizeY, ip->boxSizeZ);

	ip->nShapeROI = SHAPE_DIFF;
	ip->pairsInGroup = 2;
	if (ip->boxSizeX == ip->boxSizeY && ip->boxSizeX == ip->boxSizeZ) {
		ip->nShapeROI = SHAPE_XYZ;
		ip->pairsInGroup = 12;
	}
	else if (ip->boxSizeX == ip->boxSizeY) {
		ip->nShapeROI = SHAPE_XY;
		ip->pairsInGroup = 4;
	}
	else if (ip->boxSizeX == ip->boxSizeZ) {
		ip->nShapeROI = SHAPE_XZ;
		ip->pairsInGroup = 4;
	}
	else if (ip->boxSizeY == ip->boxSizeZ) {
		ip->nShapeROI = SHAPE_YZ;
		ip->pairsInGroup = 4;
	}

	switch (ip->nShapeROI) {
	case SHAPE_DIFF: fprintf(flog, "Shape: different\n"); break;
	case SHAPE_XY: fprintf(flog, "Shape: X = Y\n"); break;
	case SHAPE_XZ: fprintf(flog, "Shape: X = Z\n"); break;
	case SHAPE_YZ: fprintf(flog, "Shape: Y = Z\n"); break;
	case SHAPE_XYZ: fprintf(flog, "Shape: X = Y = Z\n"); break;
	}

	fprintf(flog, "Minimum signal: %f\n", ip->signalLevel);
	ip->signalLevel = myFmin(0.99f, myFmax(0.0001f, ip->signalLevel));
	fprintf(flog, "Minimum signal (corrected): %f\n", ip->signalLevel);

	fprintf(flog, "Resolution for search map: %f\n", ip->searchResolution);
	ip->searchResolution = myFmin(10.0f, myFmax(0.5f, ip->searchResolution));
	fprintf(flog, "Resolution for search map (corrected): %f\n", ip->searchResolution);

	fprintf(flog, "Mask level: %f\n", ip->maskLevel);
	ip->maskLevel = myFmin(0.9f, myFmax(0.01f, ip->maskLevel));
	fprintf(flog, "Mask level (corrected): %f\n", ip->maskLevel);
	
	fprintf(flog, "Radius for the mask (in A): %f\n", ip->maskRadius);
	ip->maskRadius = myFmin(50.0f, myFmax(0.1f, ip->maskRadius));
	fprintf(flog, "Radius for the mask (in A, corrected): %f\n", ip->maskRadius);

	fprintf(flog, "Pixel size for search/mask (in A): %f\n", ip->maskPixel);
	ip->maskPixel = myFmin(5.0f, myFmax(0.1f, ip->maskPixel));
	fprintf(flog, "Pixel size for search/mask (in A, corrected): %f\n", ip->maskPixel);

	fprintf(flog, "Output reduction: %i\n", ip->outputReduction);
	ip->outputReduction = myImin(3, myImax(1, ip->outputReduction));
	fprintf(flog, "Output reduction (corrected): %i\n", ip->outputReduction);

	fprintf(flog, "Number of output pdb files: %i\n", ip->nPDBout);

	if (ip->writeMask) {
		fprintf(flog, "Write mask to file: yes\n");
	}
	else {
		fprintf(flog, "Write mask to file: no\n");
	}

	if (ip->writeSearch) {
		fprintf(flog, "Write search map to file: yes\n");
	}
	else {
		fprintf(flog, "Write search map to file: no\n");
	}

	if (ip->writeTarget) {
		fprintf(flog, "Write target map to file: yes\n");
	}
	else {
		fprintf(flog, "Write target map to file: no\n");
	}

	if (ip->writeROI) {
		fprintf(flog, "Write ROI map to file: yes\n");
	}
	else {
		fprintf(flog, "Write ROI map to file: no\n");
	}

	if (ip->sphericalMask) {
		fprintf(flog, "Spherical mask is used: yes\n");
	}
	else {
		fprintf(flog, "Spherical mask is used: no\n");
	}

	fprintf(flog, "Minimum CC value for a peak: %f\n", ip->minCCpeak);
	fprintf(flog, "Number of peaks (extra info): %i\n", ip->nRelativeInfo);
	fprintf(flog, "\n");

	return 0;
}
