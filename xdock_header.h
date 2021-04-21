#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <string.h>
#include <math.h>
//#include <mathimf.h>
#include "gaussConsts.h"

#include "ipp.h"
#include "mkl.h"

#include <omp.h>

#define nCurves 512
#define curvePrecision 0.001

#define pi180 0.01745329251994329576923690768489
#define pi180r 57.295779513082320876798154814105

#define pi2 6.2831853071795

#define pi 3.141592653589

#define fbig 1e10

#define SHAPE_DIFF 0
#define SHAPE_XY 1
#define SHAPE_XZ 2
#define SHAPE_YZ 3
#define SHAPE_XYZ 4

struct quat {
	double q[4];
};

class MRCinfo {
public:
	MRCinfo() {};
	~MRCinfo() {};

	bool isSwap;
	int nx, ny, nz;
	int nx_start, ny_start, nz_start;
	int iSpaceGroupNumber, sizeOfExtendedHeader;
	int nVersion, numberOfLabels, nt;

	float densityMin, densityMax, densityMean;
	float cellA[3], cellB[3], rms, origin[3];
};

class InputParameters {
public:
	InputParameters();
	~InputParameters() {};

	int readInputParameters(int argc, char** argv);
	int printUsage();

	char inputTargetFile[1000], inputSearchFile[1000], outputFolder[1000], outputPrefix[100], systemFolder[1000];
	int numberOfThreads, boxSizeX, boxSizeY, boxSizeZ, iCenterX, iCenterY, iCenterZ;
	int outputReduction;
	float signalLevel, maskPixel, searchResolution, angularResolution, angleForEstimation;
	float maskRadius, maskLevel;

	int nShapeROI;

	int pairsInGroup;
	bool writeMask, writeSearch, writeTarget, writeROI;
	float minCCpeak;
	int nRelativeInfo, nPDBout;
};

class ListAtoms {
public:
	ListAtoms();
	~ListAtoms();

	int findMaxStep(double xx, double yy, double zz, double *veca, double R, double* astep);
	double distBall(double xx, double yy, double zz, int *ind);
	int createCurves(FILE *flog, float resolution);
	Ipp64f valueAB(Ipp64f x, int ind);
	Ipp64f valueXYZ(Ipp64f x, int ind, int ig);
	Ipp32f getValue(Ipp32f xx, int ind, int ig);
	int findMap(Ipp32f *vm, int nx, int ny, int nz, Ipp32f pixelSize, int icx, int icy, int icz, Ipp32f cx, Ipp32f cy, Ipp32f cz);

	int nAtoms, countAtom[100], numDifAtoms, nGauss;
	int *index, *indexAtom, *indexD;
	float *vx, *vy, *vz;
	Ipp32f **vGValues, *vRmax, *vdR, *vVmax, **vGxyz;
	Ipp32f *vT, *vX, *vY, *vZ;
	Ipp32f RadZero, aR, res;
};

class TData {
public:
	TData();
	~TData();

	int findXYZmaps();
	int rotateSM();
	int setVariables();
	int normaliseMaps();
	int runFFTStd();
	int findRStd();
	int findCC();
	
	int findBest();
	

	Ipp32f matRot[9];
	
	int pairsInGroup;

	bool isExternal;

	Ipp32u uIndex;
	int nScale;
	int nxIn, nyIn, nzIn, ntIn;
	int nxOut, nyOut, nzOut, ntOut;
	int mx, mt;
	int mrx, mrt;

	float pixM, pixT, rotShift[3];
	int iRC_mask[3];// , iRC_global[3];
	
	Ipp32f const3, const6, const2;

	Ipp32f maxMTstd, maxV, stdCut2;

	MKL_LONG status_fft;
	DFTI_DESCRIPTOR_HANDLE desc_real_for, desc_comp_back, desc_comp_back2;
	
	IpprVolume srcVolume, dstVolume;
	IpprCuboid srcVoi;
	int srcStep, srcPlaneStep, mapStep, mapPlaneStep, dstStep, dstPlaneStep;

	Ipp32fc *vCCbest, *vQbest;
	Ipp32f *mapX, *mapY, *mapZ;
	Ipp32f *mapXini, *mapYini, *mapZini;

	Ipp32fc **pmcT, **pmcT2, **pmcCCbest, **pmcQbest;

	Ipp32fc *mcTc, *mcT2c, *mcRes, *mc1;
	Ipp32fc* mcCC, *mcRStd, *mcShift;

	Ipp32f *mM, *mS, *map1;
	Ipp32fc *mcM, *mcS;
	Ipp32f *mapMS;

	int ncount;

};

class XData {
public:
	XData(InputParameters *ip_in);
	~XData();

	quat* vLocalQuat, *vInputQuat;

	int nTemp;
	int nGroupsAvailable, nGroupMain, nGroupEstimate;
	int nCornerShift;
	double* vAngleLevel;

	Ipp64f *matRot;

	InputParameters *ip;
	ListAtoms *la;

	MRCinfo *minfo;

	FILE *flog;

	char outputFileName[1000];
	char *pdbData;
	int pdbRecords, pO1, pO2, pO3;

	Ipp32f *mapT, *mapT2, *mapTemp;

	bool isExternal;

	float cmx, cmy, cmz;
	float cgx, cgy, cgz;
	float pixT, pixM, pixTOut;
	float targetOrig[3], maskOrig[3];

	int nxIn, nyIn, nzIn, ntIn;
	int nxOut, nyOut, nzOut, ntOut;
	int mx, mt;
	int mbx, mby, mbz, mbt;
	int ntx, nty, ntz, ntt;
	int nux, nuy, nuz, nut;

	int irotCentMask[3], irotCentGlobal[3], iCornerOrig[3], iCornerGlobal[3];
	int iMaskLen;

	Ipp32f maxStdValue;

	Ipp32fc **pmcT, **pmcT2, ** pmcCCbest, **pmcQbest, **pmcShift;
	Ipp32fc** pmcTemp;
	Ipp32fc** pmcRStd, ** pmcCC;
	Ipp32fc** pmcS, ** pmcM;

	Ipp32f ** pmXini, **pmYini, **pmZini;
	Ipp32f  **pmMS;
	Ipp32f *mMask, *mSearch;
	Ipp32f* vOut;

	Ipp32u* uOut;

	TData **td;

	float borderSize;

	int npeaks;

	int nRotations, nRot_final, nAngle;

	MKL_LONG status_fft, dims_in[3], dims_out[3], strides_in[4], strides_out[4];
	DFTI_DESCRIPTOR_HANDLE desc_real_for, desc_comp_for, desc_comp_back, desc_real_back, desc_comp_back2, desc_comp_for_inplace;
	float scale_for, scale_back;
	
	int correctMaskCentre();
	int getNumberOfGroups();
	int findPeaks();
	int writeOutput();
	int cropOutputData();
	int startProcessing();
	int reduceOutputData();
	int coreProcessing();
	int coreStdProcessing();
	int findSearchMap();
	int findMaskMap();
	int findRotationMatrices();
	int checkParameters();
	int printParameters();
	int readAtoms();
	int readPDB();
	int findCentreOfMass();
	int findGeomCentre();
	int writeRotatedPDB(char *fPDB, float cx, float cy, float cz, float shiftX, float shiftY, float shiftZ, quat qres);
	int readTargetMap();
	int cropTarget();
	int setFFT();
	int rotateTarget();
	int prepareMultiThread();
	int allocateThreadMemory();
};

int mapQMatrix(Ipp32f* vIn, Ipp32f* vOut, Ipp32f *v1, quat* q, int nx, int ny, int nz, int ind);

int getAtomIndex(FILE *fout, char *s, int *ind5);

int writeMRC(FILE *fl, char *fileName, Ipp32f *v, int nx, int ny, int nz, float *orig, float pixelSize);
int getMRCinfo(FILE *fl, char *fileName, MRCinfo *minfo);
int real2compFFT(Ipp32fc* vC, int nx, int ny, int nz);

int readuceFmap(Ipp32fc* vf, int nx, int ny, int nz, int nScale);
int findBestAll(Ipp32fc* vCCold, Ipp32fc* uQold, Ipp32fc* vCCnew, Ipp32fc* uQnew, int nt);

int quat2mat(quat a, Ipp64f* m);
quat quatProd(quat a, quat b);
double findAngle(quat a);
quat quatBetween(quat a, quat b);
double distQ2(quat a, quat b);

int myImin(int a, int b);
int myImax(int a, int b);
float myFmin(float a, float b);
float myFmax(float a, float b);
double myDmin(double a, double b);
double myDmax(double a, double b);

