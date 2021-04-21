#include "xdock_header.h"

TData::TData() {
	ncount = 0;
	desc_real_for = nullptr;
	desc_comp_back = nullptr;
}

TData::~TData() {
	if (desc_comp_back != nullptr) {
		DftiFreeDescriptor(&desc_comp_back); desc_comp_back = nullptr;
	}
	if (desc_real_for != nullptr) {
		DftiFreeDescriptor(&desc_real_for); desc_real_for = nullptr;
	}
}

int TData::setVariables() {
	mt = mx * mx * mx;

	mrt = mrx * mrx * mrx;

	ntIn = nxIn * nyIn * nzIn;
	ntOut = nxOut * nyOut * nzOut;

	srcVolume.width = mx;
	srcVolume.height = mx;
	srcVolume.depth = mx;

	dstVolume.width = mrx;
	dstVolume.height = mrx;
	dstVolume.depth = mrx;

	srcVoi.x = 0;
	srcVoi.y = 0;
	srcVoi.z = 0;
	srcVoi.width = mx;
	srcVoi.height = mx;
	srcVoi.depth = mx;

	srcStep = 4 * mx;
	srcPlaneStep = 4 * mx * mx;
	mapStep = 4 * mrx;
	mapPlaneStep = 4 * mrx * mrx;
	dstStep = 4 * nxIn;
	dstPlaneStep = 4 * nxIn * nyIn;

	rotShift[0] = float(iRC_mask[0]);
	rotShift[1] = float(iRC_mask[1]);
	rotShift[2] = float(iRC_mask[2]);

	isExternal = true;
	return 0;
}

int TData::normaliseMaps() {
	ippsSum_32f(mM, ntIn, &const3, IppHintAlgorithm::ippAlgHintAccurate);
	ippsDivC_32f_I(const3, mM, ntIn);
	ippsSum_32f(mS, ntIn, &const6, IppHintAlgorithm::ippAlgHintAccurate);
	ippsAddProductC_32f(mM, -const6, mS, ntIn);
	ippsSum_32f(mS, ntIn, &const6, IppHintAlgorithm::ippAlgHintAccurate);
	ippsDotProd_32f(mS, mS, ntIn, &const2);
	const2 = 1.0f/sqrt(const2*const3);
	ippsMulC_32f_I(const2, mS, ntIn);
	

	status_fft = DftiComputeForward(desc_real_for, mS, mcS);
	status_fft = DftiComputeForward(desc_real_for, mM, mcM);

	real2compFFT(mcS, nxIn, nyIn, nzIn);
	real2compFFT(mcM, nxIn, nyIn, nzIn);

	ippsMul_32fc_I(mcShift, mcS, ntIn);
	ippsMul_32fc_I(mcShift, mcM, ntIn);

	ippsConj_32fc_I(mcS, ntIn);
	ippsConj_32fc_I(mcM, ntIn);


	/*
	Ipp32f* v1, * v2;
	Ipp32fc * vc;

	vc = ippsMalloc_32fc(ntIn);
	v1 = ippsMalloc_32f(ntIn);
	v2 = ippsMalloc_32f(ntIn);

	status_fft = DftiComputeBackward(desc_comp_back, mcS, vc);
	ippsReal_32fc(vc, v1, ntIn);
	
	char outputFileName[1000];
	float morig[3];
	morig[0] = 0.0f;
	morig[1] = 0.0f;
	morig[2] = 0.0f;
	FILE* fo;
	fo = fopen("C:\\Temp2\\xDock\\0420\\out\\fout.txt", "w");
	sprintf(outputFileName, "C:\\Temp2\\xDock\\0420\\out\\_mapS.mrc");
	writeMRC(fo, outputFileName, v1, nxIn, nyIn, nzIn, morig, pixT);
	fclose(fo);

	ippsFree(vc); vc = nullptr;
	ippsFree(v1); v1 = nullptr;
	ippsFree(v2); v2 = nullptr;
	*/
	return 0;
}

int TData::findRStd() {
	ippsMul_32fc(mcM, mcTc, mcRes, ntIn);
	//ippsConj_32fc_I(mcRes, ntIn);
	if (nScale > 1) readuceFmap(mcRes, nxOut, nyOut, nzOut, nScale);
	status_fft = DftiComputeBackward(desc_comp_back, mcRes, mc1);
	ippsSqr_32f_I((Ipp32f *)mc1, 2*ntOut);
		
	ippsMul_32fc(mcM, mcT2c, mcRes, ntIn);
	//ippsConj_32fc_I(mcRes, ntIn);
	if (nScale > 1) readuceFmap(mcRes, nxOut, nyOut, nzOut, nScale);
	status_fft = DftiComputeBackward(desc_comp_back, mcRes, mcRStd);
	
	ippsSub_32f_I((Ipp32f *)mc1, (Ipp32f *)mcRStd, 2*ntOut);
	ippsThreshold_LTVal_32f_I((Ipp32f *)mcRStd, 2*ntOut, stdCut2, fbig*fbig);
	ippsSqrt_32f_I((Ipp32f *)mcRStd, 2*ntOut);
	ippsDivCRev_32f_I(1.0f, (Ipp32f*)mcRStd, 2 * ntOut);
	
	return 0;
}



int TData::findCC() {
	ippsMul_32fc(mcS, mcTc, mcRes, ntIn);
	//ippsConj_32fc_I(mcRes, ntIn);
	if (nScale > 1) readuceFmap(mcRes, nxOut, nyOut, nzOut, nScale);
	status_fft = DftiComputeBackward(desc_comp_back, mcRes, mcCC);
	ippsMul_32f_I((Ipp32f *)mcRStd, (Ipp32f *)mcCC, 2*ntOut);
	return 0;
}


int TData::runFFTStd() {
	ippsMul_32fc(mcM, mcTc, mcRes, ntIn);
	//ippsConj_32fc_I(mcRes, ntIn);
	if (nScale > 1) readuceFmap(mcRes, nxOut, nyOut, nzOut, nScale);
	status_fft = DftiComputeBackward(desc_comp_back, mcRes, mc1);
	ippsSqr_32f_I((Ipp32f*)mc1, 2 * ntOut);

	ippsMul_32fc(mcM, mcT2c, mcRes, ntIn);
	//ippsConj_32fc_I(mcRes, ntIn);
	if (nScale > 1) readuceFmap(mcRes, nxOut, nyOut, nzOut, nScale);
	status_fft = DftiComputeBackward(desc_comp_back, mcRes, mcRStd);

	ippsSub_32f_I((Ipp32f*)mc1, (Ipp32f*)mcRStd, 2 * ntOut);
	ippsMax_32f((Ipp32f *)mcRStd, 2*ntOut, &maxV);
	return 0;
}


int TData::findBest() {
	Ipp32f f1, f2;
	Ipp32f* vCCin, *vCCout, *vQ;
	unsigned int u1, u2;
	u1 = uIndex;
	u2 = u1 + 1;



	/*Ipp32f* v1, * v2;
	Ipp32fc* vc;

	vc = ippsMalloc_32fc(ntOut);
	v1 = ippsMalloc_32f(ntOut);
	v2 = ippsMalloc_32f(ntOut);

	status_fft = DftiComputeBackward(desc_comp_back, mcTc, vc);
	ippsImag_32fc(mcCC, v1, ntOut);

	char outputFileName[1000];
	float morig[3];
	morig[0] = 0.0f;
	morig[1] = 0.0f;
	morig[2] = 0.0f;
	FILE* fo;
	fo = fopen("C:\\Temp2\\xDock\\0420\\out\\fout.txt", "w");
	sprintf(outputFileName, "C:\\Temp2\\xDock\\0420\\out\\_mapCC1.mrc");
	writeMRC(fo, outputFileName, v1, nxOut, nyOut, nzOut, morig, pixT);
	fclose(fo);

	ippsFree(vc); vc = nullptr;
	ippsFree(v1); v1 = nullptr;
	ippsFree(v2); v2 = nullptr;





	*/








	vCCin = (Ipp32f*)mcCC;
	vCCout = (Ipp32f*)vCCbest;
	vQ = (Ipp32f*)vQbest;

	*(unsigned int*)(&f1) = u1;
	*(unsigned int*)(&f2) = u2;

	for (int i = 0; i < 2*ntOut; i+=2) {
		if (vCCin[i] > vCCout[i]) {
			vCCout[i] = vCCin[i];
			vQ[i] = f1;
		}
		if (vCCin[i+1] > vCCout[i+1]) {
			vCCout[i+1] = vCCin[i+1];
			vQ[i+1] = f2;
		}
	}

	return 0;
}



int TData::findXYZmaps() {
	Ipp32f *v, *mR_loc;

	for (int i = 0; i < 3; i++) {
		if (i == 0) {
			v = mapX;
		}
		else if (i == 1) {
			v = mapY;
		}
		else {
			v = mapZ;
		}
		mR_loc = matRot +  i;

		ippsMulC_32f(mapXini, mR_loc[0], v, mrt);
		ippsMulC_32f(mapYini, mR_loc[3], map1, mrt);
		ippsAdd_32f_I(map1, v, mrt);
		ippsMulC_32f(mapZini, mR_loc[6], map1, mrt);
		ippsAdd_32f_I(map1, v, mrt);
		ippsAddC_32f_I(rotShift[i], v, mrt);
	}
	return 0;
}

int TData::rotateSM() {
	ippsSet_32f(-1.0f, mS, ntIn);

	ipprRemap_32f_C1V(mapMS, srcVolume, srcStep, srcPlaneStep, srcVoi, mapX, mapY, mapZ, mapStep, mapPlaneStep, mS, dstStep, dstPlaneStep, dstVolume, IPPI_INTER_NN);

	for (int i = 0; i < ntIn; i++) {
		if (mS[i] < -0.5) {
			mS[i] = 0.0f;
			mM[i] = 0.0f;
		}
		else {
			mM[i] = 1.0f;
		}
	}
		
	return 0;
}
