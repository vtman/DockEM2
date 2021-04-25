#include "xdock_header.h"

XData::XData(InputParameters *ip_in) {
	vAngleLevel = nullptr;
	vLocalQuat = nullptr;
	vInputQuat = nullptr;

	desc_real_for = nullptr;
	desc_real_back = nullptr;
	desc_comp_for = nullptr;
	desc_comp_back = nullptr;
	desc_comp_for_inplace = nullptr;

	pdbData = nullptr;
	
	minfo = nullptr;

	ip = nullptr;
	la = nullptr;

	mapT = nullptr;
	mapT2 = nullptr;
	matRot = nullptr;
	
	mapTemp = nullptr;

	minfo = nullptr;

	pmXini = nullptr;
	pmYini = nullptr;
	pmZini = nullptr;
	pmcT = nullptr;
	pmcT2 = nullptr;
	pmcCCbest = nullptr;
	pmcQbest = nullptr;
	mMask = nullptr;
	mSearch = nullptr;
	pmMS = nullptr;

	pmcS = nullptr;
	pmcM = nullptr;
	
	pmcShift = nullptr;
	pmcCC = nullptr;
	pmcRStd = nullptr;
	pmcTemp = nullptr;

	vOut = nullptr;
	uOut = nullptr;

	td = nullptr;

	ip = ip_in;

	la = new ListAtoms();

	flog = nullptr;
	nTemp = 3;
}

XData::~XData() {
	if (vAngleLevel != nullptr) { free(vAngleLevel); vAngleLevel = nullptr; }
	if (vLocalQuat != nullptr) { free(vLocalQuat); vLocalQuat = nullptr; }
	if (vInputQuat != nullptr) { free(vInputQuat); vInputQuat = nullptr; }
	if (pdbData != nullptr) { free(pdbData); pdbData = nullptr; }
	if (minfo != nullptr) { delete minfo; minfo = nullptr; }
	if (desc_comp_for != nullptr) {
		DftiFreeDescriptor(&desc_comp_for); desc_comp_for = nullptr;
	}
	if (desc_comp_back != nullptr) {
		DftiFreeDescriptor(&desc_comp_back); desc_comp_back = nullptr;
	}
	if (desc_real_for != nullptr) {
		DftiFreeDescriptor(&desc_real_for); desc_real_for = nullptr;
	}
	if (desc_comp_for_inplace != nullptr) {
		DftiFreeDescriptor(&desc_comp_for_inplace); desc_comp_for_inplace = nullptr;
	}

	if (la != nullptr) { delete la; la = nullptr; }

	if (flog != nullptr) { fclose(flog); flog = nullptr; }

	if (mapTemp != nullptr) { ippsFree(mapTemp); mapTemp = nullptr; }
	if (mapT != nullptr) { ippsFree(mapT); mapT = nullptr; }
	if (mapT2 != nullptr) { ippsFree(mapT2); mapT2 = nullptr; }

	if (matRot != nullptr) { ippsFree(matRot); matRot = nullptr; }

	if (vOut != nullptr) { ippsFree(vOut); vOut = nullptr; }
	if (uOut != nullptr) { ippsFree(uOut); uOut = nullptr; }

	if (pmcShift != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			if (pmcShift[i] != nullptr) {
				ippsFree(pmcShift[i]); pmcShift[i] = nullptr;
			}
		}
		free(pmcShift); pmcShift = nullptr;
	}

	if (pmcT != nullptr) {
		for (int i = 0; i < ip->pairsInGroup*ip->numberOfThreads; i++) {
			if (pmcT[i] != nullptr) { ippsFree(pmcT[i]); pmcT[i] = nullptr; }

		}
		free(pmcT); pmcT = nullptr;
	}
	
	if (pmcT2 != nullptr) {
		for (int i = 0; i < ip->pairsInGroup*ip->numberOfThreads; i++) {
			if (pmcT2[i] != nullptr) { ippsFree(pmcT2[i]); pmcT2[i] = nullptr; }
		}
		free(pmcT2); pmcT2 = nullptr;
	}

	if (pmcCCbest != nullptr) {
		for (int i = 0; i < ip->pairsInGroup * ip->numberOfThreads; i++) {
			if (pmcCCbest[i] != nullptr) { ippsFree(pmcCCbest[i]); pmcCCbest[i] = nullptr; }
		}
		free(pmcCCbest); pmcCCbest = nullptr;
	}

	if (pmcQbest != nullptr) {
		for (int i = 0; i < ip->pairsInGroup * ip->numberOfThreads; i++) {
			if (pmcQbest[i] != nullptr) { ippsFree(pmcQbest[i]); pmcQbest[i] = nullptr; }
		}
		free(pmcQbest); pmcQbest = nullptr;
	}
	

	if (pmXini != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			if (pmXini[i] != nullptr) { ippsFree(pmXini[i]); pmXini[i] = nullptr; }
		}
		free(pmXini); pmXini = nullptr;
	}

	if (pmYini != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			if (pmYini[i] != nullptr) { ippsFree(pmYini[i]); pmYini[i] = nullptr; }
		}
		free(pmYini); pmYini = nullptr;
	}

	if (pmZini != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			if (pmZini[i] != nullptr) { ippsFree(pmZini[i]); pmZini[i] = nullptr; }
		}
		free(pmZini); pmZini = nullptr;
	}

	if (pmcCC != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			if (pmcCC[i] != nullptr) { ippsFree(pmcCC[i]); pmcCC[i] = nullptr; }
		}
		free(pmcCC); pmcCC = nullptr;
	}

	if (pmcRStd != nullptr) {
		for (int i = 0; i < ip->pairsInGroup * ip->numberOfThreads; i++) {
			if (pmcRStd[i] != nullptr) { ippsFree(pmcRStd[i]); pmcRStd[i] = nullptr; }
		}
		free(pmcRStd); pmcRStd = nullptr;
	}
	
	if (mMask != nullptr) { ippsFree(mMask); mMask = nullptr; }
	if (mSearch != nullptr) { ippsFree(mSearch); mSearch = nullptr; }

	if (pmMS != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			if (pmMS[i] != nullptr) { ippsFree(pmMS[i]); pmMS[i] = nullptr; }
		}
		free(pmMS); pmMS = nullptr;
	}

	if (pmcS != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			if (pmcS[i] != nullptr) { ippsFree(pmcS[i]); pmcS[i] = nullptr; }
		}
		free(pmcS); pmcS = nullptr;
	}

	if (pmcM != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			if (pmcM[i] != nullptr) { ippsFree(pmcM[i]); pmcM[i] = nullptr; }
		}
		free(pmcM); pmcM = nullptr;
	}

	if (pmcTemp != nullptr) {
		for (int i = 0; i < nTemp*ip->numberOfThreads; i++) {
			if (pmcTemp[i] != nullptr) { ippsFree(pmcTemp[i]); pmcTemp[i] = nullptr; }
		}
		free(pmcTemp); pmcTemp = nullptr;
	}

	if (td != nullptr) {
		for (int i = 0; i < ip->numberOfThreads; i++) {
			delete td[i]; td[i] = nullptr;
		}
		free(td); td = nullptr;
	}

	
}

int XData::checkParameters() {
	printf("Checking parameters\n");
	char logFile[1000];
	sprintf(logFile, "%s/%s.log", ip->outputFolder, ip->outputPrefix);

	flog = fopen(logFile, "w");
	if (flog == nullptr) {
		printf("Error: cannot open the log file \"%s\"\n", logFile);
		return -1;
	}

	printParameters();
	printf("... Done.\n");
	return 0;
}



int getAtomIndex(FILE *fout, char *s, int *ind5) {
	bool t;

	t = false;
	for (int i = 0; i < nTabAtoms5; i++) {
		if (strncmp(s, atomNames5 + 2 * i, 2) != 0) continue;
		*ind5 = i;
		t = true;
		break;
	}
	if (!t) {
		printf("Error: cannot find 5 gaussian data for \"%.2s\".\n", s);
		fprintf(fout, "Error: cannot find 5 gaussian data for \"%.2s\".\n", s);
		return -5;
	}

	return 0;
}




int XData::readAtoms(){
	printf("Reading atoms");
	int ires;
	ires = readPDB();
	printf("... Done.\n");

	return ires;
}

int XData::findCentreOfMass() {
	printf("Finding centre of mass");
	double sx, sy, sz, mass, dm;
	sx = 0.0;
	sy = 0.0;
	sz = 0.0;
	mass = 0.0;
	for (int i = 0; i < la->nAtoms; i++) {
		dm = double(la->index[i] + 1);
		mass += dm;
		sx += dm * double(la->vx[i]);
		sy += dm * double(la->vy[i]);
		sz += dm * double(la->vz[i]);
	}

	sx /= mass;
	sy /= mass;
	sz /= mass;

	cmx = (float)sx;
	cmy = (float)sy;
	cmz = (float)sz;

	fprintf(flog, "Centre of mass: %f, %f, %f\n", cmx, cmy, cmz);

	printf("... Done.\n");

	return 0;
}

int XData::findGeomCentre() {
	printf("Finding geometrical centre");
	int indb, na;
	double gx, gy, gz, dr, astep, phi, theta;
	double veca[3], da;
	gx = (double)cmx;
	gy = (double)cmy;
	gz = (double)cmz;
	
	na = 10;
	da = pi / double(na - 1);
	
	for (int ii = 0; ii < 1000; ii++) {
		phi = da * double(ii % (2 * na));
		theta = da * double((ii / (2 * na))%na);

		veca[0] = cos(phi)*sin(theta);
		veca[1] = sin(phi)*sin(theta);
		veca[2] = cos(theta);

		dr = la->distBall(gx, gy, gz, &indb);

		la->findMaxStep(gx, gy, gz, veca, dr, &astep);

		astep *= 0.5;
		gx += astep * veca[0];
		gy += astep * veca[1];
		gz += astep * veca[2];
	}

	dr = la->distBall(gx, gy, gz, &indb);

	cgx = (float)gx;
	cgy = (float)gy;
	cgz = (float)gz;

	fprintf(flog, "Geometric centre: %f, %f, %f\n", cgx, cgy, cgz);
	fprintf(flog, "Radius: %f\n", dr);

	printf("... Done.\n");

	/*cgx = cmx;
	cgy = cmy;
	cgz = cmz;
	*/
	return 0;
}

int XData::findRotationMatrices() {
	printf("Finding matrices for rotation");

	for (int i = 0; i < nGroupsAvailable; i++) {
		quat2mat(vInputQuat[i], matRot + i * 9);
	}

	printf("... Done.\n");
	return 0;
}

int XData::findSearchMap() {
	printf("Forming the search map");
	int ind1;
	float radius, rx, ry, rz;
	Ipp32f vmax;
	Ipp32f *mapS;

	la->createCurves(flog, ip->searchResolution);

	radius = (float)(la->distBall(cgx, cgy, cgz, &ind1));
	mx = 2*((int)((radius+ la->RadZero) / ip->maskPixel) + 2);
	mx = ((mx / 4) + 1) * 4;
	mt = mx * mx * mx;

	irotCentMask[0] = mx / 2;
	irotCentMask[1] = irotCentMask[0];
	irotCentMask[2] = irotCentMask[0];

	fprintf(flog, "\nSearch/mask box: %i x %i x %i\n", mx, mx, mx);
	fprintf(flog, "Pixel size (in A): %f\n", ip->maskPixel);

	pixM = ip->maskPixel;

	rx = ip->maskPixel * (mx / 2);
	ry = rx;
	rz = rx;

	maskOrig[0] = cgx - rx;
	maskOrig[1] = cgy - ry;
	maskOrig[2] = cgz - rz;

	fprintf(flog, "Mask origin (in A): %f, %f, %f\n", maskOrig[0], maskOrig[1], maskOrig[2]);

	
	mSearch = ippsMalloc_32f(mt);

	if (mSearch == nullptr) {
		fprintf(flog, "Error: cannot allocate memory for mapS.\n");
		printf("Error: cannot allocate memory for mapS.\n");
		return -1;
	}
		
	mapS = mSearch;

	ippsZero_32f(mapS, mt);
	
	la->findMap(mapS, mx, mx, mx, pixM, irotCentMask[0], irotCentMask[1], irotCentMask[2], cgx, cgy, cgz);

	ippsMax_32f(mapS, mt, &vmax);
	ippsDivC_32f_I(vmax, mapS, mt);

	if (ip->writeSearch) {
		sprintf(outputFileName, "%s/%s_mapS.mrc", ip->outputFolder, ip->outputPrefix);
		if (writeMRC(flog, outputFileName, mapS, mx, mx, mx, maskOrig, pixM) != 0) return -1;
	}

	printf("... Done.\n");
	return 0;
}




int XData::readTargetMap() {
	printf("Reading the target map");
	FILE *ft;
	Ipp32f vmax;
	minfo = new MRCinfo();
	if (getMRCinfo(flog, ip->inputTargetFile, minfo) != 0) return -1;

	ntx = minfo->nx;
	nty = minfo->ny;
	ntz = minfo->nz;
	ntt = ntx * nty * ntz;

	mapTemp = ippsMalloc_32f(ntt);
	if (mapTemp == nullptr) {
		fprintf(flog, "Error: cannot allocate memeory for mapTemp\n");
		printf("Error: cannot allocate memeory for mapTemp\n");
		return -1;
	}

	ft = fopen(ip->inputTargetFile, "rb");
	fseek(ft, 1024 + minfo->sizeOfExtendedHeader, SEEK_SET);

	fread(mapTemp, sizeof(float), ntt, ft);

	fclose(ft); ft = nullptr;

	ippsMaxAbs_32f(mapTemp, ntt, &vmax);
	
	if (vmax > 1e10) ippsSwapBytes_32u_I((Ipp32u *)mapTemp, ntt);

	pixT = minfo->cellA[0] / float(minfo->nx);

	fprintf(flog, "Pixel size, target (in A): %f\n", pixT);

	fprintf(flog, "Done.\n");

	printf("... Done.\n");

	return 0;
}

int XData::cropTarget() {
	printf("Croping the target map");
	int iTopOrig[3], iTopGlobal[3];
	int k1, k2, k3, i1, i2, i3, j1, j2, j3, nScale;
	float oT[3];


	switch (ip->outputReduction) {
	case 1: nScale = 4; break;
	case 2: nScale = 4; break;
	case 3: nScale = 12; break;
	case 4: nScale = 4; break;
	case 5: nScale = 20; break;
	}


	iCornerOrig[0] = ip->iCenterX - ip->boxSizeX / 2;
	iCornerOrig[1] = ip->iCenterY - ip->boxSizeY / 2;
	iCornerOrig[2] = ip->iCenterZ - ip->boxSizeZ / 2;

	iTopOrig[0] = iCornerOrig[0] + ip->boxSizeX - 1;
	iTopOrig[1] = iCornerOrig[1] + ip->boxSizeY - 1;
	iTopOrig[2] = iCornerOrig[2] + ip->boxSizeZ - 1;

	fprintf(flog, "ROI box: [%i, %i] x [%i, %i] x [%i, %i]\n", iCornerOrig[0], iTopOrig[0], iCornerOrig[1], iTopOrig[1], iCornerOrig[2], iTopOrig[2]);

	iMaskLen = (int)(borderSize / pixT) + 2;
	
	fprintf(flog, "Mask half-width (in pixels, Target): %i\n", iMaskLen);

	mbx = 2 * (iMaskLen + 1);

	mby = mbx;
	mbz = mbx;
	mbt = mbx * mby * mbz;

	fprintf(flog, "Mask/search box (in pixels, target): %i x %i x %i\n", mbx, mby, mbz);

	irotCentGlobal[0] = iMaskLen;
	irotCentGlobal[1] = iMaskLen;
	irotCentGlobal[2] = iMaskLen;

	fprintf(flog, "Rotation centre (in pixels, target): %i, %i, %i\n", irotCentGlobal[0], irotCentGlobal[1], irotCentGlobal[2]);
	fprintf(flog, "Rotation centre (in pixels, search): %i, %i, %i\n", irotCentMask[0], irotCentMask[1], irotCentMask[2]);

	nCornerShift = iMaskLen / nScale;
	if (iMaskLen % nScale > 0) nCornerShift++;
	nCornerShift *= nScale;

	iCornerGlobal[0] = iCornerOrig[0] - nCornerShift;
	iCornerGlobal[1] = iCornerOrig[1] - nCornerShift;
	iCornerGlobal[2] = iCornerOrig[2] - nCornerShift;

	iTopGlobal[0] = iTopOrig[0] + nCornerShift;
	iTopGlobal[1] = iTopOrig[1] + nCornerShift;
	iTopGlobal[2] = iTopOrig[2] + nCornerShift;

	nxIn = iTopGlobal[0] - iCornerGlobal[0] + 1;
	nyIn = iTopGlobal[1] - iCornerGlobal[1] + 1;
	nzIn = iTopGlobal[2] - iCornerGlobal[2] + 1;

	int nss;
	nss = ip->outputReduction;

	nux = iTopOrig[0] - iCornerOrig[0] + 1;
	if (nux % nss > 0) nux += nss;
	nux /= nss;
	nuy = iTopOrig[1] - iCornerOrig[1] + 1;
	if (nuy % nss > 0) nuy += nss;
	nuy /= nss;
	nuz = iTopOrig[2] - iCornerOrig[2] + 1;
	if (nuz % nss > 0) nuz += nss;
	nuz /= nss;
	nut = nux * nuy * nuz;

	nxOut = nxIn / nScale; 
	if (nxIn % nScale > 0) nxOut++;
	nyOut = nyIn / nScale;
	if (nyIn % nScale > 0) nyOut++;
	nzOut = nzIn / nScale;
	if (nzIn % nScale > 0) nzOut++;

	
	
	nxIn = nScale * nxOut;
	nyIn = nScale * nyOut;
	nzIn = nScale * nzOut;
	ntIn = nxIn * nyIn * nzIn;


	nxOut = nxIn / ip->outputReduction;
	nyOut = nyIn / ip->outputReduction;
	nzOut = nzIn / ip->outputReduction;

	ntOut = nxOut * nyOut * nzOut;

	iTopGlobal[0] = iCornerGlobal[0] + nxIn - 1;
	iTopGlobal[1] = iCornerGlobal[1] + nyIn - 1;
	iTopGlobal[2] = iCornerGlobal[2] + nzIn - 1;

	fprintf(flog, "Box for processing: %i x %i x %i\n", nxIn, nyIn, nzIn);
	fprintf(flog, "Box position: [%i, %i] x [%i, %i] x [%i, %i]\n", iCornerGlobal[0], iTopGlobal[0], iCornerGlobal[1], iTopGlobal[1], iCornerGlobal[2], iTopGlobal[2]);
	
	if (iCornerGlobal[0] < 0 || iCornerGlobal[1] < 0 || iCornerGlobal[2] < 0 || iTopGlobal[0] >= ntx || iTopGlobal[1] >= nty || iTopGlobal[2] >= ntz) {
		fprintf(flog, "Warning: extra padding");
	}

	mapT = ippsMalloc_32f(ntIn);
	mapT2 = ippsMalloc_32f(ntIn);

	if (mapT == nullptr || mapT2 == nullptr) {
		fprintf(flog, "Error: cannot allocate memeory for mapT.\n");
		printf("Error: cannot allocate memeory for mapT.\n");
		return -1;
	}

	ippsZero_32f(mapT, ntIn);
	k1 = myImin(ntz - 1, myImax(0, iCornerGlobal[2]));
	i1 = myImin(nty - 1, myImax(0, iCornerGlobal[1]));
	j1 = myImin(ntx - 1, myImax(0, iCornerGlobal[0]));

	k2 = myImin(ntz - 1, myImax(0, iTopGlobal[2]));
	i2 = myImin(nty - 1, myImax(0, iTopGlobal[1]));
	j2 = myImin(ntx - 1, myImax(0, iTopGlobal[0]));

	k3 = myImax(0, -iCornerGlobal[2]) - k1;
	i3 = myImax(0, -iCornerGlobal[1]) - i1;
	j3 = myImax(0, -iCornerGlobal[0]);

	fprintf(flog, "Valid box: [%i, %i] x [%i, %i] x [%i, %i]\n", j1, j2, i1, i2, k1, k2);

	for (int k = k1; k <= k2; k++) {
		for (int i = i1; i <= i2; i++) {
			ippsCopy_32f(mapTemp + (k*nty + i)*ntx + j1, mapT + ((k + k3)*nyIn + i + i3)*nxIn + j3, j2 - j1 + 1);
		}
	}

	oT[0] = minfo->origin[0] + float(iCornerGlobal[0])*pixT;
	oT[1] = minfo->origin[1] + float(iCornerGlobal[1])*pixT;
	oT[2] = minfo->origin[2] + float(iCornerGlobal[2])*pixT;

	if (ip->writeTarget) {
		sprintf(outputFileName, "%s/%s_mapT.mrc", ip->outputFolder, ip->outputPrefix);
		if (writeMRC(flog, outputFileName, mapT, nxIn, nyIn, nzIn, oT, pixT) != 0) return -2;
	}

	ippsFree(mapTemp); mapTemp = nullptr;

	vOut = ippsMalloc_32f(ntIn);
	uOut = ippsMalloc_32u(ntIn);

	printf("... Done.\n");
	return 0;
}

int XData::setFFT() {
	printf("Setting FFT");

	scale_for = 1.0f;
	scale_back = 1.0f / ntIn;

	dims_in[0] = nzIn;
	dims_in[1] = nyIn;
	dims_in[2] = nxIn;

	strides_in[0] = 0;
	strides_in[1] = nxIn * nyIn;
	strides_in[2] = nxIn;
	strides_in[3] = 1;
	
	fprintf(flog, "\nBox for FT: %i x %i x %i\n", nxIn, nyIn, nzIn);

	dims_out[0] = nzOut;
	dims_out[1] = nyOut;
	dims_out[2] = nxOut;

	strides_out[0] = 0;
	strides_out[1] = nxOut * nyOut;
	strides_out[2] = nxOut;
	strides_out[3] = 1;

	status_fft = DftiCreateDescriptor(&desc_real_for, DFTI_SINGLE, DFTI_REAL, 3, dims_in);
	status_fft = DftiSetValue(desc_real_for, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status_fft = DftiSetValue(desc_real_for, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status_fft = DftiSetValue(desc_real_for, DFTI_INPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_real_for, DFTI_OUTPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_real_for, DFTI_FORWARD_SCALE, scale_for);
	status_fft = DftiSetValue(desc_real_for, DFTI_ORDERING, DFTI_ORDERED);
	status_fft = DftiCommitDescriptor(desc_real_for);
	

	status_fft = DftiCreateDescriptor(&desc_comp_for, DFTI_SINGLE, DFTI_COMPLEX, 3, dims_in);
	status_fft = DftiSetValue(desc_comp_for, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
	status_fft = DftiSetValue(desc_comp_for, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status_fft = DftiSetValue(desc_comp_for, DFTI_INPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_comp_for, DFTI_OUTPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_comp_for, DFTI_FORWARD_SCALE, scale_for);
	status_fft = DftiSetValue(desc_comp_for, DFTI_ORDERING, DFTI_ORDERED);
	status_fft = DftiCommitDescriptor(desc_comp_for);


	status_fft = DftiCreateDescriptor(&desc_comp_back2, DFTI_SINGLE, DFTI_COMPLEX, 3, dims_in);
	status_fft = DftiSetValue(desc_comp_back2, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
	status_fft = DftiSetValue(desc_comp_back2, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status_fft = DftiSetValue(desc_comp_back2, DFTI_INPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_comp_back2, DFTI_OUTPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_comp_back2, DFTI_FORWARD_SCALE, scale_for);
	status_fft = DftiSetValue(desc_comp_back2, DFTI_BACKWARD_SCALE, scale_back);
	status_fft = DftiSetValue(desc_comp_back2, DFTI_ORDERING, DFTI_ORDERED);
	status_fft = DftiCommitDescriptor(desc_comp_back2);


	status_fft = DftiCreateDescriptor(&desc_comp_for_inplace, DFTI_SINGLE, DFTI_COMPLEX, 3, dims_in);
	status_fft = DftiSetValue(desc_comp_for_inplace, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
	status_fft = DftiSetValue(desc_comp_for_inplace, DFTI_PLACEMENT, DFTI_INPLACE);
	status_fft = DftiSetValue(desc_comp_for_inplace, DFTI_INPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_comp_for_inplace, DFTI_OUTPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_comp_for_inplace, DFTI_FORWARD_SCALE, scale_for);
	status_fft = DftiSetValue(desc_comp_for_inplace, DFTI_ORDERING, DFTI_ORDERED);
	status_fft = DftiCommitDescriptor(desc_comp_for_inplace);

	status_fft = DftiCreateDescriptor(&desc_comp_back, DFTI_SINGLE, DFTI_COMPLEX, 3, dims_out);
	status_fft = DftiSetValue(desc_comp_back, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
	status_fft = DftiSetValue(desc_comp_back, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status_fft = DftiSetValue(desc_comp_back, DFTI_INPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_comp_back, DFTI_OUTPUT_STRIDES, strides_out);
	status_fft = DftiSetValue(desc_comp_back, DFTI_BACKWARD_SCALE, scale_back);
	status_fft = DftiSetValue(desc_comp_back, DFTI_ORDERING, DFTI_ORDERED);
	status_fft = DftiCommitDescriptor(desc_comp_back);

	printf("... Done.\n");
	return 0;
}

int XData::rotateTarget() {
	printf("Rotating the target map");
	Ipp32f *map1, *map2, *mt1;
	Ipp32fc *mc1, *mc2;

	pmcShift = (Ipp32fc**)malloc(sizeof(Ipp32fc*) * ip->numberOfThreads);
	pmcT = (Ipp32fc **)malloc(sizeof(Ipp32fc *)*ip->pairsInGroup*ip->numberOfThreads);
	pmcT2 = (Ipp32fc **)malloc(sizeof(Ipp32fc *)*ip->pairsInGroup*ip->numberOfThreads);
	pmcCCbest = (Ipp32fc**)malloc(sizeof(Ipp32fc*) * ip->pairsInGroup * ip->numberOfThreads);
	pmcQbest = (Ipp32fc**)malloc(sizeof(Ipp32fc*) * ip->pairsInGroup * ip->numberOfThreads);
	vLocalQuat = (quat*)malloc(sizeof(quat) * 2 * ip->pairsInGroup);

	if (pmcT == nullptr || pmcT2 == nullptr || pmcShift == nullptr || pmcCCbest == nullptr || pmcQbest == nullptr) {
		fprintf(flog, "Error: memory for rotated target\n");
		printf("Error: memory for rotated target\n");
		return -1;
	}

	for (int i = 0; i < ip->pairsInGroup * ip->numberOfThreads; i++) {
		pmcT[i] = ippsMalloc_32fc(ntIn);
		pmcT2[i] = ippsMalloc_32fc(ntIn);
		if (pmcT[i] == nullptr || pmcT2[i] == nullptr) {
			fprintf(flog, "Error: memory for rotated target (%i)\n", i);
			printf("Error: memory for rotated target (%i)\n", i);
			return -2;
		}
	}

	for (int i = 0; i < ip->pairsInGroup * ip->numberOfThreads; i++) {
		pmcCCbest[i] = ippsMalloc_32fc(ntOut);
		pmcQbest[i] = ippsMalloc_32fc(ntOut);
		if (pmcCCbest[i] == nullptr || pmcQbest[i] == nullptr) {
			fprintf(flog, "Error: memory for rotated target CC/Q (%i)\n", i);
			printf("Error: memory for rotated target CC/Q (%i)\n", i);
			return -2;
		}
		ippsSet_32f(-1.0f, (Ipp32f*)pmcCCbest[i], 2 * ntOut);
		ippsZero_32s((Ipp32s*)pmcQbest[i], 2 * ntOut);
	}

	for (int i = 0; i < ip->numberOfThreads; i++) {
		pmcShift[i] = ippsMalloc_32fc(ntIn);
		if (pmcShift[i] == nullptr) {
			fprintf(flog, "Error: memory for rotated target/shift (%i)\n", i);
			return -3;
		}
	}

	map1 = ippsMalloc_32f(ntIn);
	map2 = ippsMalloc_32f(ntIn);
	mt1 = ippsMalloc_32f(ntIn);

	if (map1 == nullptr || map2 == nullptr || mt1 == nullptr) {
		fprintf(flog, "Error: memory for rotated target\n");
		return -2;
	}

	int nu;

	mapQMatrix(mapT, map1, mt1, vLocalQuat + 0, nxIn, nyIn, nzIn, 0);
	mapQMatrix(mapT, map2, mt1, vLocalQuat + 1, nxIn, nyIn, nzIn, 1);
	ippsRealToCplx_32f(map1, map2, pmcT[0], ntIn);
	mapQMatrix(mapT, map1, mt1, vLocalQuat + 2, nxIn, nyIn, nzIn, 2);
	mapQMatrix(mapT, map2, mt1, vLocalQuat + 3, nxIn, nyIn, nzIn, 3);
	ippsRealToCplx_32f(map1, map2, pmcT[1], ntIn);

	nu = 2;

	if (ip->nShapeROI == SHAPE_XY || ip->nShapeROI == SHAPE_XYZ) {
		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2*nu, nxIn, nyIn, nzIn, 12);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2*nu + 1, nxIn, nyIn, nzIn, 13);
		ippsRealToCplx_32f(map1, map2, pmcT[nu], ntIn);

		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2*nu + 2, nxIn, nyIn, nzIn, 14);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2*nu + 3, nxIn, nyIn, nzIn, 15);
		ippsRealToCplx_32f(map1, map2, pmcT[nu+1], ntIn);
		nu += 2;
	}

	if (ip->nShapeROI == SHAPE_XZ || ip->nShapeROI == SHAPE_XYZ) {
		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2 * nu, nxIn, nyIn, nzIn, 4);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2 * nu + 1, nxIn, nyIn, nzIn, 5);
		ippsRealToCplx_32f(map1, map2, pmcT[nu], ntIn);

		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2 * nu + 2, nxIn, nyIn, nzIn, 6);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2 * nu + 3, nxIn, nyIn, nzIn, 7);
		ippsRealToCplx_32f(map1, map2, pmcT[nu + 1], ntIn);
		nu += 2;
	}

	if (ip->nShapeROI == SHAPE_YZ || ip->nShapeROI == SHAPE_XYZ) {
		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2 * nu, nxIn, nyIn, nzIn, 20);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2 * nu + 1, nxIn, nyIn, nzIn, 21);
		ippsRealToCplx_32f(map1, map2, pmcT[nu], ntIn);

		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2 * nu + 2, nxIn, nyIn, nzIn, 22);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2 * nu + 3, nxIn, nyIn, nzIn, 23);
		ippsRealToCplx_32f(map1, map2, pmcT[nu + 1], ntIn);
		nu += 2;
	}

	if (ip->nShapeROI == SHAPE_XYZ) {
		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2 * nu, nxIn, nyIn, nzIn, 8);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2 * nu + 1, nxIn, nyIn, nzIn, 9);
		ippsRealToCplx_32f(map1, map2, pmcT[nu], ntIn);

		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2 * nu + 2, nxIn, nyIn, nzIn, 10);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2 * nu + 3, nxIn, nyIn, nzIn, 11);
		ippsRealToCplx_32f(map1, map2, pmcT[nu + 1], ntIn);
		nu += 2;

		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2 * nu, nxIn, nyIn, nzIn, 16);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2 * nu + 1, nxIn, nyIn, nzIn, 17);
		ippsRealToCplx_32f(map1, map2, pmcT[nu], ntIn);

		mapQMatrix(mapT, map1, mt1, vLocalQuat + 2 * nu + 2, nxIn, nyIn, nzIn, 18);
		mapQMatrix(mapT, map2, mt1, vLocalQuat + 2 * nu + 3, nxIn, nyIn, nzIn, 19);
		ippsRealToCplx_32f(map1, map2, pmcT[nu + 1], ntIn);
		nu += 2;
	}

	for (int i = 0; i < ip->pairsInGroup; i++) {
		ippsSqr_32f((Ipp32f*)pmcT[i], (Ipp32f*)pmcT2[i], 2 * ntIn);
		status_fft = DftiComputeForward(desc_comp_for_inplace, pmcT[i]);
		status_fft = DftiComputeForward(desc_comp_for_inplace, pmcT2[i]);
		//ippsConj_32fc_I(pmcT[i], ntIn);
		//ippsConj_32fc_I(pmcT2[i], ntIn);
	}

	fprintf(flog, "Memory for rotated target has been allocated.\n");

	ippsFree(map1); map1 = nullptr;
	ippsFree(map2); map2 = nullptr;
	ippsFree(mt1); mt1 = nullptr;
 	
	mc1 = ippsMalloc_32fc(ntIn);
	mc2 = ippsMalloc_32fc(ntIn);

	ippsZero_32f((Ipp32f *)mc1, 2 * ntIn);
	ippsZero_32f((Ipp32f *)mc2, 2 * ntIn);

	mc1[0].re = 1.0f;
	status_fft = DftiComputeForward(desc_comp_for_inplace, mc1);

	mc2[(irotCentGlobal[2]*nyIn + irotCentGlobal[1])*nxIn + irotCentGlobal[0]].re = 1.0f;
	status_fft = DftiComputeForward(desc_comp_for_inplace, mc2);

	ippsDiv_32fc(mc2, mc1, pmcShift[0], ntIn);
		
	
	ippsFree(mc1); mc1 = nullptr;
	ippsFree(mc2); mc2 = nullptr;

	for (int k = 1; k < ip->numberOfThreads; k++) {
		ippsCopy_32fc(pmcShift[0], pmcShift[k], ntIn);
		for (int i = 0; i < ip->pairsInGroup; i++) {
			ippsCopy_32fc(pmcT[i], pmcT[k*ip->pairsInGroup + i], ntIn);
			ippsCopy_32fc(pmcT2[i], pmcT2[k*ip->pairsInGroup + i], ntIn);
		}
	}

	printf("... Done.\n");
	return 0;
}

int XData::prepareMultiThread() {
	printf("Preparing multithreading");

	td = (TData **)malloc(sizeof(TData*)*ip->numberOfThreads);
	if (td == nullptr) {
		fprintf(flog, "Error: cannot create objects for multithreading.\n");
		printf("Error: cannot create objects for multithreading.\n");
		return -1;
	}

	for (int i = 0; i < ip->numberOfThreads; i++) {
		td[i] = new TData();
		td[i]->nxIn = nxIn;
		td[i]->nyIn = nyIn;
		td[i]->nzIn = nzIn;

		td[i]->isExternal = !(ip->sphericalMask);

		

		td[i]->nScale = ip->outputReduction;

		td[i]->nxOut = nxOut;
		td[i]->nyOut = nyOut;
		td[i]->nzOut = nzOut;
		
		td[i]->mx = mx;
		td[i]->mrx = mbx;

		td[i]->pairsInGroup = ip->pairsInGroup;
		td[i]->iRC_mask[0] = irotCentMask[0];
		td[i]->iRC_mask[1] = irotCentMask[1];
		td[i]->iRC_mask[2] = irotCentMask[2];

		td[i]->mapMS = pmMS[i];

		td[i]->mapXini = pmXini[i];
		td[i]->mapYini = pmYini[i];
		td[i]->mapZini = pmZini[i];

		td[i]->mapX = (Ipp32f *)(pmcTemp[nTemp*i + 0]);
		td[i]->mapY = (Ipp32f*)(pmcTemp[nTemp*i + 0]) + ntIn;
		td[i]->mapZ = (Ipp32f*)(pmcTemp[nTemp*i + 1]);

		td[i]->map1 = (Ipp32f*)(pmcTemp[nTemp * i + 1]) + ntIn;

		td[i]->mcRes = pmcTemp[nTemp * i];
		td[i]->mc1 = pmcTemp[nTemp * i + 1];
		
		td[i]->pmcT = pmcT + i * ip->pairsInGroup;
		td[i]->pmcT2 = pmcT2 + i * ip->pairsInGroup;
		td[i]->pmcCCbest = pmcCCbest + i * ip->pairsInGroup;
		td[i]->pmcQbest = pmcQbest + i * ip->pairsInGroup;
		td[i]->pmcRStd = pmcRStd + i * ip->pairsInGroup;

		td[i]->mM = (Ipp32f*)(pmcTemp[nTemp * i + 2]);
		td[i]->mS = (Ipp32f*)(pmcTemp[nTemp * i + 2]) + ntIn;

		
		td[i]->mcCC = pmcCC[i];
		td[i]->mcShift = pmcShift[i];

		td[i]->mcM = pmcM[i];
		td[i]->mcS = pmcS[i];

		DftiCopyDescriptor(desc_real_for, &(td[i]->desc_real_for));
		DftiCopyDescriptor(desc_comp_back, &(td[i]->desc_comp_back));
		DftiCopyDescriptor(desc_comp_back2, &(td[i]->desc_comp_back2));

		td[i]->setVariables();
	}

	printf("... Done.\n");
	return 0;
}

int XData::allocateThreadMemory() {
	printf("Allocating memory for multithreading");
	float pp, sx, sy, sz, val;

	// XYZini

	pmXini = (Ipp32f **)malloc(sizeof(Ipp32f*) *ip->numberOfThreads);
	pmYini = (Ipp32f **)malloc(sizeof(Ipp32f*) *ip->numberOfThreads);
	pmZini = (Ipp32f **)malloc(sizeof(Ipp32f*) *ip->numberOfThreads);


	if (pmXini == nullptr || pmYini == nullptr || pmZini == nullptr) {
		fprintf(flog, "Error: cannot allocate memory for Xini.\n");
		printf("Error: cannot allocate memory for Xini.\n");
		return -1;
	}

	for (int i = 0; i < ip->numberOfThreads; i++) {
		pmXini[i] = ippsMalloc_32f(mbt);
		pmYini[i] = ippsMalloc_32f(mbt);
		pmZini[i] = ippsMalloc_32f(mbt);

		if (pmXini[i] == nullptr || pmYini[i] == nullptr || pmZini[i] == nullptr) {
			fprintf(flog, "Error: cannot allocate memory for Xini.\n");
			printf("Error: cannot allocate memory for Xini.\n");
			return -2;
		}
	}

	pp = pixT / pixM;
	
	sx = -irotCentGlobal[0] * pp;
	sy = -irotCentGlobal[1] * pp;
	sz = -irotCentGlobal[2] * pp;
	
	// X
	ippsVectorSlope_32f(pmXini[0], mbx, sx, pp);

	for (int i = 1; i < mby; i++) {
		ippsCopy_32f(pmXini[0], pmXini[0] + i * mbx, mbx);
	}

	for (int j = 1; j < mbz; j++) {
		ippsCopy_32f(pmXini[0], pmXini[0] + j * mbx*mby, mbx*mby);
	}

	// Y
	for (int i = 0; i < mby; i++) {
		val = sy + float(i)*pp;
		ippsSet_32f(val, pmYini[0] + i * mbx, mbx);
	}

	for (int j = 1; j < mbz; j++) {
		ippsCopy_32f(pmYini[0], pmYini[0] + j * mbx*mby, mbx*mby);
	}

	//Z

	for (int i = 0; i < mbz; i++) {
		val = sz + float(i)*pp;
		ippsSet_32f(val, pmZini[0] + i * mbx*mby, mbx*mby);
	}

	for (int i = 1; i < ip->numberOfThreads; i++) {
		ippsCopy_32f(pmXini[0], pmXini[i], mbt);
		ippsCopy_32f(pmYini[0], pmYini[i], mbt);
		ippsCopy_32f(pmZini[0], pmZini[i], mbt);
	}


	pmcTemp = (Ipp32fc**)malloc(sizeof(Ipp32fc*) * ip->numberOfThreads *nTemp);
	for (int i = 0; i < ip->numberOfThreads*nTemp; i++) {
		pmcTemp[i] = ippsMalloc_32fc(ntIn);
		if (pmcTemp[i] == nullptr) {
			fprintf(flog, "Error: cannot allocate memory for cTemp.\n");
			printf("Error: cannot allocate memory for cTemp.\n");
			return -4;
		}
	}

	pmcS = (Ipp32fc **)malloc(sizeof(Ipp32fc*) * ip->numberOfThreads);
	for (int i = 0; i < ip->numberOfThreads; i++) {
		pmcS[i] = ippsMalloc_32fc(ntIn);
		if (pmcS[i] == nullptr) {
			fprintf(flog, "Error: cannot allocate memory for cS.\n");
			printf("Error: cannot allocate memory for cS.\n");
			return -4;
		}
	}

	pmcM = (Ipp32fc **)malloc(sizeof(Ipp32fc*) * ip->numberOfThreads);
	for (int i = 0; i < ip->numberOfThreads; i++) {
		pmcM[i] = ippsMalloc_32fc(ntIn);
		if (pmcM[i] == nullptr) {
			fprintf(flog, "Error: cannot allocate memory for cM.\n");
			printf("Error: cannot allocate memory for cM.\n");
			return -5;
		}
		
	}
	

	pmcCC = (Ipp32fc **)malloc(sizeof(Ipp32fc*) * ip->numberOfThreads);
	for (int i = 0; i < ip->numberOfThreads; i++) {
		pmcCC[i] = ippsMalloc_32fc(ntOut);
		if (pmcCC[i] == nullptr) {
			fprintf(flog, "Error: cannot allocate memory for cCC.\n");
			printf("Error: cannot allocate memory for cCC.\n");
			return -6;
		}
	}

	pmcRStd = (Ipp32fc**)malloc(sizeof(Ipp32fc*) * ip->numberOfThreads*ip->pairsInGroup);
	for (int i = 0; i < ip->numberOfThreads*ip->pairsInGroup; i++) {
		pmcRStd[i] = ippsMalloc_32fc(ntOut);
		if (pmcRStd[i] == nullptr) {
			fprintf(flog, "Error: cannot allocate memory for cRStd.\n");
			printf("Error: cannot allocate memory for cRStd.\n");
			return -6;
		}
	}
	
	printf("... Done.\n");
	return 0;
}

int XData::coreProcessing() {
	printf("Core processing\n");

	omp_set_num_threads(ip->numberOfThreads);
#pragma omp parallel
	{
		int tid, i;
		bool firstTime;
		TData *t;
		
		tid = omp_get_thread_num();
		t = td[tid];
		t->ncount = 0;
		firstTime = true;

#pragma omp for
		for (i = 0; i < nGroupMain; i++) {
			Ipp32u u1;
			ippsConvert_64f32f(matRot + 9 * i, t->matRot, 9);
			t->findXYZmaps();
			t->rotateSM();
			t->normaliseMapsM();
			if (t->isExternal || firstTime) {
				t->normaliseMapsFFT();
			}
			t->normaliseMapsS();
			u1 = ((unsigned int)i) << 5;
			for (int j = 0; j < ip->pairsInGroup; j++) {
				t->uIndex = u1 + 2 * j;
				t->mcTc = t->pmcT[j];
				t->mcT2c = t->pmcT2[j];
				t->mcRStd = t->pmcRStd[j];
				t->vCCbest = t->pmcCCbest[j];
				t->vQbest = t->pmcQbest[j];
				if (t->isExternal || firstTime) {
					t->findRStd();
				}
				t->findCC();
				t->findBest();
			}
			firstTime = false;
			t->ncount++;
			if (t->ncount % 10 == 0) {
				printf("Thread (%i): %f %% done\n", tid, float(t->ncount) / float(nGroupMain)*100.0);
			}
		}
	}

	printf("... Done.\n");
	return 0;
}


int XData::coreStdProcessing() {
	Ipp32f st, st2;
	int nSteps;
	if (ip->sphericalMask) {
		nSteps = 1;
	}
	else {
		nSteps = nGroupEstimate;
	}
	printf("Estimating Std for MT\n");

	omp_set_num_threads(ip->numberOfThreads);
#pragma omp parallel
	{
		int tid, i;
		TData *t;
		tid = omp_get_thread_num();
		t = td[tid];
		t->maxMTstd = 0.0f;
		t->ncount = 0;

#pragma omp for
		for (i = 0; i < nSteps; i++) {
			ippsConvert_64f32f(matRot + 9 * i, t->matRot, 9);
			t->findXYZmaps();
			t->rotateSM();
			t->normaliseMapsM();
			t->normaliseMapsFFT();
			t->normaliseMapsS();
			for (int j = 0; j < ip->pairsInGroup; j++) {
				t->mcTc = t->pmcT[j];
				t->mcT2c = t->pmcT2[j];
				t->mcRStd = t->pmcRStd[j];
				t->runFFTStd();
				if(t->maxV > t->maxMTstd)  t->maxMTstd = t->maxV;
			}
			t->ncount++;
			if (t->ncount % 10 == 0) {
				printf("MT Std Thread (%i): %f %% done\n", tid, float(t->ncount) / float(nGroupEstimate)*100.0);
			}
		}
	}

	maxStdValue = 0.0f;
	for (int i = 0; i < ip->numberOfThreads; i++) {
		if(maxStdValue < td[i]->maxMTstd){
			maxStdValue = td[i]->maxMTstd;
		}
	}

	maxStdValue = sqrt(maxStdValue);
	
	//maxStdValue = 0.290282f;

	st = maxStdValue * ip->signalLevel;
	st2 = st * st;

	fprintf(flog, "Max Std value for MT: %f\n", maxStdValue);
	fprintf(flog, "Std cut level for MT: %f\n", st);
	for (int i = 0; i < ip->numberOfThreads; i++) {
		td[i]->stdCut2 = st2;
	}

	printf("... Done.\n");
	return 0;
}




int XData::cropOutputData() {
	printf("Croping the output data");
	int ns;
	
	ns = nCornerShift / ip->outputReduction;
	nut = nux * nuy * nuz;

	for (int k = 0; k < nuz; k++) {
		for (int i = 0; i < nuy; i++) {
			ippsCopy_32f(((Ipp32f *)pmcCCbest[0]) + ((k + ns)*nyOut + i + ns)*nxOut + ns, vOut + (k*nuy + i)*nux, nux);
			ippsCopy_32s(((Ipp32s *)pmcQbest[0]) + ((k + ns) * nyOut + i + ns) * nxOut + ns, (Ipp32s *)uOut + (k * nuy + i) * nux, nux);
		}
	}
	
	printf("... Done.\n");
	return 0;
}

int XData::writeOutput() {
	printf("Writing the output data");
	targetOrig[0] = minfo->origin[0] + float(iCornerOrig[0])*pixT;
	targetOrig[1] = minfo->origin[1] + float(iCornerOrig[1])*pixT;
	targetOrig[2] = minfo->origin[2] + float(iCornerOrig[2])*pixT;

	pixTOut = pixT * ip->outputReduction;
	
	sprintf(outputFileName, "%s/%s_CC.mrc", ip->outputFolder, ip->outputPrefix);
	if(writeMRC(flog, outputFileName, vOut, nux, nuy, nuz, targetOrig, pixTOut) != 0) return -1;
	
	printf("... Done.\n");
	return 0;

}



int XData::startProcessing() {
	int idur;
	if (checkParameters() != 0) return -1;
	
	if (getNumberOfGroups() != 0) return -100;
	if (findRotationMatrices() != 0) return -101;
	if (readAtoms() != 0) return -2;
	if (findCentreOfMass() != 0) return -3;
	if (findGeomCentre() != 0) return -4;

	//return 0;
	//if (writeOrigXYZ() != 0) return -5;


	auto told = std::chrono::high_resolution_clock::now();
	if (findSearchMap() != 0) return -8;
	auto tnow = std::chrono::high_resolution_clock::now();
	idur = std::chrono::duration_cast<std::chrono::milliseconds>(tnow - told).count();
	fprintf(flog, "\nTime spent (Search map, in s): %i.%03i\n", idur / 1000, idur % 1000);

	told = std::chrono::high_resolution_clock::now();
	if (findMaskMap() != 0) return -9;
	tnow = std::chrono::high_resolution_clock::now();
	idur = std::chrono::duration_cast<std::chrono::milliseconds>(tnow - told).count();
	fprintf(flog, "\nTime spent (Mask map, in s): %i.%03i\n", idur / 1000, idur % 1000);

	if (correctMaskCentre() != 0) return -30;
	if (ip->sphericalMask) {
		if (findSpericalMask() != 0) return -31;
	}

	if (maskChoice() != 0) return -32;
	
	if (readTargetMap() != 0) return -10;
	if (cropTarget() != 0) return -11;
	if (setFFT() != 0) return -12;
	if (rotateTarget() != 0) return -13;
	if (allocateThreadMemory() != 0) return -14;
	if (prepareMultiThread() != 0) return -15;

	told = std::chrono::high_resolution_clock::now();
	if (coreStdProcessing() != 0) return -16;
	tnow = std::chrono::high_resolution_clock::now();
	idur = std::chrono::duration_cast<std::chrono::milliseconds>(tnow - told).count();
	fprintf(flog, "\nTime spent (Std, in s): %i.%03i\n", idur / 1000, idur % 1000);

	told = std::chrono::high_resolution_clock::now();
	if (coreProcessing() != 0) return -19;
	
	tnow = std::chrono::high_resolution_clock::now();
	idur = std::chrono::duration_cast<std::chrono::milliseconds>(tnow - told).count();
	fprintf(flog, "\nTime spent (Core, in s): %i.%03i\n", idur / 1000, idur % 1000);

	if (reduceOutputData() != 0) return -20;
	if (cropOutputData() != 0) return -21;
	
	if (writeOutput() != 0) return -22;
	if (findPeaks() != 0) return -23;
	
	return 0;
}
