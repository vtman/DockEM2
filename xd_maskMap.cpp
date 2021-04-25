#include "xdock_header.h"

int XData::findMaskMap() {
	printf("Forming the mask map");
	Ipp32f* vBall, * vq1, * vq2;
	Ipp32fc* mc1, * mc2, * mc3;
	Ipp32f* mapM, * v1_loc, * v2_loc;
	Ipp32f rr2, rr2max, s, dz, dy, dx, dx2, dz2, dy2, dzy2, fcount;
	Ipp32f val1, val2, v_old, v_new;
	int nfx, nfy, nfz, nft;

	int nsteps;
	

	

	mMask = ippsMalloc_32f(mt);

	mapM = mMask;

	nfx = mx / 2 + 1;
	nfy = mx;
	nfz = mx;
	nft = nfx * nfy * nfz;

	mc1 = ippsMalloc_32fc(nft);
	mc2 = ippsMalloc_32fc(nft);
	mc3 = ippsMalloc_32fc(nft);

	s = ip->maskRadius / ip->maskPixel;
	rr2max = s * s;


	vq1 = ippsMalloc_32f(mt);
	vq2 = ippsMalloc_32f(mt);
	vBall = ippsMalloc_32f(mt);

	ippsZero_32f(vBall, mt);

	for (int k = 0; k < mx; k++) {
		dz = float(myImin(k, mx - k));
		dz2 = dz * dz;
		for (int i = 0; i < mx; i++) {
			dy = float(myImin(i, mx - i));
			dy2 = dy * dy;
			dzy2 = dy2 + dz2;
			for (int j = 0; j < mx; j++) {
				dx = float(myImin(j, mx - j));
				dx2 = dx * dx;
				rr2 = dzy2 + dx2;
				if (rr2 <= rr2max) vBall[(k * mx + i) * mx + j] = 1.0f;
			}
		}
	}

	ippsSum_32f(vBall, mt, &fcount, IppHintAlgorithm::ippAlgHintAccurate);

	fcount -= 0.5f;


	ippsThreshold_LTValGTVal_32f(mSearch, vq1, mt, ip->maskLevel, 1.0f, ip->maskLevel, 0.0f);

	scale_for = 1.0f;
	scale_back = 1.0f / mt;

	dims_in[0] = mx;
	dims_in[1] = mx;
	dims_in[2] = mx;

	strides_in[0] = 0;
	strides_in[1] = mx * mx;
	strides_in[2] = mx;
	strides_in[3] = 1;

	dims_out[0] = nfz;
	dims_out[1] = nfy;
	dims_out[2] = nfx;

	strides_out[0] = 0;
	strides_out[1] = nfx * nfy;
	strides_out[2] = nfx;
	strides_out[3] = 1;

	status_fft = DftiCreateDescriptor(&desc_real_for, DFTI_SINGLE, DFTI_REAL, 3, dims_in);
	status_fft = DftiCreateDescriptor(&desc_real_back, DFTI_SINGLE, DFTI_REAL, 3, dims_in);

	status_fft = DftiSetValue(desc_real_for, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status_fft = DftiSetValue(desc_real_for, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status_fft = DftiSetValue(desc_real_for, DFTI_INPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_real_for, DFTI_OUTPUT_STRIDES, strides_out);
	status_fft = DftiSetValue(desc_real_for, DFTI_FORWARD_SCALE, scale_for);

	status_fft = DftiSetValue(desc_real_back, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	status_fft = DftiSetValue(desc_real_back, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status_fft = DftiSetValue(desc_real_back, DFTI_INPUT_STRIDES, strides_out);
	status_fft = DftiSetValue(desc_real_back, DFTI_OUTPUT_STRIDES, strides_in);
	status_fft = DftiSetValue(desc_real_back, DFTI_BACKWARD_SCALE, scale_back);

	status_fft = DftiCommitDescriptor(desc_real_for);
	status_fft = DftiCommitDescriptor(desc_real_back);

	status_fft = DftiComputeForward(desc_real_for, vq1, mc1);
	status_fft = DftiComputeForward(desc_real_for, vBall, mc2);

	ippsConj_32fc_I(mc2, nft);

	ippsMul_32fc_I(mc2, mc1, nft);

	status_fft = DftiComputeBackward(desc_real_back, mc1, vq2);

	ippsThreshold_LTValGTVal_32f(vq2, vq1, mt, fcount, 0.0f, fcount, 1.0f);

	status_fft = DftiComputeForward(desc_real_for, vq1, mc1);
	ippsMul_32fc_I(mc2, mc1, nft);
	status_fft = DftiComputeBackward(desc_real_back, mc1, vq2);

	ippsThreshold_LTValGTVal_32f(vq2, vq1, mt, 0.5f, 0.0f, 0.5f, 1.0f);

	ippsZero_32f(vBall, mt);

	vBall[0] = 1.0f;
	vBall[1] = 1.0f;
	vBall[mx] = 1.0f;
	vBall[mx * mx] = 1.0f;
	vBall[mx - 1] = 1.0f;
	vBall[(mx - 1) * mx] = 1.0f;
	vBall[(mx - 1) * mx * mx] = 1.0f;

	status_fft = DftiComputeForward(desc_real_for, vBall, mc1);
	ippsConj_32fc_I(mc1, nft);

	ippsZero_32f(vq2, mt);

	val1 = float(mx);

	for (int k = 0; k < mx; k++) {
		for (int i = 0; i < mx; i++) {
			ippsSum_32f(vq1 + (k * mx + i) * mx, mx, &val2, IppHintAlgorithm::ippAlgHintAccurate);
			if (abs(val2 - val1) < 0.5f) {
				ippsSet_32f(1.0f, vq2 + (k * mx + i) * mx, mx);
				continue;
			}
			else {

				v1_loc = vq1 + (k * mx + i) * mx;
				v2_loc = vq2 + (k * mx + i) * mx;

				for (int j = 0; j < mx; j++) {
					if (v1_loc[j] < 0.5f)break;
					v2_loc[j] = 1.0f;
				}

				for (int j = mx - 1; j > 0; j--) {
					if (v1_loc[j] < 0.5f)break;
					v2_loc[j] = 1.0f;
				}
			}
		}
	}

	for (int k = 0; k < mx; k++) {
		for (int j = 0; j < mx; j++) {
			v1_loc = vq1 + (k * mx) * mx + j;
			v2_loc = vq2 + (k * mx) * mx + j;

			for (int i = 0; i < mx; i++) {
				if (v1_loc[i * mx] < 0.5f)break;
				v2_loc[i * mx] = 1.0f;
			}

			for (int i = mx - 1; i > 0; i--) {
				if (v1_loc[i * mx] < 0.5f)break;
				v2_loc[i * mx] = 1.0f;
			}
		}
	}


	for (int i = 0; i < mx; i++) {
		for (int j = 0; j < mx; j++) {
			v1_loc = vq1 + i * mx + j;
			v2_loc = vq2 + i * mx + j;

			for (int k = 0; k < mx; k++) {
				if (v1_loc[k * mx * mx] < 0.5f)break;
				v2_loc[k * mx * mx] = 1.0f;
			}

			for (int k = mx - 1; k > 0; k--) {
				if (v1_loc[k * mx * mx] < 0.5f)break;
				v2_loc[k * mx * mx] = 1.0f;
			}
		}
	}


	printf("\nStart\n");


	v_new = 0.0f;

	nsteps = 0;

	do {
		v_old = v_new;
		status_fft = DftiComputeForward(desc_real_for, vq2, mc2);

		ippsMul_32fc_I(mc1, mc2, nft);

		status_fft = DftiComputeBackward(desc_real_back, mc2, vq2);

		ippsThreshold_LTValGTVal_32f_I(vq2, mt, 0.5, 0.0f, 0.5, 1.0f);

		ippsMul_32f_I(vq1, vq2, mt);

		ippsSum_32f(vq2, mt, &v_new, IppHintAlgorithm::ippAlgHintAccurate);
		nsteps++;
		printf("%f\n", v_new);

	} while (abs(v_old - v_new) > 0.5f);

	printf("\nNumber of steps: %i\n", nsteps);

	ippsSubCRev_32f(vq2, 1.0f, mapM, mt);

	if (!(ip->sphericalMask)) {

		DftiFreeDescriptor(&desc_real_back); desc_real_back = nullptr;
		DftiFreeDescriptor(&desc_real_for); desc_real_for = nullptr;
	}

	ippsFree(mc1); mc1 = nullptr;
	ippsFree(mc2); mc2 = nullptr;
	ippsFree(mc3); mc3 = nullptr;

	ippsFree(vq1); vq1 = nullptr;
	ippsFree(vq2); vq2 = nullptr;
	ippsFree(vBall); vBall = nullptr;

	

	

	return 0;

}

int XData::correctMaskCentre() {
	Ipp32f *mapM;
	Ipp32s* iv1, *iX, * iY, * iZ;
	
	int ind, ntot;
	int px, py, pz;
	int qk, qi, qj;
	int i_min, i_max, j_min, j_max, k_min, k_max;
	int indX, indY, indZ;

	Ipp64f sx, sy, sz, sref, s2, snew;
	Ipp64f* dRef, r1;
	
	bool t, isReset;
	
	mapM = mMask;

	dRef = ippsMalloc_64f(2 * mx);

	for (int i = 0; i < 2 * mx; i++) {
		r1 = double(i - mx);
		dRef[i] = r1 * r1;
	}

	iv1 = ippsMalloc_32s(mt);

	ntot = 0;

	for (int k = 1; k < mx - 1; k++) {
		for (int i = 1; i < mx - 1; i++) {
			for (int j = 1; j < mx - 1; j++) {
				ind = (k * mx + i) * mx + j;
				if (mapM[ind] < 0.5) continue;
				t = true;
				for (int kk = k - 1; kk <= k + 1; kk++) {
					for (int ii = i - 1; ii <= i + 1; ii++) {
						for (int jj = j - 1; jj <= j + 1; jj++) {
							if (mapM[(kk * mx + ii) * mx + jj] > 0.5) continue;
							t = false;
							break;
						}
						if (!t)break;
					}
					if (!t)break;
				}
				if (t)continue;
				iv1[ntot] = ind;
				ntot++;
			}
		}
	}


	iX = ippsMalloc_32s(ntot);
	iY = ippsMalloc_32s(ntot);
	iZ = ippsMalloc_32s(ntot);

	for (int i = 0; i < ntot; i++) {
		iX[i] = iv1[i] % mx;
		iY[i] = (iv1[i] / mx) % mx;
		iZ[i] = iv1[i] / (mx * mx);
	}

	printf("Number of border voxels: %i\n", ntot);


	j_min = iv1[0] % mx;
	i_min = (iv1[0] / mx) % mx;
	k_min = iv1[0] / (mx * mx);

	j_max = j_min;
	i_max = i_min;
	k_max = k_min;

	for (int p = 1; p < ntot; p++) {
		qj = iv1[p] % mx;
		qi = (iv1[p] / mx) % mx;
		qk = iv1[p] / (mx * mx);
		j_min = myImin(j_min, qj);
		j_max = myImax(j_max, qj);
		i_min = myImin(i_min, qi);
		i_max = myImax(i_max, qi);
		k_min = myImin(k_min, qk);
		k_max = myImax(k_max, qk);
	}

	fprintf(flog, "Mask box: [%i, %i] x [%i, %i] x [%i, %i]\n", j_min, j_max, i_min, i_max, k_min, k_max);

	indX = irotCentMask[0];
	indY = irotCentMask[1];
	indZ = irotCentMask[2];

	sref = 0.0;

	for (int s = 0; s < ntot; s++) {
		sx = dRef[iX[s] - indX + mx];
		sy = dRef[iY[s] - indY + mx];
		sz = dRef[iZ[s] - indZ + mx];
		s2 = sx + sy + sz;
		sref = myDmax(sref, s2);
	}
	fprintf(flog, "Original centre of mask (%i, %i, %i), radius (pixels, mask): %f\n", indX, indY, indZ, sqrt(sref));
	
	do {
		isReset = false;

		for (int kk = -1; kk <= 1; kk++) {
			pz = indZ + kk;
			for (int ii = -1; ii <= 1; ii++) {
				py = indY + ii;
				for (int jj = -1; jj <= 1; jj++) {
					px = indX + jj;
					if (kk == 0 && ii == 0 && jj == 0)continue;
					snew = 0.0;

					for (int s = 0; s < ntot; s++) {
						sx = dRef[iX[s] - px + mx];
						sy = dRef[iY[s] - py + mx];
						sz = dRef[iZ[s] - pz + mx];
						s2 = sx + sy + sz;
						snew = myDmax(snew, s2);
						if (snew >= sref) break;
					}
					if (snew >= sref) continue;

					indX = px;
					indY = py;
					indZ = pz;
					sref = snew;
					isReset = true;
					break;
				}
				if (isReset)break;
			}
			if (isReset)break;
		}
	} while (isReset);

	iMaskExternal[0] = indX;
	iMaskExternal[1] = indY;
	iMaskExternal[2] = indZ;

	radMaskExternal = sqrt(sref);

	fprintf(flog, "Updated external centre of mask (%i, %i, %i), radius (pixels, mask): %f\n", iMaskExternal[0], iMaskExternal[1], iMaskExternal[2], radMaskExternal);

	ippsFree(iv1); iv1 = nullptr;

	ippsFree(iX); iX = nullptr;
	ippsFree(iY); iY = nullptr;
	ippsFree(iZ); iZ = nullptr;
	ippsFree(dRef); dRef = nullptr;

	return 0;
}


int XData::findSpericalMask() {
	Ipp64f r1, r2, r3, rr;
	Ipp64f dx, dy, dz, d1, d2, d3;
	Ipp32f* vIO;
	Ipp32fc* vc1, * vc2, *vc3;
	Ipp32f v1, v2;
	int ind;

	int nfx, nfy, nfz, nft;
	nfx = mx / 2 + 1;
	nfy = mx;
	nfz = mx;
	nft = nfx * nfy * nfz;

	vc1 = ippsMalloc_32fc(nft);
	vc2 = ippsMalloc_32fc(nft);
	vc3 = ippsMalloc_32fc(nft);
	vIO = ippsMalloc_32f(mx * mx * mx);

	r1 = 0.0;
	r2 = radMaskExternal;

	status_fft = DftiComputeForward(desc_real_for, mMask, vc1);


	do {
		r3 = 0.5 * (r1 + r2);
		rr = r3 * r3;
		ippsZero_32f(vIO, mx * mx * mx);
		for (int k = 0; k < mx; k++) {
			dz = myImin(k, mx - k);
			d1 = dz * dz;
			for (int i = 0; i < mx; i++) {
				dy = myImin(i, mx - i);
				d2 = d1 + dy * dy;
				for (int j = 0; j < mx; j++) {
					dx = myImin(j, mx - j);
					d3 = d2 + dx * dx;
					if (d3 < rr) vIO[(mx * k + i) * mx + j] = 1.0f;
				}
			}
		}
		ippsDotProd_32f(vIO, vIO, mx * mx * mx, &v1);
		//printf("%f\t%f\n", r3, v1);
		v1 -= 0.5f;

		status_fft = DftiComputeForward(desc_real_for, vIO, vc2);

		ippsMul_32fc_I(vc1, vc2, nft);

		status_fft = DftiComputeBackward(desc_real_back, vc2, vIO);

		ippsMaxIndx_32f(vIO, mt, &v2, &ind);
		//printf("Max value: %f\n", v2);

		if (v2 < v1) {
			r2 = r3;
		}else{
			r1 = r3;
			iMaskInternal[0] = ind % mx;
			iMaskInternal[1] = (ind / mx) % mx;
			iMaskInternal[2] = ind / (mx * mx);
			radMaskInternal = r1;
			printf("Int radius: %f (%i, %i, %i)\n", r1, iMaskInternal[0], iMaskInternal[1], iMaskInternal[2]);
			
		}

	} while (r2 - r1 > 0.1);

	ippsZero_32f(mMask, mt);

	rr = radMaskInternal * radMaskInternal;

	for (int k = 0; k < mx; k++) {
		dz = double(k - iMaskInternal[2]);
		d1 = dz * dz;
		for (int i = 0; i < mx; i++) {
			dy = double(i - iMaskInternal[1]);
			d2 = d1 + dy * dy;
			for (int j = 0; j < mx; j++) {
				dx = double(j - iMaskInternal[0]);
				d3 = d2 + dx * dx;
				if (d3 < rr) mMask[(mx * k + i) * mx + j] = 1.0f;
			}
		}
	}

	ippsFree(vIO); vIO = nullptr;
	ippsFree(vc1); vc1 = nullptr;
	ippsFree(vc2); vc2 = nullptr;
	ippsFree(vc3); vc3 = nullptr;

	DftiFreeDescriptor(&desc_real_back); desc_real_back = nullptr;
	DftiFreeDescriptor(&desc_real_for); desc_real_for = nullptr;

	fprintf(flog, "\nInternal sphere (%i, %i, %i), radius (pixels, mask): %f\n\n", iMaskInternal[0], iMaskInternal[1], iMaskInternal[2], radMaskInternal);

	return 0;
}

int XData::maskChoice() {
	if (ip->writeMask) {
		sprintf(outputFileName, "%s/%s_mapM.mrc", ip->outputFolder, ip->outputPrefix);
		if (writeMRC(flog, outputFileName, mMask, mx, mx, mx, maskOrig, pixM) != 0) return -1;
	}

	pmMS = (Ipp32f**)malloc(sizeof(Ipp32f*) * ip->numberOfThreads);

	for (int i = 0; i < ip->numberOfThreads; i++) {
		pmMS[i] = ippsMalloc_32f(mt);
		if (pmMS[i] == nullptr) {
			fprintf(flog, "Error: cannot allocate memory for mapM.\n");
			printf("Error: cannot allocate memory for mapM.\n");
			return -1;
		}
	}


	ippsMul_32f_I(mMask, mSearch, mt);

	ippsSubCRev_32f_I(1.0f, mMask, mt);

	ippsSub_32f(mMask, mSearch, pmMS[0], mt);

	for (int i = 1; i < ip->numberOfThreads; i++) {
		ippsCopy_32f(pmMS[0], pmMS[i], mt);
	}

	if (ip->sphericalMask) {
		irotCentMask[0] = iMaskInternal[0];
		irotCentMask[1] = iMaskInternal[1];
		irotCentMask[2] = iMaskInternal[2];
		borderSize = radMaskInternal * pixM;
	}
	else {
		irotCentMask[0] = iMaskExternal[0];
		irotCentMask[1] = iMaskExternal[1];
		irotCentMask[2] = iMaskExternal[2];
		borderSize = radMaskExternal * pixM;
	}


	fprintf(flog, "Rotation centre (mask, pixels): %i, %i, %i\n", irotCentMask[0], irotCentMask[1], irotCentMask[2]);
	fprintf(flog, "Border size (in A): %f\n", borderSize);

	printf("... Done.\n");
	return 0;
}

