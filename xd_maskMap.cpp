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

	DftiFreeDescriptor(&desc_real_back); desc_real_back = nullptr;
	DftiFreeDescriptor(&desc_real_for); desc_real_for = nullptr;

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
	Ipp32s* iv1, *ivX, *ivY, *ivZ;
	double rbest, uxyz2, uz2, uy2, ux2, uyz2;
	int idx, idy, idz, ind, ntot, nL;
	bool t;
	int i_min, i_max, j_min, j_max, k_min, k_max;

	nL = 100;

	mapM = mMask;

	ivX = ippsMalloc_32s(nL);
	ivY = ippsMalloc_32s(nL);
	ivZ = ippsMalloc_32s(nL);

	iv1 = ippsMalloc_32s(mt);

	/*ntot = 0;

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

	printf("Number of border voxels: %i\n", ntot);

	ippsSortAscend_32s_I(iv1, ntot);

	int qk, qi, qj;

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

	printf("mx: %i\n", mx);

	printf("Mask box: [%i, %i] x [%i, %i] x [%i, %i]\n", j_min, j_max, i_min, i_max, k_min, k_max);

	ivX[0] = iv1[0] % mx;
	ivY[0] = (iv1[0] / mx) % mx;
	ivZ[0] = iv1[0] / (mx * mx);

	Ipp64f dx, dy, dz, r2, rmax, rmin, w, rbig;

	rbig = 9.e10;

	for (int s = 1; s < nL; s++) {
		ind = -1;
		rmax = -1.0;
		for (int p = 0; p < ntot; p++) {
			qj = iv1[p] % mx;
			qi = (iv1[p] / mx) % mx;
			qk = iv1[p] / (mx * mx);

			w = rbig;

			for (int q = 0; q < s; q++) {
				dx = double(ivX[q] - qj);
				dy = double(ivY[q] - qi);
				dz = double(ivZ[q] - qk);
				r2 = dx * dx + dy * dy + dz * dz;
				if (s == 5 && p == 17 && q == 4) {
					int iuo = 0;
				}
				w = myDmin(w, r2);
			}
			if (w <= rmax) continue;
			rmax = w;
			ind = p;
		}

		ivX[s] = iv1[ind] % mx;
		ivY[s] = (iv1[ind] / mx) % mx;
		ivZ[s] = iv1[ind] / (mx * mx);

		printf("\t%i\t %i, %i, %i\t%lf\t%i\n", s, ivX[s], ivY[s], ivZ[s], sqrt(rmax), ind);
	}

	int indX, indY, indZ;

	indX = 0;
	indY = 0;
	indZ = 0;
	rmin = rbig;

	for (int k = 0; k < mx; k++) {
		for (int i = 0; i < mx; i++) {
			for (int j = 0; j < mx; j++) {
				w = 0.0;
				
				for (int s = 0; s < nL; s++) {
					dx = double(ivX[s] - j);
					dy = double(ivY[s] - i);
					dz = double(ivZ[s] - k);
					r2 = dx * dx + dy * dy + dz * dz;
					w = myDmax(w, r2);
				}

				if (w >= rmin)continue;
				rmin = w;
				indZ = k;
				indY = i;
				indX = j;
			}
		}
	}

	printf("rmin: %f\n", sqrt(rmin));
	
	irotCentMask[0] = indX;
	irotCentMask[1] = indY;
	irotCentMask[2] = indZ;

	*/

	rbest = -1.0;

	for (int k = 0; k < mx; k++) {
		idz = k - irotCentMask[2];
		uz2 = double(idz * idz);
		for (int i = 0; i < mx; i++) {
			idy = i - irotCentMask[1];
			uy2 = double(idy * idy);
			uyz2 = uz2 + uy2;
			for (int j = 0; j < mx; j++) {
				idx = j - irotCentMask[0];
				ux2 = double(idx * idx);
				uxyz2 = uyz2 + ux2;
				if (mapM[(k * mx + i) * mx + j] < 0.8) continue;
				rbest = myDmax(rbest, uxyz2);
			}
		}
	}


	//printf("New center: %i, %i, %i (%f)\n", indX, indY, indZ, sqrt(rbest));
	//printf("old center: %i, %i, %i (%f)\n", irotCentMask[0], irotCentMask[1], irotCentMask[2], sqrt(rbest));

	ippsFree(iv1); iv1 = nullptr;
	ippsFree(ivX); ivX = nullptr;
	ippsFree(ivY); ivY = nullptr;
	ippsFree(ivZ); ivZ = nullptr;


	if (ip->writeMask) {
		sprintf(outputFileName, "%s/%s_mapM.mrc", ip->outputFolder, ip->outputPrefix);
		if (writeMRC(flog, outputFileName, mapM, mx, mx, mx, maskOrig, pixM) != 0) return -1;
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

	borderSize = (float)sqrt(rbest) * pixM;

	fprintf(flog, "Border size (in A): %f\n", borderSize);

	printf("... Done.\n");

	return 0;
}

