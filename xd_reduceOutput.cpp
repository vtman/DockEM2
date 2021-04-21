#include "xdock_header.h"

int XData::reduceOutputData() {
	printf("Rotating the output data");
	Ipp32f* mt1;
	Ipp32f* v1, * v2, * v3, * v4;
	quat* vQuat;


	v1 = (Ipp32f*)pmcTemp[0];
	v2 = (Ipp32f*)pmcTemp[1];
	v3 = (Ipp32f*)pmcS[0];
	v4 = (Ipp32f*)pmcM[0];
	mt1 = (Ipp32f*)pmcRStd[0];


	for (int k = 0; k < ip->pairsInGroup; k++) {
		for (int i = 1; i < ip->numberOfThreads; i++) {
			findBestAll(pmcCCbest[i * ip->pairsInGroup + k], pmcQbest[i * ip->pairsInGroup + k], pmcCCbest[k], pmcQbest[k], ntOut);
		}
	}

	vQuat = (quat*)malloc(sizeof(quat) * ip->pairsInGroup * 2);



	int nu;

	ippsCplxToReal_32fc(pmcCCbest[0], v1, v2, ntOut);
	mapQMatrix(v1, v3, mt1, vQuat + 0, nxOut, nyOut, nzOut, 0);
	mapQMatrix(v2, v4, mt1, vQuat + 1, nxOut, nyOut, nzOut, 1);
	ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[0], ntOut);
	ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[0]) + ntOut, ntOut);

	ippsCplxToReal_32fc(pmcQbest[0], v1, v2, ntOut);
	mapQMatrix(v1, v3, mt1, vQuat + 0, nxOut, nyOut, nzOut, 0);
	mapQMatrix(v2, v4, mt1, vQuat + 1, nxOut, nyOut, nzOut, 1);
	ippsCopy_32f(v3, (Ipp32f*)pmcQbest[0], ntOut);
	ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[0]) + ntOut, ntOut);

	ippsCplxToReal_32fc(pmcCCbest[1], v1, v2, ntOut);
	mapQMatrix(v1, v3, mt1, vQuat + 2, nxOut, nyOut, nzOut, 2);
	mapQMatrix(v2, v4, mt1, vQuat + 3, nxOut, nyOut, nzOut, 3);
	ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[1], ntOut);
	ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[1]) + ntOut, ntOut);

	ippsCplxToReal_32fc(pmcQbest[1], v1, v2, ntOut);
	mapQMatrix(v1, v3, mt1, vQuat + 2, nxOut, nyOut, nzOut, 2);
	mapQMatrix(v2, v4, mt1, vQuat + 3, nxOut, nyOut, nzOut, 3);
	ippsCopy_32f(v3, (Ipp32f*)pmcQbest[1], ntOut);
	ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[1]) + ntOut, ntOut);


	nu = 2;

	if (ip->nShapeROI == SHAPE_XY || ip->nShapeROI == SHAPE_XYZ) {
		ippsCplxToReal_32fc(pmcCCbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 12);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 14);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 12);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 14);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcCCbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 13);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 15);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu + 1]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 13);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 15);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu + 1]) + ntOut, ntOut);

		nu += 2;
	}

	if (ip->nShapeROI == SHAPE_XZ || ip->nShapeROI == SHAPE_XYZ) {
		ippsCplxToReal_32fc(pmcCCbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 4);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 7);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 4);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 7);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcCCbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 6);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 5);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu + 1]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 6);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 5);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu + 1]) + ntOut, ntOut);

		nu += 2;
	}

	if (ip->nShapeROI == SHAPE_YZ || ip->nShapeROI == SHAPE_XYZ) {
		ippsCplxToReal_32fc(pmcCCbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 20);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 21);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 20);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 21);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcCCbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 23);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 22);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu + 1]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 23);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 22);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu + 1]) + ntOut, ntOut);

		nu += 2;
	}

	if (ip->nShapeROI == SHAPE_XYZ) {
		ippsCplxToReal_32fc(pmcCCbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 16);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 19);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 16);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 19);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcCCbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 17);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 18);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu + 1]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 17);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 18);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu + 1]) + ntOut, ntOut);

		nu += 2;

		ippsCplxToReal_32fc(pmcCCbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 8);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 10);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu, nxOut, nyOut, nzOut, 8);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 1, nxOut, nyOut, nzOut, 10);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcCCbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 11);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 9);
		ippsCopy_32f(v3, (Ipp32f*)pmcCCbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcCCbest[nu + 1]) + ntOut, ntOut);

		ippsCplxToReal_32fc(pmcQbest[nu + 1], v1, v2, ntOut);
		mapQMatrix(v1, v3, mt1, vQuat + 2 * nu + 2, nxOut, nyOut, nzOut, 11);
		mapQMatrix(v2, v4, mt1, vQuat + 2 * nu + 3, nxOut, nyOut, nzOut, 9);
		ippsCopy_32f(v3, (Ipp32f*)pmcQbest[nu + 1], ntOut);
		ippsCopy_32f(v4, ((Ipp32f*)pmcQbest[nu + 1]) + ntOut, ntOut);

		nu += 2;
	}

	free(vQuat); vQuat = nullptr;

	v1 = (Ipp32f*)pmcCCbest[0];
	v2 = (Ipp32f*)pmcQbest[0];

	for (int k = 1; k < ip->pairsInGroup; k++) {
		v3 = (Ipp32f*)pmcCCbest[k];
		v4 = (Ipp32f*)pmcQbest[k];
		for (int i = 0; i < 2 * ntOut; i++) {
			if (v3[i] > v1[i]) {
				v1[i] = v3[i];
				v2[i] = v4[i];
			}
		}
	}

	v3 = (Ipp32f*)pmcCCbest[0] + ntOut;
	v4 = (Ipp32f*)pmcQbest[0] + ntOut;

	for (int i = 0; i < ntOut; i++) {
		if (v3[i] > v1[i]) {
			v1[i] = v3[i];
			v2[i] = v4[i];
		}
	}

	printf("... Done.\n");
	return 0;
}
