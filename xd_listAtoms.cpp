#include "xdock_header.h"

ListAtoms::ListAtoms() {
	nGauss = 4;
	nAtoms = 0;
	index = nullptr;
	vx = nullptr;
	vy = nullptr;
	vz = nullptr;
	indexAtom = nullptr;
	indexD = nullptr;

	vGValues = nullptr;
	vGxyz = nullptr;
	vRmax = nullptr;
	vdR = nullptr;
	vVmax = nullptr;
	vT = nullptr;
	vX = nullptr;
	vY = nullptr;
	vZ = nullptr;
}

ListAtoms::~ListAtoms() {
	if (indexD != nullptr) { free(indexD); indexD = nullptr;}
	if (index != nullptr) { free(index); index = nullptr; }
	if (vx != nullptr) { free(vx); vx = nullptr; }
	if (vy != nullptr) { free(vy); vy = nullptr; }
	if (vz != nullptr) { free(vz); vz = nullptr; }
	if (indexAtom != nullptr) { free(indexAtom); indexAtom = nullptr; }
	if (vRmax != nullptr) {
		ippsFree(vRmax); vRmax = nullptr;
	}
	if (vdR != nullptr) {
		ippsFree(vdR); vdR = nullptr;
	}
	if (vVmax != nullptr) {
		ippsFree(vVmax); vVmax = nullptr;
	}
	if (vGValues != nullptr) {
		for (int i = 0; i < numDifAtoms; i++) {
			ippsFree(vGValues[i]);
		}
		free(vGValues); vGValues = nullptr;
	}

	if (vGxyz != nullptr) {
		for (int i = 0; i < numDifAtoms; i++) {
			ippsFree(vGxyz[i]);
		}
		free(vGxyz); vGxyz = nullptr;
	}

	if (vT != nullptr) {
		ippsFree(vT); vT = nullptr;
	}
	if (vX != nullptr) {
		ippsFree(vX); vX = nullptr;
	}
	if (vY != nullptr) {
		ippsFree(vY); vY = nullptr;
	}
	if (vZ != nullptr) {
		ippsFree(vZ); vZ = nullptr;
	}
}

double ListAtoms::distBall(double xx, double yy, double zz, int *ind) {
	double dx, dy, dz, dr2, drbest;
	drbest = -1.0;
	*ind = -1;

	for (int i = 0; i < nAtoms; i++) {
		dx = double(vx[i]) - xx;
		dy = double(vy[i]) - yy;
		dz = double(vz[i]) - zz;

		dr2 = dx * dx + dy * dy + dz * dz;
		if (dr2 > drbest) {
			drbest = dr2;
			*ind = i;
		}
	}

	return sqrt(drbest);
}


int ListAtoms::findMaxStep(double xx, double yy, double zz, double *veca, double R, double* astep) {
	double dx, dy, dz, D;
	double ca, cb, cc, aa, amin;

	amin = 99999999.0;

	for (int i = 0; i < nAtoms; i++) {
		dx = double(vx[i]) - xx;
		dy = double(vy[i]) - yy;
		dz = double(vz[i]) - zz;

		ca = 0.0;
		cb = 0.0;
		cc = 0.0;

		ca = veca[0] * veca[0] + veca[1] * veca[1] + veca[2] * veca[2];
		cb = dx * veca[0] + dy * veca[1] + dz * veca[2];
		cc = dx * dx + dy * dy + dz * dz - R * R;

		D = cb * cb - ca * cc;
		if (D < 0.0) continue;

		aa = (cb + sqrt(D)) / ca;
		amin = myDmin(amin, aa);

	}

	*astep = amin;
	return 0;
}

Ipp32f ListAtoms::getValue(Ipp32f xx, int ind, int ig) {
	Ipp32f *vg, w, w0, r, x;
	int ix;
	if (xx < 0.0f) {
		x = -xx;
	}
	else {
		x = xx;
	}
	vg = vGxyz[ind] + ig * nCurves;
	w0 = x / vdR[ind];
	ix = (int)(floorf(w0));
	w = w0 - float(ix);
	if (ix <= 0) return vg[0];
	if (ix >= nCurves - 1) return 0.0f;
	r = (1.0f - w)*vg[ix] + w * vg[ix + 1];
	return r;
}

int ListAtoms::findMap(Ipp32f *vm, int nx, int ny, int nz, Ipp32f pixelSize, int icx, int icy, int icz, Ipp32f cx, Ipp32f cy, Ipp32f cz) {
	Ipp32f x1, y1, z1, x, y, z;
	int ix_left, ix_right, iy_left, iy_right, iz_left, iz_right;
	int idx, idy, idz, indA, nt;
	float ss, vLeft, vRight;

	nt = nx * ny * nz;

	x1 = cx - pixelSize * float(icx);
	y1 = cy - pixelSize * float(icy);
	z1 = cz - pixelSize * float(icz);
		
	vX = ippsMalloc_32f(nx);
	vT = ippsMalloc_32f(nx);
	vY = ippsMalloc_32f(ny);
	vZ = ippsMalloc_32f(nz);
	ippsZero_32f(vm, nt);

	for (int k = 0; k < nAtoms; k++) {
		indA = indexD[k];
		x = vx[k] - x1;
		y = vy[k] - y1;
		z = vz[k] - z1;

		vLeft = (x - RadZero) / pixelSize;
		vRight = (x + RadZero) / pixelSize;
		ix_left = (int)(floorf(vLeft));
		ix_right = (int)(ceilf(vRight));
		if (ix_left < 0) ix_left = 0;
		if (ix_right > nx - 1) ix_right = nx - 1;

		vLeft = (y - RadZero) / pixelSize;
		vRight = (y + RadZero) / pixelSize;
		iy_left = (int)(floorf(vLeft));
		iy_right = (int)(ceilf(vRight));
		if (iy_left < 0) iy_left = 0;
		if (iy_right > ny - 1) iy_right = ny - 1;

		vLeft = (z - RadZero) / pixelSize;
		vRight = (z + RadZero) / pixelSize;
		iz_left = (int)(floorf(vLeft));
		iz_right = (int)(ceilf(vRight));
		if (iz_left < 0) iz_left = 0;
		if (iz_right > nz - 1) iz_right = nz - 1;

		idx = ix_right - ix_left + 1;
		idy = iy_right - iy_left + 1;
		idz = iz_right - iz_left + 1;

		for (int ig = 0; ig < nGauss; ig++) {

			for (int i = 0; i < idx; i++) {
				ss = float(ix_left + i)*pixelSize - x;
				vX[i] = getValue(ss, indA, ig);
			}
			for (int i = 0; i < idy; i++) {
				ss = float(iy_left + i)*pixelSize - y;
				vY[i] = getValue(ss, indA, ig);
			}
			for (int i = 0; i < idz; i++) {
				ss = float(iz_left + i)*pixelSize - z;
				vZ[i] = getValue(ss, indA, ig);
			}

			for (int i = 0; i < idz; i++) {
				for (int j = 0; j < idy; j++) {
					ippsMulC_32f(vX, vZ[i] * vY[j], vT, idx);
					ippsAdd_32f_I(vT, vm + ((iz_left + i)*ny + iy_left + j)*nx + ix_left, idx);
				}
			}
		}

	}

	ippsFree(vX); vX = nullptr;
	ippsFree(vY); vY = nullptr;
	ippsFree(vZ); vZ = nullptr;
	ippsFree(vT); vT = nullptr;


	return 0;
}



int ListAtoms::createCurves(FILE *flog, float resolution) {
	int iu;
	Ipp64f vAbsMax, vAB;
	Ipp64f s1, s2, s3, s4;
	Ipp64f x1, x2, x3;
	Ipp32f* vg;

	res = resolution;
	aR = 2.772588f / (res*res);

	for (int i = 0; i < 100; i++) {
		countAtom[i] = 0;
	}

	for (int i = 0; i < nAtoms; i++) {
		countAtom[index[i]]++;
	}

	numDifAtoms = 0;

	for (int i = 0; i < 100; i++) {
		if (countAtom[i] > 0) numDifAtoms++;
	}

	fprintf(flog, "Number of different atoms: %i\n", numDifAtoms);

	indexAtom = (int*)malloc(sizeof(int)*numDifAtoms);
	indexD = (int *)malloc(sizeof(int)*nAtoms);
		
	iu = 0;
	for (int i = 0; i < 100; i++) {
		if (countAtom[i] == 0) continue;
		fprintf(flog, "%.2s: %i\n", atomNames5 + 2 * i, countAtom[i]);
		indexAtom[iu] = i;
		iu++;
	}

	for (int i = 0; i < nAtoms; i++) {
		for (int j = 0; j < numDifAtoms; j++) {
			if (index[i] != indexAtom[j])continue;
			indexD[i] = j;
			break;
		}
	}

	vGValues = (Ipp32f **)malloc(sizeof(Ipp32f*)*numDifAtoms);

	vGxyz = (Ipp32f **)malloc(sizeof(Ipp32f*)*numDifAtoms);

	for (int i = 0; i < numDifAtoms; i++) {
		vGValues[i] = ippsMalloc_32f(nCurves);
		vGxyz[i] = ippsMalloc_32f(nCurves*nGauss);
	}

	vRmax = ippsMalloc_32f(numDifAtoms);
	vdR = ippsMalloc_32f(numDifAtoms);
	vVmax = ippsMalloc_32f(numDifAtoms);

	vAbsMax = 1.e9;

	for (int i = 0; i < numDifAtoms; i++) {
		x1 = 0.0;
		vAB = valueAB(x1, i);
		if (vAbsMax > vAB) vAbsMax = vAB;
	}

	s4 = curvePrecision * vAbsMax;

	RadZero = 0.0;

	for (int i = 0; i < numDifAtoms; i++) {

		x1 = 0.0;
		s1 = valueAB(x1, i);
		vVmax[i] = s1;
		x2 = 2.0;
		s2 = valueAB(x2, i);
		while (s2 > s4) {
			x1 = x2;
			s1 = s2;
			x2 *= 2.0;
			s2 = valueAB(x2, i);
		}
		while (x2 - x1 > 1e-6) {
			x3 = 0.5*(x1 + x2);
			s3 = valueAB(x3, i);
			if (s3 > s4) {
				x1 = x3;
				s1 = s3;
			}
			else {
				x2 = x3;
				s2 = s3;
			}
		}
		vRmax[i] = (float)x3;
		
		vdR[i] = vRmax[i] / float(nCurves);

		if (RadZero < vRmax[i]) RadZero = vRmax[i];

		for (int k = 0; k < nGauss; k++) {
			vg = vGxyz[i] + k * nCurves;
			for (int j = 0; j < nCurves; j++) {
				x1 = vdR[i] * double(j);
				vg[j] = (float)valueXYZ(x1, i, k);
			}
		}


	}

	fprintf(flog, "Max zero radius: %f\n", RadZero);



	return 0;
}

Ipp64f ListAtoms::valueAB(Ipp64f x, int ind) {
	Ipp64f s1, s2, s3, cR, w, qR;
	s3 = 0.0;
	for (int j = 0; j < nGauss; j++) {
		s1 = 4.0*pi / (B5values[4 * indexAtom[ind] + j]);
		qR = s1 * pi;
		cR = aR * qR / (aR + qR);
		w = sqrt(pi / (aR + qR));
		s2 = A5values[4 * indexAtom[ind] + j] * s1*sqrt(s1) *w*exp(-cR * x*x);
		s3 += s2;
	}

	return s3;
}

Ipp64f ListAtoms::valueXYZ(Ipp64f x, int ind, int ig) {
	Ipp64f s1, s2, qR, cR, w;

	//s1 = 4.0*pi / (B5values[4 * indexAtom[ind] + ig] + Br);
	//s2 = cbrt(A5values[4 * indexAtom[ind] + ig]) * sqrt(s1) *exp(-s1 * pi*x*x);

	s1 = 4.0*pi / (B5values[4 * indexAtom[ind] + ig]);

	qR = s1 * pi;
	cR = aR * qR / (aR + qR);
	w = sqrt(pi / (aR + qR));

	s2 = cbrt(A5values[4 * indexAtom[ind] + ig] * w) * sqrt(s1) *exp(-cR * x*x);

	return s2;
}

