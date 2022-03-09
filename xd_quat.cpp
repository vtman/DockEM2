#include "xdock_header.h"

int quat2mat(quat a, Ipp64f *m) {
	m[0] = a.q[0] * a.q[0] + a.q[1] * a.q[1] - a.q[2] * a.q[2] - a.q[3] * a.q[3];
	m[4] = a.q[0] * a.q[0] - a.q[1] * a.q[1] + a.q[2] * a.q[2] - a.q[3] * a.q[3];
	m[8] = a.q[0] * a.q[0] - a.q[1] * a.q[1] - a.q[2] * a.q[2] + a.q[3] * a.q[3];

	m[1] = 2.0 * (a.q[1] * a.q[2] - a.q[0] * a.q[3]);
	m[2] = 2.0 * (a.q[1] * a.q[3] + a.q[0] * a.q[2]);
	m[3] = 2.0 * (a.q[1] * a.q[2] + a.q[0] * a.q[3]);
	m[5] = 2.0 * (a.q[2] * a.q[3] - a.q[0] * a.q[1]);
	m[6] = 2.0 * (a.q[1] * a.q[3] - a.q[0] * a.q[2]);
	m[7] = 2.0 * (a.q[2] * a.q[3] + a.q[0] * a.q[1]);
	return 0;
}

double findAngle(quat a) {
	double s;
	s = acos(a.q[0]) * 2.0 * pi180r;
	return s;
}

double distQ2(quat a, quat b) {
	double d, ang;
	d = a.q[0] * b.q[0] + a.q[1] * b.q[1] + a.q[2] * b.q[2] + a.q[3] * b.q[3];
	ang = acos(2.0 * d * d - 1.0)*pi180r;

	return ang;
}

quat quatBetween(quat a, quat b) {
	quat c;
	double A[16];

	A[0] = a.q[0]; A[1] = a.q[1]; A[2] = a.q[2]; A[3] = a.q[3];
	A[4] = -a.q[1]; A[5] = a.q[0]; A[6] = a.q[3]; A[7] = -a.q[2];
	A[8] = -a.q[2]; A[9] = -a.q[3]; A[10] = a.q[0]; A[11] = a.q[1];
	A[12] = -a.q[3]; A[13] = a.q[2]; A[14] = -a.q[1]; A[15] = a.q[0];

	c.q[0] = A[0] * b.q[0] + A[1] * b.q[1] + A[2] * b.q[2] + A[3] * b.q[3];
	c.q[1] = A[4] * b.q[0] + A[5] * b.q[1] + A[6] * b.q[2] + A[7] * b.q[3];
	c.q[2] = A[8] * b.q[0] + A[9] * b.q[1] + A[10] * b.q[2] + A[11] * b.q[3];
	c.q[3] = A[12] * b.q[0] + A[13] * b.q[1] + A[14] * b.q[2] + A[15] * b.q[3];

	return c;
}


quat quatProd(quat a, quat b) {
	quat c;
	double A[16];

	A[0] = a.q[0]; A[1] = -a.q[1]; A[2] = -a.q[2]; A[3] = -a.q[3];
	A[4] = a.q[1]; A[5] = a.q[0]; A[6] = -a.q[3]; A[7] = a.q[2];
	A[8] = a.q[2]; A[9] = a.q[3]; A[10] = a.q[0]; A[11] = -a.q[1];
	A[12] = a.q[3]; A[13] = -a.q[2]; A[14] = a.q[1]; A[15] = a.q[0];

	c.q[0] = A[0] * b.q[0] + A[1] * b.q[1] + A[2] * b.q[2] + A[3] * b.q[3];
	c.q[1] = A[4] * b.q[0] + A[5] * b.q[1] + A[6] * b.q[2] + A[7] * b.q[3];
	c.q[2] = A[8] * b.q[0] + A[9] * b.q[1] + A[10] * b.q[2] + A[11] * b.q[3];
	c.q[3] = A[12] * b.q[0] + A[13] * b.q[1] + A[14] * b.q[2] + A[15] * b.q[3];

	return c;

}


int XData::getNumberOfGroups() {
	char fileQuat[1000], fileNameQ[1000];
	double* qD;
	long long fSize;
	FILE* fa;

	switch (ip->nShapeROI) {
	case SHAPE_DIFF: sprintf(fileNameQ, "quat.bin"); break;
	case SHAPE_XY: sprintf(fileNameQ, "quatXY.bin"); break;
	case SHAPE_XZ: sprintf(fileNameQ, "quatXZ.bin"); break;
	case SHAPE_YZ: sprintf(fileNameQ, "quatYZ.bin"); break;
	case SHAPE_XYZ: sprintf(fileNameQ, "quatXYZ.bin"); break;
	}

	if (strlen(ip->systemFolder) == 0) {
		sprintf(fileQuat, fileNameQ);
	}
	else {
		sprintf(fileQuat, "%s/%s", ip->systemFolder, fileNameQ);
	}

	fprintf(flog, "Binary file for optimal quaternions: \"%s\"\n", fileQuat);
	fa = fopen(fileQuat, "rb");
	if (fa == nullptr) {
		fprintf(flog, "Error: cannot open the file.\n");
		return -1;
	}

	fseeko(fa, 0, SEEK_END);
	fSize = ftello(fa);
	fprintf(flog, "File size: %lld\n", fSize);

	nGroupsAvailable = fSize / (40);
	fprintf(flog, "Number of groups available: %i\n", nGroupsAvailable);

	qD = (double*)malloc(sizeof(double) * 5*nGroupsAvailable);
	vAngleLevel = (double*)malloc(sizeof(double) * nGroupsAvailable);
	vInputQuat = (quat*)malloc(sizeof(quat) * nGroupsAvailable);
	matRot = ippsMalloc_64f(9 * nGroupsAvailable);

	rewind(fa);
	fread(qD, sizeof(double), nGroupsAvailable * 5, fa);
	fclose(fa); fa = nullptr;

	for (int i = 0; i < nGroupsAvailable; i++) {
		memcpy(vInputQuat + i, qD + 5 * i, sizeof(double) * 4);
		vAngleLevel[i] = qD[5 * i + 4];
	}
	free(qD); qD = nullptr;
	

	nGroupEstimate = 0;
	nGroupMain = 0;

	while (vAngleLevel[nGroupMain] > ip->angularResolution) {
		nGroupMain++;
	}

	while (vAngleLevel[nGroupEstimate] > ip->angleForEstimation) {
		nGroupEstimate++;
	}

	fprintf(flog, "Groups of rotations for main processing: %i\n", nGroupMain);
	fprintf(flog, "Groups of rotations to estimate Std: %i\n", nGroupEstimate);

	fprintf(flog, "\n");

	return 0;
}
