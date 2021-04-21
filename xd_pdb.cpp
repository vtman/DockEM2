#include "xdock_header.h"

int XData::readPDB() {
	FILE* fi;
	char qw[10], * vi;
	long fS;
	int i1, i2, na;

	fi = fopen(ip->inputSearchFile, "rb");
	if (fi == nullptr) {
		fprintf(flog, "Error: cannot open search file \"%s\".\n", ip->inputSearchFile);
		printf("Error: cannot open search file \"%s\".\n", ip->inputSearchFile);
		return -1;
	}

	fseek(fi, 0, SEEK_END);
	fS = ftell(fi);
	rewind(fi);

	vi = (char*)malloc(sizeof(char) * fS);

	fread(vi, sizeof(char), fS, fi);

	fclose(fi); fi = nullptr;

	pdbRecords = 0;
	for (int i = 0; i < fS; i++) {
		if (vi[i] == '\n') pdbRecords++;
	}
	pdbData = (char*)malloc(sizeof(char) * 81 * pdbRecords);
	int irec, jj;
	irec = 0;
	jj = 0;
	for (int i = 0; i < fS; i++) {
		if (vi[i] == '\n' || vi[i] == '\r') {
			for (; jj < 80; jj++) {
				pdbData[irec * 81 + jj] = ' ';
			}
			pdbData[irec * 81 + 80] = '\n';
			irec++;
			if (irec == pdbRecords) break;
			jj = 0;
			i++;
			while (vi[i] == '\n' || vi[i] == '\r') {
				i++;
			}

			i--;
			continue;
		}
		pdbData[irec * 81 + jj] = vi[i];
		jj++;
	}

	pO1 = -1;
	pO2 = -1;
	pO3 = -1;

	for (int i = 0; i < pdbRecords; i++) {
		if (strncmp(pdbData + 81 * i, "ORIGX1", 6) == 0) pO1 = i;
		if (strncmp(pdbData + 81 * i, "ORIGX2", 6) == 0) pO2 = i;
		if (strncmp(pdbData + 81 * i, "ORIGX3", 6) == 0) pO3 = i;
	}

	/*FILE *fop;
	fop = fopen("C:\\Temp2\\xDock\\Nov\\in\\out.pdb", "wb");
	fwrite(pdbData, sizeof(char), 81 * pdbRecords, fop);
	fclose(fop);
	*/


	na = 0;

	for (i1 = 0; i1 < fS; i1++) {
		if (vi[i1] == '\n') continue;
		if (vi[i1] == '\r') continue;
		for (i2 = i1 + 1; i2 < fS; i2++) {
			if (vi[i2] == '\n' || vi[i2] == '\r') break;
		}
		if (strncmp(vi + i1, "ATOM  ", 6) == 0 || strncmp(vi + i1, "HETATM", 6) == 0) {
			na++;
			continue;
		}
		i1 = i2;
	}


	fprintf(flog, "Number of atoms: %i\n", na);

	la->nAtoms = na;
	la->index = (int*)malloc(sizeof(int) * na);
	la->vx = (float*)malloc(sizeof(float) * na);
	la->vy = (float*)malloc(sizeof(float) * na);
	la->vz = (float*)malloc(sizeof(float) * na);

	na = 0;
	qw[8] = '\0';

	for (i1 = 0; i1 < fS; i1++) {
		if (vi[i1] == '\n') continue;
		if (vi[i1] == '\r') continue;
		for (i2 = i1 + 1; i2 < fS; i2++) {
			if (vi[i2] == '\n' || vi[i2] == '\r') break;
		}
		if (strncmp(vi + i1, "ATOM  ", 6) == 0 || strncmp(vi + i1, "HETATM", 6) == 0) {
			strncpy(qw, vi + i1 + 30, 8); la->vx[na] = (float)atof(qw);
			strncpy(qw, vi + i1 + 38, 8); la->vy[na] = (float)atof(qw);
			strncpy(qw, vi + i1 + 46, 8); la->vz[na] = (float)atof(qw);
			if (getAtomIndex(flog, vi + i1 + 76, la->index + na) != 0) return -2;
			na++;
			continue;
		}
		i1 = i2;
	}

	free(vi); vi = nullptr;

	return 0;
}

int XData::writeRotatedPDB(char* fPDB, float cx, float cy, float cz, float shiftX, float shiftY, float shiftZ, quat qres) {
	FILE* fo;
	Ipp64f matABC[9], xyz_in[3], xyz_out[3];
	Ipp64f sum;
	char* vcO, qw[80];

	vcO = (char*)malloc(sizeof(char) * 81 * pdbRecords);

	memcpy(vcO, pdbData, 81 * pdbRecords);

	fo = fopen(fPDB, "wb");
	if (fo == nullptr) {
		fprintf(flog, "Error: cannot write to \"%s\"\n", fPDB);
		printf("Error: cannot write to \"%s\"\n", fPDB);
		return -1;
	}

	quat2mat(qres, matABC);


	for (int i = 0; i < pdbRecords; i++) {
		if (strncmp(pdbData + 81 * i, "ATOM  ", 6) != 0 && strncmp(pdbData + 81 * i, "HETATM", 6) != 0) continue;

		strncpy(qw, pdbData + 81 * i + 30, 8); xyz_in[0] = atof(qw) - cx;
		strncpy(qw, pdbData + 81 * i + 38, 8); xyz_in[1] = atof(qw) - cy;
		strncpy(qw, pdbData + 81 * i + 46, 8); xyz_in[2] = atof(qw) - cz;

		for (int k = 0; k < 3; k++) {
			sum = 0.0;
			for (int j = 0; j < 3; j++) {
				sum += matABC[3 * k + j] * xyz_in[j];
			}
			xyz_out[k] = sum;
		}
		xyz_out[0] += (cx + shiftX);
		xyz_out[1] += (cy + shiftY);
		xyz_out[2] += (cz + shiftZ);
		sprintf(qw, "%8.3f%8.3f%8.3f", xyz_out[0], xyz_out[1], xyz_out[2]);

		strncpy(vcO + 81 * i + 30, qw, 24);
	}

	fwrite(vcO, sizeof(char), pdbRecords * 81, fo);
	free(vcO); vcO = nullptr;

	fclose(fo); fo = nullptr;

	return 0;
}
