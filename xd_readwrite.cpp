#include "xdock_header.h"

int getMRCinfo(FILE *fl, char *fileName, MRCinfo *minfo) {
	char header[1030], mapString[4], labels[800];

	int icount, mode;
	int mx, my, mz, mapC, mapR, mapS;
	int iSpaceGroupNumber;
	
	FILE *fi;

	fprintf(fl, "\nFile: \"%s\"\n", fileName);
	fi = fopen(fileName, "rb");
	if (fi == nullptr) {
		fprintf(fl, "Error: cannot open the file.\n");
		printf("Error cannot open file \"%s\"\n", fileName);
		return -1;
	}

	icount = (int)fread(header, sizeof(char), 1024, fi);
	if (icount != 1024) {
		printf("Only %i bytes have been read", icount);
		return -2;
	}

	mapString[0] = header[208];
	mapString[1] = header[209];
	mapString[2] = header[210];
	mapString[3] = header[211];

	if (mapString[0] != 'M' || mapString[1] != 'A' || mapString[2] != 'P' || mapString[3] != ' ') {
		printf("Error: not MRC file.\n");
		fprintf(fl, "Error: not MRC file.\n");
		return -3;
	}
		
	minfo->nx = *(int *)(header + 0);
	minfo->ny = *(int *)(header + 4);
	minfo->nz = *(int *)(header + 8);
	
	mode = *(int *)(header + 12);

	if (mode != 2) {
		printf("Error: we can process only 32-bit signed real numbers");
		fprintf(fl, "Error: we can process only 32-bit signed real numbers");
		return -4;
	}

	minfo->nx_start = *(int *)(header + 16);
	minfo->ny_start = *(int *)(header + 20);
	minfo->nz_start = *(int *)(header + 24);
	mx = *(int *)(header + 28);
	my = *(int *)(header + 32);
	mz = *(int *)(header + 36);
	minfo->cellA[0] = *(float *)(header + 40);
	minfo->cellA[1] = *(float *)(header + 44);
	minfo->cellA[2] = *(float *)(header + 48);
	minfo->cellB[0] = *(float *)(header + 52);
	minfo->cellB[1] = *(float *)(header + 56);
	minfo->cellB[2] = *(float *)(header + 60);

	mapC = *(int *)(header + 64);
	mapR = *(int *)(header + 68);
	mapS = *(int *)(header + 72);

	minfo->densityMin = *(float *)(header + 76);
	minfo->densityMax = *(float *)(header + 80);
	minfo->densityMean = *(float *)(header + 84);

	iSpaceGroupNumber = *(int *)(header + 88);

	minfo->sizeOfExtendedHeader = *(int *)(header + 92);

	minfo->nVersion = *(int *)(header + 108);

	minfo->origin[0] = *(float *)(header + 196);
	minfo->origin[1] = *(float *)(header + 200);
	minfo->origin[2] = *(float *)(header + 204);

	minfo->rms = *(float *)(header + 216);

	strncpy(labels, (char *)header + 224, 800);

	if (header[212] == 0x44 && header[213] == 0x44) {
		minfo->isSwap = false;
	}
	else {
		minfo->isSwap = true;
	}
	

	fprintf(fl, "Number of columns, rows, sections: %i x %i x %i\n", minfo->nx, minfo->ny, minfo->nz);
	fprintf(fl, "Location of first column:          %i\n", minfo->nx_start);
	fprintf(fl, "Location of first row:             %i\n", minfo->ny_start);
	fprintf(fl, "Location of first section:         %i\n", minfo->nz_start);
	fprintf(fl, "Sampling x, y, z:                  %i x %i x %i\n", mx, my, mz);
	fprintf(fl, "Cell size (in A):                  %f x %f x %f\n", minfo->cellA[0], minfo->cellA[1], minfo->cellA[2]);
	fprintf(fl, "Cell size (in degrees):            %f x %f x %f\n", minfo->cellB[0], minfo->cellB[1], minfo->cellB[2]);
	fprintf(fl, "Axes for columns, rows, sections:  %i x %i x %i\n", mapC, mapR, mapS);
	fprintf(fl, "Density (min):                     %f\n", minfo->densityMin);
	fprintf(fl, "Density (max):                     %f\n", minfo->densityMax);
	fprintf(fl, "Density (mean):                    %f\n", minfo->densityMean);
	fprintf(fl, "Density (rms):                     %f\n", minfo->rms);
	fprintf(fl, "Origin x, y, z:                    %f x %f x %f\n", minfo->origin[0], minfo->origin[1], minfo->origin[2]);
	fprintf(fl, "Space group number:                %i\n", iSpaceGroupNumber);
	fprintf(fl, "Size of extended header:           %i\n", minfo->sizeOfExtendedHeader);
	fprintf(fl, "MRC version:                       %i\n\n", minfo->nVersion);

	minfo->origin[0] += (float(minfo->nx_start)*minfo->cellA[0] / float(minfo->nx));
	minfo->origin[1] += (float(minfo->ny_start)*minfo->cellA[1] / float(minfo->ny));
	minfo->origin[2] += (float(minfo->nz_start)*minfo->cellA[2] / float(minfo->nz));

	fprintf(fl, "Origin (absolute) x, y, z:         %f x %f x %f\n", minfo->origin[0], minfo->origin[1], minfo->origin[2]);

	
	fclose(fi); fi = nullptr;

	return 0;
}


int writeMRC(FILE *fl, char *fileName, Ipp32f *v, int nx, int ny, int nz, float *orig, float pixelSize) {

	FILE *fo;
	Ipp32f vMin, vMax, vMean, vStd;
	char header[1200];
	int nt;
	nt = nx * ny * nz;

	ippsZero_8u((Ipp8u*)header, 1024);
	
	ippsMeanStdDev_32f(v, nt, &vMean, &vStd, IppHintAlgorithm::ippAlgHintAccurate);
	ippsMinMax_32f(v, nt, &vMin, &vMax);

	*(Ipp32s*)(header + 0) = nx;
	*(Ipp32s*)(header + 4) = ny;
	*(Ipp32s*)(header + 8) = nz;
	*(Ipp32s*)(header + 12) = 2;
	*(Ipp32s*)(header + 28) = nx;
	*(Ipp32s*)(header + 32) = ny;
	*(Ipp32s*)(header + 36) = nz;

	*(Ipp32f*)(header + 40) = pixelSize * float(nx);
	*(Ipp32f*)(header + 44) = pixelSize * float(ny);
	*(Ipp32f*)(header + 48) = pixelSize * float(nz);

	*(Ipp32f*)(header + 52) = 90.0f;
	*(Ipp32f*)(header + 56) = 90.0f;
	*(Ipp32f*)(header + 60) = 90.0f;

	*(Ipp32s*)(header + 64) = 1;
	*(Ipp32s*)(header + 68) = 2;
	*(Ipp32s*)(header + 72) = 3;

	*(Ipp32f*)(header + 76) = vMin;
	*(Ipp32f*)(header + 80) = vMax;
	*(Ipp32f*)(header + 84) = vMean;

	*(Ipp32s*)(header + 88) = 1;
	*(Ipp32s*)(header + 92) = 0;

	*(Ipp32s*)(header + 108) = 20140;

	*(Ipp32f*)(header + 196) = orig[0];
	*(Ipp32f*)(header + 200) = orig[1];
	*(Ipp32f*)(header + 204) = orig[2];

	header[208] = 'M';
	header[209] = 'A';
	header[210] = 'P';
	header[211] = ' ';

	header[212] = 0x44;
	header[213] = 0x44;

	*(Ipp32f*)(header + 216) = vStd;
	*(Ipp32s*)(header + 220) = 0;


	fprintf(fl, "Output MRC file: \"%s\".\n", fileName);
	fo = fopen(fileName, "wb");
	if (fo == nullptr) {
		fprintf(fl, "Error: cannot open output file \"%s\".\n", fileName);
		printf("Error: cannot open output file \"%s\".\n", fileName);
		return -2;

	}
	fwrite(header, sizeof(char), 1024, fo);
	fwrite(v, sizeof(float), nt, fo);

	fclose(fo);
	fprintf(fl, "Done.\n");
	
	return 0;
}
