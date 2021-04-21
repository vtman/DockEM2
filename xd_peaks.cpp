#include "xdock_header.h"

int XData::findPeaks() {
	int iq, kq, jq, kqq, iqq, jqq, np, np2;
	Ipp64f sk, sj, si;
	unsigned int u1, u2;

	Ipp32s* vPeak, * vindex;
	Ipp32u* uIndex, * uTemp;
	Ipp32f* v_loc, * v_loc2, *wCC;
	Ipp64f* vAng;

	quat* vQ;
	quat qa, qb, qc, qres;

	bool t;

	npeaks = 0;
	vPeak = (Ipp32s*)pmcTemp[0];

	for (int k = 1; k < nuz - 1; k++) {
		for (int i = 1; i < nuy - 1; i++) {
			for (int j = 1; j < nux - 1; j++) {
				v_loc = vOut + (k * nuy + i) * nux + j;

				t = false;
				for (int kk = -1; kk <= 1; kk++) {
					for (int ii = -1; ii <= 1; ii++) {
						for (int jj = -1; jj <= 1; jj++) {
							v_loc2 = v_loc + (kk * nuy + ii) * nux + jj;
							if (v_loc2 == v_loc) continue;
							if (v_loc[0] < v_loc2[0]) {
								t = true;
								break;
							}
						}
						if (t) break;
					}
					if (t) break;
				}
				if (t) continue;

				vPeak[npeaks] = (k * nuy + i) * nux + j;
				npeaks++;
			}
		}
	}

	np2 = npeaks * npeaks;

	vQ = (quat*)malloc(sizeof(quat) * npeaks);

	vAng = ippsMalloc_64f(np2);

	wCC = ippsMalloc_32f(npeaks);
	uIndex = ippsMalloc_32u(npeaks);

	uTemp = ippsMalloc_32u(npeaks);
	vindex = ippsMalloc_32s(npeaks);

	for (int i = 0; i < npeaks; i++) {
		wCC[i] = vOut[vPeak[i]];
		uIndex[i] = uOut[vPeak[i]];
	}

	ippsSortIndexDescend_32f_I(wCC, vindex, npeaks);

	for (int i = 0; i < npeaks; i++) {
		uTemp[i] = uIndex[vindex[i]];
	}
	ippsCopy_32s((Ipp32s*)uTemp, (Ipp32s*)uIndex, npeaks);

	for (int i = 0; i < npeaks; i++) {
		uTemp[i] = *((Ipp32u*)vPeak + vindex[i]);
	}
	ippsCopy_32s((Ipp32s*)uTemp, (Ipp32s*)vPeak, npeaks);

	fprintf(flog, "\nNumber of local peaks: %i\n\n", npeaks);

	fprintf(flog, "    #        CC      ix    iy    iz   ix(T) iy(T) iz(T)          fx        fy        fz          q_0         q_1         q_2         q_4\n");
	fprintf(flog, "----------------------------------------------------------------------------------------------------------------------------------------\n");

	if (npeaks >= 9999) npeaks = 9999;

	int nn;
	nn = ip->outputReduction;

	for (int i = 0; i < npeaks; i++) {
		u1 = uIndex[i];
		u2 = u1 & 31;
		u1 >>= 5;
		qa = vLocalQuat[u2];
		qc = qa;
		qc.q[1] *= -1.0;
		qc.q[2] *= -1.0;
		qc.q[3] *= -1.0;
		qb = vInputQuat[u1];

		qres = quatProd(qc, qb);
		vQ[i] = qres;

		kq = vPeak[i] / (nux * nuy);
		iq = (vPeak[i] / nux) % nuy;
		jq = vPeak[i] % nux;
		if (wCC[i] >= ip->minCCpeak) {
			fprintf(flog, "%5i  %8.5f | %5i %5i %5i | %5i %5i %5i | %9.5f %9.5f %9.5f | %10.7f  %10.7f  %10.7f  %10.7f\n", i + 1, wCC[i], jq, iq, kq,
				nn * jq + iCornerOrig[0], nn * iq + iCornerOrig[1], nn * kq + iCornerOrig[2],
				targetOrig[0] + float(jq) * pixTOut, targetOrig[1] + float(iq) * pixTOut, targetOrig[2] + float(kq) * pixTOut, qres.q[0], qres.q[1], qres.q[2], qres.q[3]);
		}
		if (i < ip->nPDBout) {
			sprintf(outputFileName, "%s/%s_%05i.pdb", ip->outputFolder, ip->outputPrefix, i + 1);
			writeRotatedPDB(outputFileName, cgx, cgy, cgz, targetOrig[0] + float(jq) * pixTOut - cgx, targetOrig[1] + float(iq) * pixTOut - cgy,
				targetOrig[2] + float(kq) * pixTOut - cgz, qres);
		}
	}
	fprintf(flog, "----------------------------------------------------------------------------------------------------------------------------------------\n");

	np = myImin(ip->nRelativeInfo, npeaks);

	for (int i = 0; i < np; i++) {
		kqq = vPeak[i] / (nux * nuy);
		iqq = (vPeak[i] / nux) % nuy;
		jqq = vPeak[i] % nux;
		vAng[i * np + i] = 0.0;
		for (int j = i + 1; j < np; j++) {
			kq = vPeak[j] / (nux * nuy) - kqq;
			iq = (vPeak[j] / nux) % nuy - iqq;
			jq = vPeak[j] % nux - jqq;

			sk = double(kq);
			si = double(iq);
			sj = double(jq);

			vAng[i * np + j] = pixTOut * sqrt(sk * sk + si * si + sj * sj);
			vAng[j * np + i] = vAng[i * np + j];
		}
	}

	fprintf(flog, "\nDistance between centres (in A)\n");

	for (int i = 0; i < (np + 1) * 9; i++) {
		fprintf(flog, "-");
	}
	fprintf(flog, "\n      ");

	for (int i = 0; i < np; i++) {
		fprintf(flog, "%9i", i + 1);
	}
	fprintf(flog, "\n");
	for (int i = 0; i < (np + 1) * 9; i++) {
		fprintf(flog, "-");
	}
	fprintf(flog, "\n");

	for (int i = 0; i < np; i++) {
		fprintf(flog, "%7i |", i + 1);
		for (int j = 0; j < np; j++) {
			fprintf(flog, "%8.4f ", vAng[np * i + j]);
		}
		fprintf(flog, "\n");
	}

	for (int i = 0; i < (np + 1) * 9; i++) {
		fprintf(flog, "-");
	}
	fprintf(flog, "\n");



	for (int i = 0; i < np; i++) {
		vAng[i * np + i] = 0.0;
		for (int j = i + 1; j < np; j++) {
			vAng[i * np + j] = distQ2(vQ[i], vQ[j]);
			vAng[j * np + i] = vAng[i * np + j];
		}
	}

	fprintf(flog, "\nAngular distance between rotations (in degrees)\n");


	for (int i = 0; i < (np + 1) * 7; i++) {
		fprintf(flog, "-");
	}
	fprintf(flog, "\n      ");

	for (int i = 0; i < np; i++) {
		fprintf(flog, "%7i", i + 1);
	}
	fprintf(flog, "\n");
	for (int i = 0; i < (np + 1) * 7; i++) {
		fprintf(flog, "-");
	}
	fprintf(flog, "\n");

	for (int i = 0; i < np; i++) {
		fprintf(flog, "%5i |", i + 1);
		for (int j = 0; j < np; j++) {
			fprintf(flog, "%6.2f ", vAng[np * i + j]);
		}
		fprintf(flog, "\n");
	}

	for (int i = 0; i < (np + 1) * 7; i++) {
		fprintf(flog, "-");
	}
	fprintf(flog, "\n");







	free(vQ); vQ = nullptr;
	ippsFree(vAng); vAng = nullptr;

	ippsFree(wCC); wCC = nullptr;
	ippsFree(uIndex); uIndex = nullptr;
	ippsFree(uTemp); uTemp = nullptr;
	ippsFree(vindex); vindex = nullptr;

	return 0;
}
