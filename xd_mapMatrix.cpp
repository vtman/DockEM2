#include "xdock_header.h"

int findBestAll(Ipp32fc* vCCold, Ipp32fc* uQold, Ipp32fc* vCCnew, Ipp32fc* uQnew, int nt) {
	Ipp32f* v1in, * v2in, * v1out, * v2out;
	v1in = (Ipp32f*)vCCold;
	v2in = (Ipp32f*)uQold;
	v1out = (Ipp32f*)vCCnew;
	v2out = (Ipp32f*)uQnew;
	
	for (int i = 0; i < 2*nt; i++) {
		if (v1in[i] > v1out[i]) {
			v1out[i] = v1in[i];
			v2out[i] = v2in[i];
		}
	}

	return 0;
}


int real2compFFT(Ipp32fc* vC, int nx, int ny, int nz) {
	int nn, ixc;
	nn = (nx + 1) / 2 - 1;
	ixc = nx - nn;
	int k1, i1;

	for (int k = 0; k < nz; k++) {
		k1 = nz - k;
		if (k1 == nz) k1 = 0;
		for (int i = 0; i < ny; i++) {
			i1 = ny - i;
			if (i1 == ny) i1 = 0;
			ippsConjFlip_32fc(vC + k1 * nx * ny + i1 * nx + 1, vC + k * nx * ny + i * nx + ixc, nn);
		}
	}
	return 0;
}

int readuceFmap(Ipp32fc *vf, int nx, int ny, int nz, int nScale) {
	int nmx, nmy, nmxy;

	nmx = nx * nScale;
	nmy = ny * nScale;
	nmxy = nmx * nmy;

	for (int k = 1; k < nScale; k++) {
		ippsAdd_32fc_I(vf + (k * nz * nmxy), vf, nz * nmxy);
	}

	for (int k = 0; k < nz; k++) {
		for (int i = 1; i < nScale; i++) {
			ippsAdd_32fc_I(vf + (k * nmxy + i * ny * nmx), vf + (k * nmxy), nmx * ny);
		}
	}

	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < ny; i++) {
			for (int j = 1; j < nScale; j++) {
				ippsAdd_32fc_I(vf + k * nmxy + i * nmx + j * nx, vf + k * nmxy + i * nmx, nx);
			}
		}
	}
	return 0;
}

int mapQMatrix(Ipp32f* vIn, Ipp32f* vOut, Ipp32f* v1, quat* q, int nx, int ny, int nz, int ind){
	int kk, ii, jj;

	switch (ind) {
	case 0: // (j, i, k)
		ippsCopy_32f(vIn, vOut, nx*ny*nz);
		q->q[0] = 1.0; q->q[1] = 0.0; q->q[2] = 0.0; q->q[3] = 0.0;
		break;

	case 1: // (j, -i, -k)
		for (int k = 0; k < nz; k++) {
			ippsCopy_32f(vIn + k * nx * ny, v1 + (nz - 1 - k) * nx * ny, nx * ny);
		}
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				ippsCopy_32f(v1 + (k * ny + i) * nx, vOut + (k * ny + (ny - 1 - i)) * nx, nx);
			}
		}
		q->q[0] = 0.0; q->q[1] = 1.0; q->q[2] = 0.0; q->q[3] = 0.0;
		break;

	case 2: // (-j, i, -k)
		for (int k = 0; k < nz*ny; k++) {
			ippsFlip_32f(vIn + k * nx, v1 + k * nx, nx);
		}
		for (int k = 0; k < nz; k++) {
			ippsCopy_32f(v1 + k * nx * ny, vOut + (nz - 1 - k) * nx * ny, nx * ny);
		}
		q->q[0] = 0.0; q->q[1] = 0.0; q->q[2] = 1.0; q->q[3] = 0.0;
		break;

	case 3: // (-j, -i, k)
		for (int k = 0; k < nz * ny; k++) {
			ippsFlip_32f(vIn + k * nx, v1 + k * nx, nx);
		}
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				ippsCopy_32f(v1 + (k * ny + i) * nx, vOut + (k * ny + (ny - 1 - i)) * nx, nx);
			}
		}
		q->q[0] = 0.0; q->q[1] = 0.0; q->q[2] = 0.0; q->q[3] = 1.0;
		break;

	case 4: // (-k, -i, -j)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = nz - 1 - k;
					ii = ny - 1 - i; 
					kk = nx - 1 - j;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = 0.0; q->q[1] = sqrt(0.5); q->q[2] = 0.0; q->q[3] = -sqrt(0.5);
		break;

	case 5: // (-k, i, j)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = nz - 1 - k;
					ii = i;
					kk = j;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -sqrt(0.5); q->q[1] = 0.0; q->q[2] = sqrt(0.5); q->q[3] = 0.0;
		break;

	case 6: // (k, -i, j)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = k;
					ii = ny - 1 - i;
					kk = j;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = 0.0; q->q[1] = -sqrt(0.5); q->q[2] = 0.0; q->q[3] = -sqrt(0.5);
		break;

	case 7: // (k, i, -j)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = k;
					ii = i;
					kk = nx - 1 -j;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = sqrt(0.5); q->q[1] = 0.0; q->q[2] = sqrt(0.5); q->q[3] = 0.0;
		break;

	case 8: // (k, j, i)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = k;
					ii = j;
					kk = i;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = 0.5; q->q[1] = 0.5; q->q[2] = 0.5; q->q[3] = 0.5;
		break;

	case 9: // (k, -j, -i)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = k;
					ii = nx - 1 - j;
					kk = ny - 1 - i;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -0.5; q->q[1] = 0.5; q->q[2] = -0.5; q->q[3] = 0.5;
		break;

	case 10: // (-k, j, -i)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = nz - 1 - k;
					ii = j;
					kk = ny - 1 - i;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -0.5; q->q[1] = 0.5; q->q[2] = 0.5; q->q[3] = -0.5;
		break;

	case 11: // (-k, -j, i)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = nz - 1 - k;
					ii = nx - 1 - j;
					kk = i;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -0.5; q->q[1] = -0.5; q->q[2] = 0.5; q->q[3] = 0.5;
		break;

	case 12: // (-i, -j, -k)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = ny - 1 - i; 
					ii = nx - 1 - j; 
					kk = nz - 1 - k; 
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = 0.0; q->q[1] = sqrt(0.5); q->q[2] = -sqrt(0.5); q->q[3] = 0.0;
		break;

	case 13: // (-i, j, k)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = ny - 1 - i;
					ii = j;
					kk = k;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -sqrt(0.5); q->q[1] = 0.0; q->q[2] = 0.0; q->q[3] = -sqrt(0.5);
		break;

	case 14: // (i, -j, k)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = i;
					ii = nx - 1 - j;
					kk = k;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = sqrt(0.5); q->q[1] = 0.0; q->q[2] = 0.0; q->q[3] = -sqrt(0.5);
		break;

	case 15: // (i, j, -k)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = i;
					ii = j;
					kk = nz - 1 - k;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = 0.0; q->q[1] = sqrt(0.5); q->q[2] = sqrt(0.5); q->q[3] = 0.0;
		break;

	case 16: // (i, k, j)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = i;
					ii = k;
					kk = j;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -0.5; q->q[1] = 0.5; q->q[2] = 0.5; q->q[3] = 0.5;
		break;

	case 17: // (i, -k, -j)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = i;
					ii = nz - 1 - k;
					kk = nx - 1 - j;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -0.5; q->q[1] = -0.5; q->q[2] = -0.5; q->q[3] = 0.5;
		break;

	case 18: // (-i, k, -j)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = ny - 1 - i;
					ii = k;
					kk = nx - 1 - j;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -0.5; q->q[1] = 0.5; q->q[2] = -0.5; q->q[3] = -0.5;
		break;

	case 19: // (-i, -k, j)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = ny - 1 - i;
					ii = nz - 1 - k;
					kk = j;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -0.5; q->q[1] = -0.5; q->q[2] = 0.5; q->q[3] = -0.5;
		break;

	case 20: // (-j, -k, -i)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = nx - 1 - j;
					ii = nz - 1 - k;
					kk = ny - 1 - i;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = 0.0; q->q[1] = 0.0; q->q[2] = sqrt(0.5); q->q[3] = -sqrt(0.5);
		break;

	case 21: // (-j, k, i)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = nx - 1 - j;
					ii = k;
					kk = i;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = 0.0; q->q[1] = 0.0; q->q[2] = sqrt(0.5); q->q[3] = sqrt(0.5);
		break;

	case 22: // (j, -k, i)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = j;
					ii = nz - 1 - k;
					kk = i;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = -sqrt(0.5); q->q[1] = -sqrt(0.5); q->q[2] = 0.0; q->q[3] = 0.0;
		break;

	case 23: // (j, k, -i)
		for (int k = 0; k < nz; k++) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++) {
					jj = j;
					ii = k;
					kk = ny - 1 - i;
					vOut[(kk * ny + ii) * nx + jj] = vIn[(k * ny + i) * nx + j];
				}
			}
		}
		q->q[0] = sqrt(0.5); q->q[1] = -sqrt(0.5); q->q[2] = 0.0; q->q[3] = 0.0;
		break;

	};
	return 0;
}
