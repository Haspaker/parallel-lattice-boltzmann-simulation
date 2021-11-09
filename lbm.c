// mpicc lbm.c -std=gnu11 -o lbm_out -O3

#include <stdio.h>
#include <stdlib.h> // malloc
#include <string.h> //memcpy
#include <math.h>

/* Use MPI */
#include "mpi.h"

#define RED 0
#define BLACK 1

#define K 9 // number of velocity directions

#define IDX_YX(y, x) ((int) (y*M + x))
#define IDX_YXI(y, x, i) ((int) (i*M*M + y*M + x))
#define GET_MACRO(_1,_2,_3,NAME,...) NAME
#define idx(...) GET_MACRO(__VA_ARGS__, IDX_YXI, IDX_YX)(__VA_ARGS__)

#define OUTPUT_FILE "csim_data.bin"

#define dprintf if(DEBUG)printf

void write_data_to_row(double* A, double* data, int y, int i, int M) {
	for (int x = 0; x < M; x++) {
		A[idx(y, x, i)] = data[x];
	}
}

void write_data_to_col(double* A, double* data, int x, int i, int M) {
	for (int y = 0; y < M; y++) {
		A[idx(y, x, i)] = data[y];
	}
}

void write_row_to_data(double* A, double* data, int y, int i, int M) {
	for (int x = 0; x < M; x++) {
		data[x] = A[idx(y, x, i)];
	}
}

void write_col_to_data(double* A, double* data, int x, int i, int M) {
	for (int y = 0; y < M; y++) {
		data[y] = A[idx(y, x, i)];
	}
}

int cyclic_shift(int px, int shift, int sqrtP) {
	px = (px + shift) % sqrtP;
	if (px < 0) px += sqrtP;
	return px;
}

int get_p(int py, int px, int dy, int dx, int sqrtP) {
	return cyclic_shift(py, dy, sqrtP) * sqrtP + cyclic_shift(px, dx, sqrtP);
}

int main(int argc, char *argv[]) {

    MPI_Status status;
    int ghost_tag = 100;
    int print_tag = 200;
    int P, p;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);


    if (p == 0) {
    	if (argc < 2) printf("Missing length L (position 2)!\n" );
    	if (argc < 3) printf("Missing time T (position 3)!\n" );
    	if (argc < 4) printf("Missing strength G (position 4)!\n" );
    	if (argc < 5) printf("Missing DEBUG flag (position 5)!\n" );
    	if (argc < 6) printf("Missing SAVE flag (position 6)!\n" );
    	if (argc < 6) exit(1);
    }

    int L = atoi(argv[1]);
    int T = atoi(argv[2]);
    double G = atof(argv[3]);
    int DEBUG = atoi(argv[4]);
    int SAVE = atoi(argv[5]);


    if (p == 0) {
    	printf("L = %d\n", L);
    	printf("T = %d\n", T);
    	printf("G = %f\n", G);
    	printf("DEBUG = %d\n", DEBUG);
    	printf("SAVE = %d\n", SAVE);
    }

    if (floor(sqrt(P)) != sqrt(P)) {
        fprintf(stdout, "Number of processors must be a square number.\n");
        exit(1);
    }

    int sqrtP = (int) sqrt(P);

    if (L % sqrtP != 0) {
        fprintf(stdout, "sqrt(P) must evenly divide L.\n");
        exit(1);
    }

    int px = p % sqrtP;
    int py = (p - px) / sqrtP;
    //dprintf("Indx: p:%d px:%d py:%d sqrtP:%d\n", p, px, py, sqrtP);
    int colorx = px % 2 == 0 ? RED : BLACK; 
    int colory = py % 2 == 0 ? RED : BLACK; 
    int M = L/sqrtP + 2;

	double ex[] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
	double ey[] = {0, 0, 1,  0, -1, 1,  1, -1, -1};
	double  w[] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};

	int updown[] = {2, 5, 6, 4, 7, 8};
	int rightleft[] = {1, 5, 8, 3, 6, 7};
	int corners[] = {5, 6, 7, 8};

	double* f = (double *) malloc(M*M*K*sizeof(double));
	double* rho = (double *) malloc(M*M*sizeof(double));
	double* ux = (double *) malloc(M*M*sizeof(double));
	double* uy = (double *) malloc(M*M*sizeof(double));
	double* psi = (double *) malloc(M*M*sizeof(double));

	double* psisend = (double *) malloc(M*sizeof(double));
	double* psirecv = (double *) malloc(M*sizeof(double));
	double* fsend = (double *) malloc(M*3*sizeof(double));
	double* frecv = (double *) malloc(M*3*sizeof(double));

	srandom(p);
	for (int x = 1; x < M-1; x++) {
		for (int y = 1; y < M-1; y++) {
			double r = ((double) random())/RAND_MAX;
			double initial_rho = 2.0 + r/100.0;
			for (int i = 0; i < K; i++) {
				f[idx(y, x, i)] = w[i] * initial_rho;
			}
		}
	}

	if (SAVE && p == 0) {
		char write_mode = 'w';
		FILE *file = fopen(OUTPUT_FILE, &write_mode);
		fwrite(&L, sizeof(int), 1, file);
		fwrite(&T, sizeof(int), 1, file);
		fwrite(&P, sizeof(int), 1, file);
		fclose (file);
	}

	for (int t = 0; t < T; t++) {
		for (int x = 1; x < M-1; x++) {
			for (int y = 1; y < M-1; y++) {
				rho[idx(y, x)] = 0;
				ux[idx(y, x)] = 0;
				uy[idx(y, x)] = 0;
				for (int i = 0; i < K; i++) {
					rho[idx(y, x)] += f[idx(y, x, i)];
					ux[idx(y, x)] += f[idx(y, x, i)] * ex[i];
					uy[idx(y, x)] += f[idx(y, x, i)] * ey[i];
				}
				ux[idx(y, x)] /= rho[idx(y, x)];
				uy[idx(y, x)] /= rho[idx(y, x)];
				psi[idx(y, x)] = 2*exp(-2/rho[idx(y, x)]);
			}
		}

		if (SAVE && t % 10 == 0) {
			char dummy;
			if (p > 0) {
			    MPI_Recv(&dummy, 1, MPI_BYTE, p-1, print_tag, MPI_COMM_WORLD, &status);
			}
			char append_mode = 'a';

			FILE *file = fopen(OUTPUT_FILE, &append_mode);
			fwrite(rho, sizeof(double), M*M, file);
			fclose (file);

			if (p < P-1) {
			    MPI_Send(&dummy, 1, MPI_BYTE, p+1, print_tag, MPI_COMM_WORLD);
			}
		}

        // up-down
        for (int d = 0; d < 2; d++) {
        	int dir = ey[updown[d*3]];
			if (colory == RED) {
				int send_row_idx = dir < 0 ? 1 : M-2;
				int recv_row_idx = dir < 0 ? 0 : M-1; 
				int partner = get_p(py, px, -dir, 0, sqrtP);
				write_row_to_data(psi, psisend, send_row_idx, 0, M);
				MPI_Send(psisend, M, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				MPI_Recv(psirecv, M, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				write_data_to_row(psi, psirecv, recv_row_idx, 0, M);
			} else {
				int send_row_idx = dir < 0 ? M-2 : 1; 
				int recv_row_idx = dir < 0 ? M-1 : 0; 
				int partner = get_p(py, px, dir, 0, sqrtP);
				write_row_to_data(psi, psisend, send_row_idx, 0, M);
				MPI_Recv(psirecv, M, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				MPI_Send(psisend, M, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				write_data_to_row(psi, psirecv, recv_row_idx, 0, M);
			}
		}

        // right-left
        for (int d = 0; d < 2; d++) {
        	int dir = ex[rightleft[d*3]];
			if (colorx == RED) {
				int send_col_idx = dir < 0 ? 1 : M-2;
				int recv_col_idx = dir < 0 ? 0 : M-1; 
				int partner = get_p(py, px, 0, dir, sqrtP);
				write_col_to_data(psi, psisend, send_col_idx, 0, M);
				MPI_Send(psisend, M, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				MPI_Recv(psirecv, M, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				write_data_to_col(psi, psirecv, recv_col_idx, 0, M);
			} else {
				int send_col_idx = dir < 0 ? M-2 : 1; 
				int recv_col_idx = dir < 0 ? M-1 : 0; 
				int partner = get_p(py, px, 0, -dir, sqrtP);
				write_col_to_data(psi, psisend, send_col_idx, 0, M);
				MPI_Recv(psirecv, M, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				MPI_Send(psisend, M, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				write_data_to_col(psi, psirecv, recv_col_idx, 0, M);
			}
		}

		// corners
        for (int d = 0; d < 4; d++) {
        	int dirx = ex[corners[d]];
        	int diry = ey[corners[d]];
			if (colorx == RED) {
				int send_x = dirx < 0 ? 1 : M-2;
				int recv_x = dirx < 0 ? 0 : M-1;
				int send_y = diry < 0 ? 1 : M-2;
				int recv_y = diry < 0 ? 0 : M-1;
				int partner = get_p(py, px, -diry, dirx, sqrtP);
				MPI_Send(psi + idx(send_y, send_x), 1, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				MPI_Recv(psi + idx(recv_y, recv_x), 1, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
			} else {
				int send_x = dirx < 0 ? M-2 : 1;
				int recv_x = dirx < 0 ? M-1 : 0;
				int send_y = diry < 0 ? M-2 : 1;
				int recv_y = diry < 0 ? M-1 : 0;
				int partner = get_p(py, px, diry, -dirx, sqrtP);
				MPI_Recv(psi + idx(recv_y, recv_x), 1, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				MPI_Send(psi + idx(send_y, send_x), 1, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);			
			}
		}

		for (int x = 1; x < M-1; x++) {
			for (int y = 1; y < M-1; y++) {
				for (int i = 0; i < K; i++) {
					int neighbor_x = x - (int)ex[i];
					int neighbor_y = y - (int)ey[i];
					double Fxi = G * w[i] * psi[idx(y, x)] * psi[idx(neighbor_y, neighbor_x)] * ex[i];
					double Fyi = G * w[i] * psi[idx(y, x)] * psi[idx(neighbor_y, neighbor_x)] * ey[i];
					ux[idx(y, x)] = ux[idx(y, x)] + Fxi / rho[idx(y, x)];
					uy[idx(y, x)] = uy[idx(y, x)] + Fyi / rho[idx(y, x)];
				}

				for (int i = 0; i < K; i++) {
					double edotu = ux[idx(y, x)] * ex[i] + uy[idx(y, x)] * ey[i];
					double udotu = ux[idx(y, x)] * ux[idx(y, x)] + uy[idx(y, x)] * uy[idx(y, x)];
					int stream_x = x + (int)ex[i];
					int stream_y = y + (int)ey[i];
					f[idx(stream_y, stream_x, i)] = w[i] * rho[idx(y, x)] * (1 + 3*edotu + 4.5*edotu*edotu - 1.5*udotu);
				}

			}
		}

		// up-down
        for (int d = 0; d < 2; d++) {
        	int dir = ey[updown[d*3]];
			if (colory == RED) {
				int send_row_idx = dir < 0 ? 0 : M-1;
				int recv_row_idx = dir < 0 ? 1 : M-2; 
				int partner = get_p(py, px, -dir, 0, sqrtP);
				for (int j = 0; j < 3; j++) 
					write_row_to_data(f, fsend+M*j, send_row_idx, updown[j+d*3], M);
				MPI_Send(fsend, M*3, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				MPI_Recv(frecv, M*3, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				for (int j = 0; j < 3; j++)
					write_data_to_row(f, frecv+M*j, recv_row_idx, updown[j+(1-d)*3], M);
			} else {
				int send_row_idx = dir < 0 ? M-1 : 0; 
				int recv_row_idx = dir < 0 ? M-2 : 1; 
				int partner = get_p(py, px, dir, 0, sqrtP);
				for (int j = 0; j < 3; j++) 
					write_row_to_data(f, fsend+M*j, send_row_idx, updown[j+(1-d)*3], M);
				MPI_Recv(frecv, M*3, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				MPI_Send(fsend, M*3, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				for (int j = 0; j < 3; j++)
					write_data_to_row(f, frecv+M*j, recv_row_idx, updown[j+d*3], M);
			}
		}


		// right-left
        for (int d = 0; d < 2; d++) {
        	int dir = ex[rightleft[d*3]];
			if (colorx == RED) {
				int send_col_idx = dir < 0 ? 0 : M-1;
				int recv_col_idx = dir < 0 ? 1 : M-2; 
				int partner = get_p(py, px, 0, dir, sqrtP);
				for (int j = 0; j < 3; j++) 
					write_col_to_data(f, fsend+M*j, send_col_idx, rightleft[j+d*3], M);
				MPI_Send(fsend, M*3, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				MPI_Recv(frecv, M*3, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				for (int j = 0; j < 3; j++) 
					write_data_to_col(f, frecv+M*j, recv_col_idx, rightleft[j+(1-d)*3], M);
			} else {
				int send_col_idx = dir < 0 ? M-1 : 0; 
				int recv_col_idx = dir < 0 ? M-2 : 1; 
				int partner = get_p(py, px, 0, -dir, sqrtP);
				for (int j = 0; j < 3; j++)
					write_col_to_data(f, fsend+M*j, send_col_idx, rightleft[j+(1-d)*3], M);
				MPI_Recv(frecv, M*3, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD, &status);
				MPI_Send(fsend, M*3, MPI_DOUBLE, partner, ghost_tag, MPI_COMM_WORLD);
				for (int j = 0; j < 3; j++) 
					write_data_to_col(f, frecv+M*j, recv_col_idx, rightleft[j+d*3], M);
			}
		}

		dprintf("%d\n", t);

	}

    MPI_Finalize();
    return(0);

}