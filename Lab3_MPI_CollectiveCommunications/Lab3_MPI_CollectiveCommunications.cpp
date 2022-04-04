#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <time.h>

// It is necessary to find the determinant of a matrix with a given accuracy using LU decomposition

// Messages tags
const int TAG_LINE = 10;

// Timer
clock_t start, endclock;

// Files
const char* input_file = "input.txt";
const char* output_file = "output.txt";

const double epsilon = 0.00001;

// Calculate line
void calcLine(double* working, double* self, int n, int k) {
	double l = k[self] / k[working];
	for (int i = k; i < n; i++) {
		i[self] -= l * i[working];
	}
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double* ma = new double;
	int n;

	if (rank == 0)
	{
		start = clock();

		FILE* input;
		fopen_s(&input, input_file, "r");

		if (!input)
		{
			printf("Failed to open the file\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		fscanf_s(input, "%d", &n);

		ma = new double[n * n];

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				fscanf_s(input, "%lf", &ma[i * n + j]);
			}
		}
		fclose(input);
	}

	// Send matrix size to all tasks
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank != 0) double* ma = new double[n * n];

	// WORK ONLY IF n MULTUPLE TO np
	int groups = n / np;
	double* lines = new double[n * groups];

	for (int i = 0; i < groups; i++) {
		MPI_Scatter(ma + i * n * np, n, MPI_DOUBLE, lines + n * i, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	double* line = new double[n];

	for (int k = 0; k < n - 1; k++) {
		int group = k / np;
		int owner = k % np;
		int nextGroup = (k + 1) / np;
		int nextOwner = (k + 1) % np;
		int start = (group == groups - 1 ? owner + 1 : 0);
		if (rank == 0) {
			MPI_Request* requests = new MPI_Request[np - 1];
			for (int i = 0; i < np - 1; i++)
				requests[i] = MPI_REQUEST_NULL;
			for (int i = start; i < np; i++) {
				if (i != owner) {
					if (i == 0) {
						std::copy(ma + k + k * n, ma + n + k * n, line + k);
					}
					else {
						MPI_Isend(ma + k + k * n, n - k, MPI_DOUBLE, i, TAG_LINE, MPI_COMM_WORLD, &requests[i - 1]);
					}
				}
			}
			MPI_Waitall(np - 1, requests, MPI_STATUSES_IGNORE);
		}

		if (rank >= start) {
			if (rank == owner) {
				std::copy(lines + k + group * n, lines + n + group * n, line + k);
			}
			else if (rank != 0) {
				MPI_Recv(line + k, n - k, MPI_DOUBLE, 0, TAG_LINE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			for (int i = 0; i < groups; i++) {
				if (!(rank == owner && i == group)) {
					calcLine(line, lines + i * n, n, k);
				}
			}
		}
		if (nextOwner == 0 && rank == 0) {
			std::copy(lines + nextGroup * n, lines + n + nextGroup * n, ma + (k + 1) * n);
		}
		else {
			if (rank == nextOwner) {
				MPI_Send(lines + k + nextGroup * n, n - k, MPI_DOUBLE, 0, TAG_LINE, MPI_COMM_WORLD);
			}
			if (rank == 0) {
				MPI_Recv(ma + k + (k + 1) * n, n - k, MPI_DOUBLE, nextOwner, TAG_LINE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}

	if (rank == 0) {
		double det = 1;
		for (int i = 0; i < n; i++) {
			det *= ma[i + i * n];
		}

		FILE* output;
		fopen_s(&output, output_file, "w");

		if (!output)
		{
			printf("No output file!\n");
			MPI_Abort(MPI_COMM_WORLD, 3);
		}

		fprintf(output, "%24.16lf", det);
		fclose(output);

		endclock = clock();
		printf("Time: %20.16lf second(s)\n", ((double)endclock - start) / ((double)CLOCKS_PER_SEC));

	}

	MPI_Finalize();
	return 0;
}
