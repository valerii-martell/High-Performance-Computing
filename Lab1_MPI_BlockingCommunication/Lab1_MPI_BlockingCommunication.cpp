#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

// It is necessary to calculate the area hyperbolic tangent from a given number with a given accuracy

// Timers
clock_t start;
clock_t end;

// Accuracy
const double epsilon = 1E-11;

// Filenames
const char* input_fl = "input.txt";
const char* output_fl = "output.txt";

// Messages tags
const int valueTag = 10;
const int statusTag = 11;
const int stepTag = 12;
const int resultTag = 13;

// Main function
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	double number;
	bool status = true;

	// Main task
	if (rank == 0) {
		start = clock();

		// Read õ
		FILE* input;
		fopen_s(&input, input_fl, "r");

		if (!input)
		{
			printf("No input file!\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}

		fscanf_s(input, "%lf", &number);
		fclose(input);
		//printf("%d", number);

		// Check the range
		if (fabs(number) >= 1) {
			printf("Argument must be in range (-1;1) !\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}

		double result = number;

		for (int i = 1; i <= np - 1; i++) {
			// Send x to workers
			MPI_Send(&number, 1, MPI_DOUBLE, i, valueTag, MPI_COMM_WORLD);
		}

		int step = 1;

		while (status) {
			for (int i = 1; i <= np - 1; i++) {
				// Send current step of calculation to workers
				MPI_Send(&step, 1, MPI_INT, i, stepTag, MPI_COMM_WORLD);
				step++;
			}

			for (int i = 1; i <= np - 1; i++) {
				double res;
				// Receive results from workers
				MPI_Recv(&res, 1, MPI_DOUBLE, i, resultTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				result += res;

				// Check stop criteria
				if (fabs(res) < epsilon) {
					status = false;
				}
			}

			for (int i = 1; i <= np - 1; i++) {
				// Send status about desired accuracy achievment  
				MPI_Send(&status, 1, MPI_C_BOOL, i, statusTag, MPI_COMM_WORLD);
			}

		}

		FILE* output;
		fopen_s(&output, output_fl, "w");

		if (!output)
		{
			printf("No output file!\n");
			MPI_Abort(MPI_COMM_WORLD, 3);
			return 3;
		}
		fprintf(output, "%24.16lf", result);
		fclose(output);

	}
	// Workers 
	else {
		// Receive õ
		MPI_Recv(&number, 1, MPI_DOUBLE, 0, valueTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		while (status) {
			int step;
			// Receive the number of current calculation step for this task
			MPI_Recv(&step, 1, MPI_INT, 0, stepTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			double res = pow(number, ((step * 2) + 1)) / ((step * 2) + 1);

			// Send the result to master task
			MPI_Send(&res, 1, MPI_DOUBLE, 0, resultTag, MPI_COMM_WORLD);
			// Receive the status
			MPI_Recv(&status, 1, MPI_C_BOOL, 0, statusTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	if (rank == 0) {
		end = clock();
		printf("Time: %20.16lf second(s)\n", ((double)end - start) / ((double)CLOCKS_PER_SEC));
	}

	MPI_Finalize();
	return 0;
}
