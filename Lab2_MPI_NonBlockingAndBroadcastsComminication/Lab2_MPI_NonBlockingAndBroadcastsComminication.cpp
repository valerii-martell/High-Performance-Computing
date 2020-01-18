#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>

// It is necessary to calculate the value of the integral of the input function 
// in the interval between 0 and 15 with a given accuracy using the Simpson's Rule

// Timers
clock_t start;
clock_t end;

// Filenames
const char* input_fl = "input.txt";
const char* output_fl = "output.txt";

// Function for integration
double function(double x)
{
	return cos(x)*pow(x, 3);
	//return x;
}

// Runge function for accuracy
bool check_Runge(double I2, double I, double epsilon)
{
	return (fabs(I2 - I) / 15) < epsilon;
}

// The Simpson's Rule
double parabola(double start, double finish, double epsilon)
{
	int num_iterations = 1;
	double last_res = 0;
	double res = -1, res1 = -1;
	double h = 0;
	while (!check_Runge(res, last_res, epsilon))
	{
		//printf("%d %lg | ", num_iterations, res);
		num_iterations *= 2;

		last_res = res;
		res = 0;
		res1 = 0;

		h = (finish - start) / num_iterations;
		/*
		for (int i = 0; i < num_iterations - 1; i++)
		{
			if (i != 0) 
			{
				res += function(start + i * h);
			}
			res1 += function(start + i * h + h / 2);
		}
		res = h / 3 * ((function(start) / 2) + res + res1 + (function(finish) / 2));		
		*/
		for (int i = 0; i < num_iterations - 1; i++)
		{
			double a = start + i * h;
			double b = start + i * h + h;
			res += ((b - a) / 6) * (function(a) + 4 * function((a + b) / 2) + function(b));
		}
	}
	return res;
}

// Write one double number in *filename* file
void write_double_to_file(const char* filename, double data)
{
	// Open file for writing
	FILE* fp;
	fopen_s(&fp, filename, "w");
	// Check the file
	if (!fp)
	{
		printf("Failed to open the file\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	fprintf(fp, "%lg\n", data);
	fclose(fp);
}

// Main function
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	

	// Enter data from a file into an array of 3 variables.
	// Occurs in Problem 0.
	// input [0] is the lower bound of integration
	// input [1] is the upper bound of integration
	// input [2] is a valid absolute error
	double input[3];
	if (rank == 0)
	{
		start = clock();
		// Open file input.txt in only-read mode
		FILE* fp;
		fopen_s(&fp, "input.txt", "r");
		// Check the file
		if (!fp)
		{
			printf("Failed to open the file\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		// Read 3 double numbers
		for (int i = 0; i < 3; i++)
			fscanf_s(fp, "%lg", &input[i]);
		fclose(fp);

		//for (int i = 0; i < 3; i++)
		//	printf("%lg ", input[i]);
		/*FILE* inp;
		fopen_s(&inp, input_fl, "r");

		if (!inp)
		{
			printf("No input file!\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}

		/* Зчитування 3 чисел типу double */
		/*for (int i = 0; i < 3; i++)
			fscanf_s(inp, "%lg", &input[i]);
		fclose(inp);
		/*printf("%lf,", input[0]);
		printf("%lf,", input[1]);
		printf("%lf,", input[2]);*/
	}
	// Send the entered data from task 0 to all other tasks
	MPI_Bcast(input, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	double start = input[0];
	double finish = input[1];
	double eps = input[2];
	double step = (finish - start) / np;
	double res = parabola(start + rank * step, start + (rank + 1) * step, eps / np);
	double result = res;
	// Send an intermediate partial result from all tasks except task 0 to task 0
	MPI_Reduce(&res, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	// Writing the result in file "output.txt" 
	if (rank == 0) 
	{
		/*FILE* output;
		fopen_s(&output, output_fl, "w");

		if (!output)
		{
			printf("No output file!\n");
			MPI_Abort(MPI_COMM_WORLD, 3);
			return 3;
		}
		fprintf(output, "%24.16lf", result);
		fclose(output);*/
 
		write_double_to_file("output.txt", result);
		printf("%lf", result);

		end = clock();
		printf("Time: %24.16lf second(s)\n", ((double)end - start) / ((double)CLOCKS_PER_SEC));
	}
	MPI_Finalize();
	return 0;
}
