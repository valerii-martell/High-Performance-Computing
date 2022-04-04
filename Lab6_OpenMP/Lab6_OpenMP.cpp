#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "omp.h"

// It is necessary to find all (both simpleand complex) divisors of a natural number A.
// For this purpose, the number A is divided by all positive integers[1; A].
// The task is to write a function that implements a parallel version of this algorithm.
// The dividers are output using printf().

// Number
const int A = 1000000;
// Number of threads
const int P = 4;

// Main function
int main()
{
	// Timer start
	int start_time = clock();


	// Simple solution. For the case when not necessary to store dividers
	/*
#pragma omp parallel for default(none) shared(A)
	for (int i = 1; i < A; i++) 
	{
		if (A % i == 0) printf("%d ", i);
	}
	*/

	// A more complicated solution. For the case when necessary to save dividers in sorted order

	// Count of dividers
	int count = 0;
	// Dividers flags array
	bool isDiv[A];

#pragma omp parallel num_threads(P) default(none) shared(A, isDiv) reduction(+:count)
	{
		int task_id = omp_get_thread_num();
		int H = A / P;
		int start_index = task_id * H;
		int end_index;
		if (task_id == P - 1) 			// If this is the last task by identifier
			end_index = A - 1; 			// Then it works with the last part of array 
		else							// If this isn't the last task by identifier
			end_index = start_index + H - 1; // Then it works with the fixed (H) part of array 
		if (start_index == 0) start_index=1;
		for (int i = start_index; i <= end_index; i++)
		{
			if (A % i == 0)
			{
				isDiv[i] = true;
				count++;
			}
		}
	}
	
	// Print dividers count
	printf("Count: %d\n", count);

	// Allocate memory for array of dividers 
	int* res;
	res = (int*)malloc(sizeof(int) * count);
	count = 0;
	// Filling the array of dividers
	for (int i = 1; i < A; i++)
	{
		if (isDiv[i] == true)
		{
			res[count] = i;
			count++;
		}
	}
	// Print all dividers
	for (int i = 0; i < count; i++)
	{
		printf("%d ", res[i]);
	}
	//free(res);
	//free(isDiv);
	
	// Timer stop
	int end_time = clock();
	
	// Get the duration and print
	double time = (double)(end_time - start_time) / 1000;
	printf("Time: %f ms.\n", time);

	return 0;
}
