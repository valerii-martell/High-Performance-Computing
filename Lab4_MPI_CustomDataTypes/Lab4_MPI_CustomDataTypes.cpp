#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// It is necessary to solve the system of linear algebraic equations with a given accuracy using the Jacobi method

clock_t start;
clock_t end;

// Accuracy
const double epsilon = 0.0001;

// Messages tags
const int BLOCK_TAG = 11;
const int FREE_VECT_TAG = 12;
const int RESUIDE_NORM_PART_TAG = 13;
const int RESUIDE_NORM_TAG = 14;
const int NEWX_TAG = 15;
const int OLDX_TAG = 16;


struct BaseVector {
	int size;
	double data[1]; //Base for dynamic array
};
typedef struct BaseVector Vector;

struct BaseMatrix {
	int rows;
	int cols;
	double data[1]; //Base for dynamic array
};
typedef struct BaseMatrix Matrix;

Matrix* readMatrix(const char* path) {
	FILE* file;
	fopen_s(&file, path, "r");

	int size[2];

	if (!file)
	{
		printf("Failed to open the file\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	for (int i = 0; i < 2; i++)
		fscanf_s(file, "%d", &size[i]);

	//Dynamically create array
	Matrix* res;
	res = (Matrix*)malloc(sizeof(res) + sizeof(double) * (size[0] * size[1] - 1)); 
	res->rows = size[0];
	res->cols = size[1];

	for (int i = 0; i < size[0]; i++) {
		for (int j = 0; j < size[1]; j++) {
			fscanf_s(file, "%lf", &res->data[i * size[0] + j]);
		}
	}

	fclose(file);

	return res;
}

Vector* readVector(const char* path) {
	FILE* file;
	fopen_s(&file, path, "r");

	int size;

	if (!file)
	{
		printf("Failed to open the file\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	fscanf_s(file, "%d", &size);

	Vector* res;
	res = (Vector*)malloc(sizeof(res) + sizeof(double) * (size - 1)); //Dynamically create array
	res->size = size;

	for (int i = 0; i < size; i++) {
		fscanf_s(file, "%lf", &res->data[i]);
	}

	fclose(file);

	return res;
}

void writeVector(const char* path, Vector* vec) {
	FILE* file;
	fopen_s(&file, path, "w");

	if (!file)
	{
		printf("Failed to open the file\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	for (int i = 0; i < vec->size; i++) {
		fprintf(file, "%lf ", vec->data[i]);
	}
	fclose(file);
}

void Build_derived_Matrix_type(Matrix matrix, int len, MPI_Datatype* MPI_NODE_ptr)
{
	int block_lengths[3];
	MPI_Aint displacements[3];

	MPI_Datatype typelists[3];
	MPI_Aint start_address;
	MPI_Aint address;
	//fill block lengths
	block_lengths[0] = block_lengths[1] = 1;
	block_lengths[2] = len;
	//Fill Typelists
	typelists[0] = MPI_INT;
	typelists[1] = MPI_INT;
	typelists[2] = MPI_DOUBLE;

	//Calculate Displacements
	displacements[0] = 0;
	MPI_Get_address(&(matrix.rows), &start_address);
	MPI_Get_address(&(matrix.cols), &address);
	displacements[1] = address - start_address;
	MPI_Get_address(&(matrix.data), &address);
	displacements[2] = address - start_address;
	//Create and Commit new mpi type
	MPI_Type_create_struct(3, block_lengths, displacements, typelists, MPI_NODE_ptr);
	MPI_Type_commit(MPI_NODE_ptr);
}

void Build_derived_Vector_type(Vector vector, int len, MPI_Datatype* MPI_NODE_ptr)
{
	int block_lengths[2];
	MPI_Aint displacements[2];

	MPI_Datatype typelists[2];
	MPI_Aint start_address;
	MPI_Aint address;
	//fill block lengths
	block_lengths[0] = 1;
	block_lengths[1] = len;
	//Fill Typelists
	typelists[0] = MPI_INT;
	typelists[1] = MPI_DOUBLE;

	//Calculate Displacements
	displacements[0] = 0;
	MPI_Get_address(&(vector.size), &start_address);
	MPI_Get_address(&(vector.data), &address);
	displacements[1] = address - start_address;
	//Create and Commit new mpi type
	MPI_Type_create_struct(2, block_lengths, displacements, typelists, MPI_NODE_ptr);
	MPI_Type_commit(MPI_NODE_ptr);
}

void jacobi_iteration(
	int start_row, 				// The number of the first row of the matrix part
	Matrix* mat_A_part, 		// Part of the rows of the coefficient matrix
	Vector* b, 					// Vector free members
	Vector* vec_prev_x, 		// Preliminary vector
	Vector* vec_next_x_part, 	// Part of the vector of the next approximation (set in function)
	double* residue_norm_part) 	// The value of the norm in the previous step
{
	// The accumulator value of the norm of this part of the calculations
	double my_residue_norm_part = 0.0;
	// Elemental computation of part of the vector of the following approximation
	for (int i = 0; i < vec_next_x_part->size; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < mat_A_part->cols; j++)
		{
			if (i + start_row != j)
			{
				sum += mat_A_part->data[i * mat_A_part->cols + j] * vec_prev_x->data[j];
			}
		}
		sum = b->data[i + start_row] - sum;
		vec_next_x_part->data[i] = sum / mat_A_part->data[i * mat_A_part->cols + i + start_row];
		// Calculate the norm in the previous step
		sum = -sum + mat_A_part->data[i * mat_A_part->cols + i + start_row] *
			vec_prev_x->data[i + start_row];
		my_residue_norm_part += sum * sum;
	}
	*residue_norm_part = my_residue_norm_part;
}

// Main function 
int main(int argc, char* argv[])
{
	const char* input_file_MA = "MA.txt";
	const char* input_file_b = "b.txt";
	const char* output_file_x = "x.txt";

	Matrix* MA = new Matrix();
	Vector* b = new Vector();

	MPI_Init(&argc, &argv);

	int np, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int N;
	if (rank == 0)
	{
		start = clock();	
		
		MA = readMatrix(input_file_MA);
		
		b = readVector(input_file_b);
		
		if (MA->rows != MA->cols) {
			printf("Matrix is not square\n");
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
		if (MA->rows != b->size) {
			printf("Matrix and vector doesn't match\n");
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
		N = b->size;
		//for (int i = 0; i < N * N; i++) 
		//	printf("%f ", MA->data[i]);
		//writeVector(output_file_x, b);
	}
	
	// Sending to all matrix and vector dimension problems
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// Allocate memory to store vector of free members
	if (rank != 0)
	{
		b = (Vector*)malloc(sizeof(b) + sizeof(double) * (N - 1)); //Dynamically create array
		b->size = N;
		for (int i = 0; i < N; i++) 
			b->data[i] = 0.0;
	}

	// Calculate the part of vectors and matrix that will be stored in each problem, we assume that N = k * np. 
	// Allocate memory to store parts of vectors and matrices in each task and set their initial values
	int part = N / np;
	Matrix* MAh;
	MAh = (Matrix*)malloc(sizeof(MAh) + sizeof(double) * (N * part - 1)); //Dynamically create array
	MAh->rows = part;
	MAh->cols = N;	
	for (int i = 0; i < part * N; i++)
			MAh->data[i] = 0.0;
	
	Vector* oldx;
	oldx = (Vector*)malloc(sizeof(oldx) + sizeof(double) * (N - 1)); //Dynamically create array
	oldx->size = N;
	for (int i = 0; i < N; i++)
		oldx->data[i] = 0.0;
	
	Vector* newx;
	newx = (Vector*)malloc(sizeof(newx) + sizeof(double) * (part - 1)); //Dynamically create array
	newx->size = part;
	for (int i = 0; i < part; i++)
		newx->data[i] = 0.0;
	
	MPI_Datatype MPI_Matrix_part;
	MPI_Datatype MPI_Vector;
	MPI_Datatype MPI_Vector_part;
	Build_derived_Matrix_type(*MAh, N*part, &MPI_Matrix_part); //Builds derived type
	Build_derived_Vector_type(*oldx, N, &MPI_Vector); //Builds derived type
	Build_derived_Vector_type(*newx, part, &MPI_Vector_part); //Builds derived type
	
	// Divide the original matrix MA into parts by part rows and send parts to all problems. 
	// Free up memory allocated for the MA matrix.
	if (rank == 0)
	{
		int offset = N*part;
		for (int i = 1; i < np; i++)
		{
			for (int j = 0; j < N * part; j++)
				MAh->data[j] = MA->data[offset + j];
			offset += N * part;
			MPI_Send(MAh, 1, MPI_Matrix_part, i, BLOCK_TAG, MPI_COMM_WORLD);
		}
		for (int j = 0; j < N * part; j++)
			MAh->data[j] = MA->data[j];

		// Deallocation
		//free(MA);
		//delete MA;
	}
	else
	{
		MPI_Recv(MAh, 1, MPI_Matrix_part, 0, BLOCK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);		
	}

	// Sending a vector of free members and freeing memory from it in the main task
	MPI_Bcast(b, 1, MPI_Vector, 0, MPI_COMM_WORLD);
	
	// Calculate the norm of the vector of free terms in problem 0 and send it to the workers
	double b_norm = 0.0;
	if (rank == 0)
	{
		for (int i = 0; i < b->size; i++)
		{
			b_norm = pow(b->data[i], 2);
		}
		b_norm = sqrt(b_norm);
	}
	MPI_Bcast(&b_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// The value of the stop iteration criterion
	double last_stop_criteria;
	// Jacobi main iteration cycle
	for (int i = 0; i < 1000; i++)
	{
		double residue_norm_part = 0;
		double residue_norm = 0;
		jacobi_iteration(rank * part, MAh, b, oldx, newx, &residue_norm_part);

		MPI_Allreduce(&residue_norm_part, &residue_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		residue_norm = sqrt(residue_norm);
		// Check the iteration stop criteria.
		last_stop_criteria = residue_norm / b_norm;
		if (last_stop_criteria < epsilon)
		{
			break;
		}
		// Receive values of current approximation of vector unknown
		MPI_Allgather(newx->data, part, MPI_DOUBLE, oldx->data, part, MPI_DOUBLE, MPI_COMM_WORLD);
	}


	// Output result
	if (rank == 0)
	{
		writeVector(output_file_x, oldx);
		end = clock();
		printf("Time: %20.16lf second(s)\n", ((double)end - start) / ((double)CLOCKS_PER_SEC));
	}

	MPI_Finalize();
	return 0;
}