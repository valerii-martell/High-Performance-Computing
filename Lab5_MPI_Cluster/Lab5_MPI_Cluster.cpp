#define _GNU_SOURCE
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdbool.h>
/* ���� ����, �������� �� ������� */
struct my_matrix 		/* ���������, �� ����� ������� */
{
	int rows; 			/* ʳ������ ����� � ������� */
	int cols; 			/* ʳ������ �������� � ������� */
	double* data; 		/* �������� �� �������, ���������� � ��������� �����*/
};
MPI_Datatype node_mat_t; /* ��� ��� �������������� ������ ���� �������������
�� ������� �������*/
MPI_Datatype mat_col; 	/* ��� ��� �������� ������� ������� */
void inline evaluate(int);	/*�������, �� ������ ��������� ������� */
double func(double X, double Y); /* ������� �������� ����� ������� ��������������
 ������� � ����� X,Y*/
void init_grid(); 		/* �������� �������� ���� ������������� � ������ */
void fatal_error(const char* message, int errorcode);	/* ������� ��������� �������
 �����-������ */
struct my_matrix* matrix_alloc(int rows, int cols, double initial);	/* ������� ���
 �������� ����� �� ����������� �������*/
void matrix_print(const char* filename, struct my_matrix* mat);	/* ������� ���
��������� ������� � ����*/
struct my_matrix* read_matrix(const char* filename);			/* ������� ���
 ���������� ������� �� �����*/
 /* ************************************************************************ */
 /* ���� ������ */
int np; 				/* ʳ������ ������� */
int rank; 			/* ���� ������� � MPI_COMM_WORLD */
int nodes_width; 		/* ������ ������� ������� */
int node_X; 			/* ���������� � � ������� ������� */
int node_Y; 			/* ���������� � � ������� ������� */
int total_width;		 /* ������ ������� ����� ���� �������������*/
int points_per_node; 		/* ������ ������� ������� ����� ���� �������������,
�� ������� �� ����� ������*/
const char* inpfileX = �inpX.csv�; /* ��� ����� � ����� ��������� ������������
���������� ���� ������������� */
const char* inpfileY = �inpY.csv�; /* ��� ����� � ����� ��������� ����������
���������� ���� ������������� */
const char* inpfileInit = �inpInit.csv�;/*��� ����� � ����� ��������� �������� �� �������
 ����� ��� ���� �������������*/
const char* outpfile = �outp.csv�;/* ��� ����� ���� ���������� ��������� */
my_matrix* allcoords_X;  /* ������� �������� �������������� ��������� ����� ����
 ������������� */
my_matrix* allcoords_Y;  /* ������� �������� ������������ ��������� ����� ����
������������� */
my_matrix* all_grid_init;  /*�� �������� �� ������� ����� �������*/
my_matrix* total_mat;  /* ����� ���� �������������, �� ���������� � ������ 0*/
my_matrix* grid_init;   /*������� ���������� ���� ��� �������*/
my_matrix* coords_X;  /* ������� ������� �������� ��������� ���� �������������,
��� ���������� � ������� � �������*/
my_matrix* coords_Y;  /* ������� ������� �������� ��������� ���� �������������,
��� ���������� � ������� � �������*/
my_matrix* node_mat;  /* ������� ���� �������������, ��� ���������� � ������� �
�������*/
my_matrix* f_xy;  	/* ����� ������� ������� � ������ ���� �������������*/:
MPI_Status stat;  	/* �����, � ��� ����������� ������ ������� ���������� */
double epsilon = 0.01;  /* ������� ������*/
double omega = 0.4;   	/* ���������� ����������*/
int max_iter = 1000000; /*����������� ������� �������� � ������� ������� ������� */
double h;		 /* ���� ���� ������������� */
bool* local_stop;	/* ����� ��� ��������������� ��������� ���� ������� */
int* node_iter;
short node_stop;	/* ������ ���������� ��������� � ������*/
short global_stop;	/* ��������� ������ ���������� ���������*/
double* left_col;	/* �������, �� ������ ���� �������� � ����� �� ������� �������
 �� ��� �������� */
double* right_col;
double* top_row;	/* �����, �� ������ ���� �������� � ��������� �� ��������
������� �� ��� �������� */
double* bot_row;
/* ************************************************************************ */
/* ������� �������� */
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv); 			/* ���������� ���������� */
	MPI_Comm_size(MPI_COMM_WORLD, &np); /* �������� ����� MPI_COMM_WORLD */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* �������� ���� ������� � MPI_COMM_WORLD*/
	nodes_width = sqrt(np);	 	/* �������� ������ ������� �������*/
	node_X = rank % nodes_width; 	/* �������� ������������� ���������� ������� �
 ������� */
	node_Y = rank / nodes_width; 	/* �������� ����������� ���������� ������� �
 ������� */
	if (rank == 0)
	{
		allcoords_X = read_matrix(inpfileX);		/* ������� ���������� */
		allcoords_Y = read_matrix(inpfileY);		/* ������� ���������� */
		all_grid_init = read_matrix(inpfileInit);	/* ������� ���������� */
		total_width = allcoords_X->rows; 	/* �������� ������ ���� ������������� */
	}
	MPI_Bcast(&total_width, 1, MPI_INT, 0, MPI_COMM_WORLD); *�������������
		������ ���� �������������* /
		h = 1.0 / total_width;		/* ��������� ���� ���� ������������� */
	points_per_node = total_width / nodes_width;
	/* ��������� ���, �� ������� ������ �������, ��� ���� ������������ */
	MPI_Type_vector(points_per_node, points_per_node, total_width, MPI_DOUBLE, &node_mat_t);
	MPI_Type_commit(&node_mat_t);	 /* �������� ��� */
	coords_X = matrix_alloc(points_per_node, points_per_node, 0.0);
	coords_Y = matrix_alloc(points_per_node, points_per_node, 0.0);
	grid_init = matrix_alloc(points_per_node, points_per_node, 0.0);
	if (rank == 0)
	{	/* ������ 0 ������� ������ ���������� ���� �������������*/
		for (int i = 0; i < nodes_width; i++)
		{
			for (int j = 0; j < nodes_width; j++)
			{
				if (!(i == 0 && j == 0))
				{
					MPI_Send(allcoords_X->data + i * total_width * points_per_node + j * points_per_node
						, 1, node_mat_t, i * nodes_width + j, 0, MPI_COMM_WORLD);
					MPI_Send(allcoords_Y->data + i * total_width * points_per_node + j * points_per_node
						, 1, node_mat_t, i * nodes_width + j, 0, MPI_COMM_WORLD);
					MPI_Send(all_grid_init->data + i * total_width * points_per_node + j * points_per_node
						, 1, node_mat_t, i * nodes_width + j, 0, MPI_COMM_WORLD);
				}
			}
		}
		for (int i = 0; i < points_per_node; i++)
		{
			for (int j = 0; j < points_per_node; j++)
			{
				coords_X->data[i * points_per_node + j] = allcoords_X->data[i * total_width + j];
				coords_Y->data[i * points_per_node + j] = allcoords_Y->data[i * total_width + j];
				grid_init->data[i * points_per_node + j] = all_grid_init->data[i * total_width + j];
			}
		}
	}
	else
	{		/* ���� ������� �� ��������� */
		my_matrix* temp_X = matrix_alloc(total_width, points_per_node, 0);
		my_matrix* temp_Y = matrix_alloc(total_width, points_per_node, 0);
		my_matrix* temp_init = matrix_alloc(total_width, points_per_node, 0);
		MPI_Recv(temp_X->data, 1, node_mat_t, 0, 0, MPI_COMM_WORLD, &stat);
		MPI_Recv(temp_Y->data, 1, node_mat_t, 0, 0, MPI_COMM_WORLD, &stat);
		MPI_Recv(temp_init->data, 1, node_mat_t, 0, 0, MPI_COMM_WORLD, &stat);
		for (int i = 0; i < points_per_node; i++)
		{
			for (int j = 0; j < points_per_node; j++)
			{
				coords_X->data[i * points_per_node + j] = temp_X->data[i * total_width + j];
				coords_Y->data[i * points_per_node + j] = temp_Y->data[i * total_width + j];
				grid_init->data[i * points_per_node + j] = temp_init->data[i * total_width + j];
			}
		}
	}
	/* ���������� ����� �������� ������� ���������*/
	local_stop = (bool*)malloc(points_per_node * points_per_node * sizeof(bool));
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		local_stop[i] = false;
	}
	for (int i = 0; i < points_per_node; i++)	/* ������� ����� ����������� ���������*/
	{
		if (node_Y == 0)
		{
			local_stop[i] = true;
		}
		if (node_X == 0)
		{
			local_stop[i * points_per_node] = true;
		}
		if (node_X == nodes_width - 1)
		{
			local_stop[(i + 1) * points_per_node - 1] = true;
		}
		if (node_Y == nodes_width - 1)
		{
			local_stop[points_per_node * (points_per_node - 1) + i] = true;
		}
	}
	/* ���������� � ������������ �������� �������� ����������� �������� */
	node_mat = matrix_alloc(points_per_node, points_per_node, 0);
	/* ������������ �������� �������� ��������� ����� ���� ������������� */
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		node_mat->data[i] = grid_init->data[i];
	}
	/* �������� ����� ������� ������� � ��������� ������ ���� �������������*/
	f_xy = matrix_alloc(points_per_node, points_per_node, 0.0);
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		f_xy->data[i] = func(coords_X->data[i], coords_Y->data[i]);
	}
	/*��������� ���� ��� ������� �� ��� ��������*/
	/* ���, �� ������� �������� �������*/
	MPI_Type_vector(points_per_node, 1, points_per_node, MPI_DOUBLE, &mat_col);
	MPI_Type_commit(&mat_col);			/* �������� ��� */
	left_col = (double*)malloc(points_per_node * points_per_node * sizeof(double));
	right_col = (double*)malloc(points_per_node * points_per_node * sizeof(double));
	top_row = (double*)malloc(points_per_node * points_per_node * sizeof(double));
	bot_row = (double*)malloc(points_per_node * points_per_node * sizeof(double));
	node_iter = (int*)malloc(points_per_node * points_per_node * sizeof(int));
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		node_iter[i] = 0;
	}
	/* ������� �������� */
	do
	{
		/* ���� �� ��������� � ����*/
		if (node_X != 0)
		{
			/* �������� �������� � ����� ������� */
			MPI_Recv(left_col, 1, mat_col, rank - 1, rank - 1, MPI_COMM_WORLD, &stat);
		}
		if (node_X != nodes_width - 1)
		{
			/* ³������ �������� �� ������� ������� */
			MPI_Send(node_mat->data + points_per_node - 1, 1, mat_col,
				rank + 1, rank, MPI_COMM_WORLD);
			/* �������� �������� � ������� ������� */
			MPI_Recv(right_col, 1, mat_col, rank + 1, rank + 1, MPI_COMM_WORLD, &stat);
		}
		if (node_X != 0)
		{
			/* ³������ �������� �� ����� ������� */
			MPI_Send(node_mat->data, 1, mat_col, rank - 1, rank, MPI_COMM_WORLD);
		}
		if (node_Y != 0)
		{
			/* �������� ������ � ��������� ������� */
			MPI_Recv(top_row, points_per_node, MPI_DOUBLE, rank - nodes_width,
				rank - nodes_width, MPI_COMM_WORLD, &stat);
		}
		if (node_Y != nodes_width - 1)
		{
			/* ³������ ������ �� �������� ������� */
			MPI_Send(node_mat->data + (points_per_node - 1) * points_per_node, points_per_node
				, MPI_DOUBLE, rank + nodes_width, rank, MPI_COMM_WORLD);
			/* �������� ������ � �������� ������� */
			MPI_Recv(bot_row, points_per_node, MPI_DOUBLE, rank + nodes_width,
				rank + nodes_width, MPI_COMM_WORLD, &stat);
		}
		if (node_Y != 0)
		{
			/* ³������ ������ �� ��������� ������� */
			MPI_Send(node_mat->data, points_per_node, MPI_DOUBLE,
				rank - nodes_width, rank, MPI_COMM_WORLD);
		}
		/* ������ �������, �� ��������� ���������� �������*/
		for (int i = 0; i < points_per_node; i++)
		{
			for (int j = 0; j < points_per_node; j++)
			{
				evaluate(i * points_per_node + j);
			}
		}
		/*���� �� ����� ���� ������������� � ������ ��������� ���������� �
		 ������������ ������ ���������� ��������� � ������*/
		node_stop = 1;
		for (int i = 0; i < points_per_node * points_per_node; i++)
		{
			if (!local_stop[i])
			{
				node_stop = 0;
				break;
			}
		}
		global_stop = 1;
		/*���������� �� �� ������� �������� ����������*/
		MPI_Reduce(&node_stop, &global_stop, 1, MPI_SHORT, MPI_BAND, 0, MPI_COMM_WORLD);
		/*���� �� - �� �������� ������ ��� ����������*/
		MPI_Bcast(&global_stop, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
	} while (!global_stop);
	/* ������� ��������� */
	total_mat = matrix_alloc(total_width, total_width, 0.0);
	if (rank == 0)
	{
		for (int i = 0; i < nodes_width; i++)
		{
			for (int j = 0; j < nodes_width; j++)
			{
				if (!(i == 0 && j == 0))
				{
					my_matrix* temp = matrix_alloc(points_per_node, points_per_node, 0.0);
					/* �������� ������� ������� ���� ������������� � ������� ������� */
					MPI_Recv(temp->data, points_per_node * points_per_node,
						MPI_DOUBLE, i * nodes_width + j, i * nodes_width + j, MPI_COMM_WORLD, &stat);
					/* �� ����������� �� � ����� ������� */
					for (int k = 0; k < points_per_node; k++)
					{
						for (int l = 0; l < points_per_node; l++)
						{
							total_mat->data[i * total_width * points_per_node +
								j * points_per_node + k * total_width + l]
								= temp->data[k * points_per_node + l];
						}
					}
				}
			}
		}
		for (int i = 0; i < points_per_node; i++)
		{
			for (int j = 0; j < points_per_node; j++)
			{
				total_mat->data[i * total_width + j] = node_mat->data[i * points_per_node + j];
			}
		}
	}
	else
	{
		MPI_Send(node_mat->data, points_per_node * points_per_node
			, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	if (rank == 0)
	{
		matrix_print(outpfile, total_mat);;
	}
	MPI_Finalize();
}
/* ************************************************************************* */
/* ������� ������� */
double func(double X, double Y)
{
	return X * Y;
}
/*�������, �� ������ ��������� ������� */
void inline evaluate(int index)
{
	/*���� ��������� ������� �������, �� �� ������� ��������� ���������� */
	if (local_stop[index])
	{
		return;
	}
	/*���� ����� ���� ����������� �� ������� �� ���������������
	 �������, �� ������� ����� �������� � ������ ��������*/
	double left = index % points_per_node == 0
		? left_col[index] : node_mat->data[index - 1];
	double right = (index + 1) % points_per_node == 0
		? right_col[index + 1 - points_per_node] : node_mat->data[index + 1];
	double top = index / points_per_node == 0 ? top_row[index % points_per_node]
		: node_mat->data[index - points_per_node];
	double bottom = index / points_per_node == points_per_node - 1
		? bot_row[index % points_per_node]
		: node_mat->data[index + points_per_node];
	/*������ ����������*/
	double old_node = node_mat->data[index];
	double LU = old_node - (left + right + top + bottom - 4 * old_node) / (h * h);
	node_mat->data[index] = old_node - omega * h * h / 4 * (LU - f_xy->data[index]);
	/*���� ��������� ������� ��� ���������� ����������� ������� ��������,
	�� ��������� ���������� � ����� ���� ����*/
	if (fabs(node_mat->data[index] - old_node) <=
		fa bs(node_mat->data[index] * epsilon) || ++node_iter[index] > max_iter)
	{
		local_stop[index] = true;
	}
}
/* �������� �������� ���� ������������� � ������ */
void init_grid()
{
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		node_mat->data[i] = grid_init->data[i];
	}
}
void fatal_error(const char* message, int errorcode)
{
	printf(�fatal error : code % d, % s\n�, errorcode, message);
	fflush(stdout);
	MPI_Abort(MPI_COMM_WORLD, errorcode);
}
struct my_matrix* matrix_alloc(int rows, int cols, double initial)
{
	struct my_matrix* result = (my_matrix*)malloc(sizeof(struct my_matrix));
	result->rows = rows;
	result->cols = cols;
	result->data = (double*)malloc(sizeof(double) * rows * cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			result->data[i * cols + j] = initial;
		}
	}
	return result;
}
void matrix_print(const char* filename, struct my_matrix* mat)
{
	FILE* f = fopen(filename, �w�);
	if (f == NULL)
	{
		fatal_error(�cant write to file�, 2);
	}
	for (int i = 0; i < mat->rows; i++)
	{
		for (int j = 0; j < mat->cols; j++)
		{
			fprintf(f, � % lf �, mat->data[i * mat->cols + j]);
		}
		fprintf(f, �\n�);
	}
	fclose(f);
}
struct my_matrix* read_matrix(const char* filename)
{
	FILE* mat_file = fopen(filename, �r�);
	if (mat_file == NULL)
	{
		fatal_error(�can�t open matrix file�, 1);
	}
	int rows;
	int cols;
	fscanf(mat_file, � % d % d�, &rows, &cols);
	struct my_matrix* result = matrix_alloc(rows, cols, 0.0);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			fscanf(mat_file, � % lf�, &result->data[i * cols + j]);
		}
	}
	fclose(mat_file);
	return result;
}
