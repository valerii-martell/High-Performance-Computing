#define _GNU_SOURCE
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdbool.h>
/* Опис типів, структур та функцій */
struct my_matrix 		/* Структура, що описує матриці */
{
	int rows; 			/* Кількість строк у матриці */
	int cols; 			/* Кількість стовпців у матриці */
	double* data; 		/* Покажчик на матрицю, розвернуту у одномірний масив*/
};
MPI_Datatype node_mat_t; /* Тип для розповсюдження частин сітки дискретизації
по матриці процесів*/
MPI_Datatype mat_col; 	/* Тип для передачі стовпця матриці */
void inline evaluate(int);	/*Функція, що містить ітераційні формули */
double func(double X, double Y); /* Повертає значення правої частини диференційного
 рівняння у точці X,Y*/
void init_grid(); 		/* Ініціалізує значення сітки дискретизації у процесі */
void fatal_error(const char* message, int errorcode);	/* Функція виведення помилок
 вводу-виводу */
struct my_matrix* matrix_alloc(int rows, int cols, double initial);	/* Функція для
 виділення пам’яті та ініціалізацію матриці*/
void matrix_print(const char* filename, struct my_matrix* mat);	/* Функція для
виведення матриці у файл*/
struct my_matrix* read_matrix(const char* filename);			/* Функція для
 зчитування матриці із файлу*/
 /* ************************************************************************ */
 /* Опис змінних */
int np; 				/* Кількість процесів */
int rank; 			/* Ранг процесу у MPI_COMM_WORLD */
int nodes_width; 		/* Ширина матриці процесів */
int node_X; 			/* Координата Х у матриці процесів */
int node_Y; 			/* Координата У у матриці процесів */
int total_width;		 /* Ширина матриці вузлів сітки дискретизації*/
int points_per_node; 		/* Ширина частини матриці вузлів сітки дискретизації,
що припадає на кожен процес*/
const char* inpfileX = ”inpX.csv”; /* Ім’я файлу з якого вводяться горизонтальні
координати сітки дискретизації */
const char* inpfileY = ”inpY.csv”; /* Ім’я файлу з якого вводяться вертикальні
координати сітки дискретизації */
const char* inpfileInit = ”inpInit.csv”;/*Ім’я файлу з якого вводяться початков та граничні
 умови для сітки дискретизації*/
const char* outpfile = ”outp.csv”;/* Ім’я файлу куди виводиться результат */
my_matrix* allcoords_X;  /* Матриця реальних горизонтальних координат вузлів сітки
 дискретизації */
my_matrix* allcoords_Y;  /* Матриця реальних вертикальних координат вузлів сітки
дискретизації */
my_matrix* all_grid_init;  /*Всі початкові та граничні умови функції*/
my_matrix* total_mat;  /* Повна сітка дискретизації, що виводиться у процесі 0*/
my_matrix* grid_init;   /*Частина початкових умов для процесу*/
my_matrix* coords_X;  /* Частина матриці реальних координат сітки дискретизації,
яка зберігається у кожному з процесів*/
my_matrix* coords_Y;  /* Частина матриці реальних координат сітки дискретизації,
яка зберігається у кожному з процесів*/
my_matrix* node_mat;  /* Частина сітки дискретизації, яка зберігається у кожному з
процесів*/
my_matrix* f_xy;  	/* Масив значень функції в вузлах сітки дискретизації*/:
MPI_Status stat;  	/* Змінна, у яку повертається статус прийому повідомлень */
double epsilon = 0.01;  /* Точність рішення*/
double omega = 0.4;   	/* Коефіцієнт релаксації*/
int max_iter = 1000000; /*Максимальна кількість ітерацій у випадку нев’язки функції */
double h;		 /* Крок сітки дискретизації */
bool* local_stop;	/* Масив для запам’ятовування локальних умов зупинки */
int* node_iter;
short node_stop;	/* Ознака завершення обчислень у процесі*/
short global_stop;	/* Глобальна ознака завершення обчислень*/
double* left_col;	/* Стовпці, що процес буде приймати з лівого та правого процесів
 під час ітерацій */
double* right_col;
double* top_row;	/* Рядки, що процес буде приймати з верхнього та нижнього
процесів під час ітерацій */
double* bot_row;
/* ************************************************************************ */
/* Головна програма */
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv); 			/* Ініціалізуємо середовище */
	MPI_Comm_size(MPI_COMM_WORLD, &np); /* Отримуємо розмір MPI_COMM_WORLD */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Отримуємо ранг процесу у MPI_COMM_WORLD*/
	nodes_width = sqrt(np);	 	/* Отримуємо ширину матриці процесів*/
	node_X = rank % nodes_width; 	/* Отримуємо горизонтальну координату процесу у
 матриці */
	node_Y = rank / nodes_width; 	/* Отримуємо вертикальну координату процесу у
 матриці */
	if (rank == 0)
	{
		allcoords_X = read_matrix(inpfileX);		/* Вводимо координати */
		allcoords_Y = read_matrix(inpfileY);		/* Вводимо координати */
		all_grid_init = read_matrix(inpfileInit);	/* Вводимо координати */
		total_width = allcoords_X->rows; 	/* Отримуємо ширину сітки дискретизації */
	}
	MPI_Bcast(&total_width, 1, MPI_INT, 0, MPI_COMM_WORLD); *Розповсюджуємо
		ширину сітки дискретизації* /
		h = 1.0 / total_width;		/* Визначаємо крок сітки дискретизації */
	points_per_node = total_width / nodes_width;
	/* Створюємо тип, що відповідає частині матриці, яка буде передаватися */
	MPI_Type_vector(points_per_node, points_per_node, total_width, MPI_DOUBLE, &node_mat_t);
	MPI_Type_commit(&node_mat_t);	 /* Реєструємо тип */
	coords_X = matrix_alloc(points_per_node, points_per_node, 0.0);
	coords_Y = matrix_alloc(points_per_node, points_per_node, 0.0);
	grid_init = matrix_alloc(points_per_node, points_per_node, 0.0);
	if (rank == 0)
	{	/* Процес 0 розсилає реальні координати сітки дискретизації*/
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
	{		/* Інші процеси їх приймають */
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
	/* Ініціалізуємо умови локальної зупинки алгоритму*/
	local_stop = (bool*)malloc(points_per_node * points_per_node * sizeof(bool));
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		local_stop[i] = false;
	}
	for (int i = 0; i < points_per_node; i++)	/* Граничні точки залишаються постійними*/
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
	/* Ініціалізуємо і встановлюємо початкові значення внутрішнього діапазону */
	node_mat = matrix_alloc(points_per_node, points_per_node, 0);
	/* Встановлюємо початкові значення граничних вузлів сітки дискретизації */
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		node_mat->data[i] = grid_init->data[i];
	}
	/* Створимо масив значень функції в локальних вузлах сітки дискретизації*/
	f_xy = matrix_alloc(points_per_node, points_per_node, 0.0);
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		f_xy->data[i] = func(coords_X->data[i], coords_Y->data[i]);
	}
	/*Визначимо типи для передач під час ітерацій*/
	/* Тип, що визначає стовпець матриці*/
	MPI_Type_vector(points_per_node, 1, points_per_node, MPI_DOUBLE, &mat_col);
	MPI_Type_commit(&mat_col);			/* Реєструємо тип */
	left_col = (double*)malloc(points_per_node * points_per_node * sizeof(double));
	right_col = (double*)malloc(points_per_node * points_per_node * sizeof(double));
	top_row = (double*)malloc(points_per_node * points_per_node * sizeof(double));
	bot_row = (double*)malloc(points_per_node * points_per_node * sizeof(double));
	node_iter = (int*)malloc(points_per_node * points_per_node * sizeof(int));
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		node_iter[i] = 0;
	}
	/* Початок ітерацій */
	do
	{
		/* Обмін між процесами у сітці*/
		if (node_X != 0)
		{
			/* Приймемо стовпець з лівого процесу */
			MPI_Recv(left_col, 1, mat_col, rank - 1, rank - 1, MPI_COMM_WORLD, &stat);
		}
		if (node_X != nodes_width - 1)
		{
			/* Відішлемо стовпець до правого процесу */
			MPI_Send(node_mat->data + points_per_node - 1, 1, mat_col,
				rank + 1, rank, MPI_COMM_WORLD);
			/* Приймемо стовпець з правого процесу */
			MPI_Recv(right_col, 1, mat_col, rank + 1, rank + 1, MPI_COMM_WORLD, &stat);
		}
		if (node_X != 0)
		{
			/* Відішлемо стовпець до лівого процесу */
			MPI_Send(node_mat->data, 1, mat_col, rank - 1, rank, MPI_COMM_WORLD);
		}
		if (node_Y != 0)
		{
			/* Приймемо строку з верхнього процесу */
			MPI_Recv(top_row, points_per_node, MPI_DOUBLE, rank - nodes_width,
				rank - nodes_width, MPI_COMM_WORLD, &stat);
		}
		if (node_Y != nodes_width - 1)
		{
			/* Відішлемо строку до нижнього процесу */
			MPI_Send(node_mat->data + (points_per_node - 1) * points_per_node, points_per_node
				, MPI_DOUBLE, rank + nodes_width, rank, MPI_COMM_WORLD);
			/* Приймемо строку з нижнього процесу */
			MPI_Recv(bot_row, points_per_node, MPI_DOUBLE, rank + nodes_width,
				rank + nodes_width, MPI_COMM_WORLD, &stat);
		}
		if (node_Y != 0)
		{
			/* Відішлемо строку до верхнього процесу */
			MPI_Send(node_mat->data, points_per_node, MPI_DOUBLE,
				rank - nodes_width, rank, MPI_COMM_WORLD);
		}
		/* Виклик функції, що проводить ітераційну обробку*/
		for (int i = 0; i < points_per_node; i++)
		{
			for (int j = 0; j < points_per_node; j++)
			{
				evaluate(i * points_per_node + j);
			}
		}
		/*Якщо усі вузли сітки дискретизації у процесі завершили обчислення –
		 встановлюємо ознаку завершення обчислень у процесі*/
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
		/*Перевіряємо чи всі процеси закінчили обчислення*/
		MPI_Reduce(&node_stop, &global_stop, 1, MPI_SHORT, MPI_BAND, 0, MPI_COMM_WORLD);
		/*Якщо всі - то посилаємо сигнал про завершення*/
		MPI_Bcast(&global_stop, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
	} while (!global_stop);
	/* Збираємо результат */
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
					/* Приймаємо частини матриці сітки дискретизації з кожного процесу */
					MPI_Recv(temp->data, points_per_node * points_per_node,
						MPI_DOUBLE, i * nodes_width + j, i * nodes_width + j, MPI_COMM_WORLD, &stat);
					/* Та розташовуємо їх у повній матриці */
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
/* Допоміжні функції */
double func(double X, double Y)
{
	return X * Y;
}
/*Функція, що містить ітераційні формули */
void inline evaluate(int index)
{
	/*Якщо досягнута потрібна точність, то не потрібно проводити обчислення */
	if (local_stop[index])
	{
		return;
	}
	/*Якщо вузол сітки знаходиться на границі між обчислювальними
	 вузлами, то потрібно взяти значення з буферу передачі*/
	double left = index % points_per_node == 0
		? left_col[index] : node_mat->data[index - 1];
	double right = (index + 1) % points_per_node == 0
		? right_col[index + 1 - points_per_node] : node_mat->data[index + 1];
	double top = index / points_per_node == 0 ? top_row[index % points_per_node]
		: node_mat->data[index - points_per_node];
	double bottom = index / points_per_node == points_per_node - 1
		? bot_row[index % points_per_node]
		: node_mat->data[index + points_per_node];
	/*Власне обчислення*/
	double old_node = node_mat->data[index];
	double LU = old_node - (left + right + top + bottom - 4 * old_node) / (h * h);
	node_mat->data[index] = old_node - omega * h * h / 4 * (LU - f_xy->data[index]);
	/*Якщо досягнута точність або перевищена максимальна кількість ітерацій,
	то завершуємо обчислення у цьому вузлі сітки*/
	if (fabs(node_mat->data[index] - old_node) <=
		fa bs(node_mat->data[index] * epsilon) || ++node_iter[index] > max_iter)
	{
		local_stop[index] = true;
	}
}
/* Ініціалізує значення сітки дискретизації у процесі */
void init_grid()
{
	for (int i = 0; i < points_per_node * points_per_node; i++)
	{
		node_mat->data[i] = grid_init->data[i];
	}
}
void fatal_error(const char* message, int errorcode)
{
	printf(”fatal error : code % d, % s\n”, errorcode, message);
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
	FILE* f = fopen(filename, ”w”);
	if (f == NULL)
	{
		fatal_error(”cant write to file”, 2);
	}
	for (int i = 0; i < mat->rows; i++)
	{
		for (int j = 0; j < mat->cols; j++)
		{
			fprintf(f, ” % lf ”, mat->data[i * mat->cols + j]);
		}
		fprintf(f, ”\n”);
	}
	fclose(f);
}
struct my_matrix* read_matrix(const char* filename)
{
	FILE* mat_file = fopen(filename, ”r”);
	if (mat_file == NULL)
	{
		fatal_error(”can’t open matrix file”, 1);
	}
	int rows;
	int cols;
	fscanf(mat_file, ” % d % d”, &rows, &cols);
	struct my_matrix* result = matrix_alloc(rows, cols, 0.0);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			fscanf(mat_file, ” % lf”, &result->data[i * cols + j]);
		}
	}
	fclose(mat_file);
	return result;
}
