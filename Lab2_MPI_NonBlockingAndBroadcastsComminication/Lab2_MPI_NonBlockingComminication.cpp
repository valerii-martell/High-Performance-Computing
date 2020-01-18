#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
const int MY_TAG = 42;
/*
* Програма розрахована рівно на 2 процеси. Від процесу 0 до
* процесу 1 передається масив з 5 цілих чисел за допомогою
* неблокуючої передачі. При цьому процес 1 починає прийом
* даних раніше ніж процес 0 починає їх передачу.
*/
int main(int argc, char* argv[])
{
	int np;
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		sleep(5);				 /* Затримка. В реальній програмі -- обчислення. */
		int data[] = { 1, 2, 3, 4, 5 };		 /* Дані для передачі */
		MPI_Request send_req; 		/* Дескриптор передачі даних */
		printf(”Task 0. MPI_Isend()\n”);
		MPI_Isend(&data, 5, MPI_INT, 1, MY_TAG, MPI_COMM_WORLD, &send_req);
		/* Тепер буфер data недоступний */
		printf(”Task 0. Sending data...\n”);
		sleep(5);					 /* Затримка. В реальній програмі -- обчислення. */
		printf(”Task 0. Computation finished.\n”);
		/* Очікування завершення передачі */
		MPI_Wait(&send_req, MPI_STATUS_IGNORE);
	}
	if (rank == 1)
	{
		int data[5];					 /* Буфер для прийому даних */
		MPI_Request recv_req;			 /* Дескриптор прийому даних */
		printf(”Task 1. MPI_Irecv()\n”);
		MPI_Irecv(&data, 5, MPI_INT, 0, MY_TAG, MPI_COMM_WORLD, &recv_req);
		/* Тепер буфер data недоступний */
		printf(”Task 1. Recieving data...\n”);
		sleep(2);					 /* Затримка. В реальній програмі -- обчислення. */
		/* Очікування на завершення передачі */
		MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
		printf(”Task 1. Recieved data : ”);
		/* Вивід даних */
		for (int i = 0; i < 5; i++)
		{
			printf(” % d”, data[i]);
		}
		printf(”\n”);
	}
	MPI_Finalize();
	return 0;
}
