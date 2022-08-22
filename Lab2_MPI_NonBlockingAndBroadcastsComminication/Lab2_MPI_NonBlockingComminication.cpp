#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
const int MY_TAG = 42;
/*
* �������� ����������� ���� �� 2 �������. ³� ������� 0 ��
* ������� 1 ���������� ����� � 5 ����� ����� �� ���������
* ���������� ��������. ��� ����� ������ 1 ������ ������
* ����� ����� �� ������ 0 ������ �� ��������.
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
		sleep(5);				 /* ��������. � ������� ������� -- ����������. */
		int data[] = { 1, 2, 3, 4, 5 };		 /* ��� ��� �������� */
		MPI_Request send_req; 		/* ���������� �������� ����� */
		printf(�Task 0. MPI_Isend()\n�);
		MPI_Isend(&data, 5, MPI_INT, 1, MY_TAG, MPI_COMM_WORLD, &send_req);
		/* ����� ����� data ����������� */
		printf(�Task 0. Sending data...\n�);
		sleep(5);					 /* ��������. � ������� ������� -- ����������. */
		printf(�Task 0. Computation finished.\n�);
		/* ���������� ���������� �������� */
		MPI_Wait(&send_req, MPI_STATUS_IGNORE);
	}
	if (rank == 1)
	{
		int data[5];					 /* ����� ��� ������� ����� */
		MPI_Request recv_req;			 /* ���������� ������� ����� */
		printf(�Task 1. MPI_Irecv()\n�);
		MPI_Irecv(&data, 5, MPI_INT, 0, MY_TAG, MPI_COMM_WORLD, &recv_req);
		/* ����� ����� data ����������� */
		printf(�Task 1. Recieving data...\n�);
		sleep(2);					 /* ��������. � ������� ������� -- ����������. */
		/* ���������� �� ���������� �������� */
		MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
		printf(�Task 1. Recieved data : �);
		/* ���� ����� */
		for (int i = 0; i < 5; i++)
		{
			printf(� % d�, data[i]);
		}
		printf(�\n�);
	}
	MPI_Finalize();
	return 0;
}
