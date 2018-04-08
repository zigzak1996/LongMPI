#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "a_div_b.h"


int main(int argc, char* argv[]){
	int *A = NULL;
	int *B = NULL;
	int n,m;
	int comm_sz, my_rank;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	if (my_rank == 0){
		FILE *f;
		f = fopen("input.txt","r");
		
		char a[200000];
		char b[200000];
			

		fscanf(f, "%s", a);
		fscanf(f, "%s", b);
			
		A = malloc(sizeof(int) * (strlen(a)/3+5) );
		B = malloc(sizeof(int) * (strlen(b)/3+5) );
		
		memset(A,0,strlen(a)/3+5);
		memset(B,0,strlen(b)/3+5);	
		
		n = 0;
		m = 0;
		
		for (int i = strlen(a) - 1; i >= 0; i--, n++){
			A[n / 3] += (a[i] - '0') * power(10, n % 3);
		}
		
		if (n % 3 == 0) n /= 3;
		else n = n / 3 + 1;
		
		for (int i = strlen(b) - 1; i >= 0; i--, m++){
			B[m / 3] += (b[i] - '0') * power(10, m % 3);
		}
			
		if (m % 3 == 0) m /= 3;
		else m = m / 3 + 1;
	}

	division(A, B, n, m, my_rank, comm_sz); 
	subtract(A, B, n, m, my_rank, comm_sz);
	addition(A, B, n, m, my_rank, comm_sz);
	multiplication(A, B, n, m, my_rank, comm_sz);
	
	MPI_Finalize();
	
	return 0;
}