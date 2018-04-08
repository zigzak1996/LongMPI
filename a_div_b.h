#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

const int base = 1000;

int compare(int *A, int *B, int len_A, int len_B){
	while(A[len_A-1] == 0) len_A--;
	int t = len_A - len_B;
	
	if (t) return t;

	for (int i = len_A - 1; i >= 0; i--){
		t = A[i] - B[i];
		if (t) return t;
	}
	return 0;
}

int power(int x,int y){
	int res = 1;
	for (int i = 1; i <= y; i++) res *= x;
	return res;
}

int* multiply(int *A, int *sz, int q, int comm_sz, int my_rank){
	int *partA = NULL;
	partA = malloc(sizeof(int)*(*sz / comm_sz+2));
	MPI_Scatter(A, *sz/comm_sz, MPI_INT, partA, *sz/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	int *C = NULL;
	C = malloc(sizeof(int) * (*sz / comm_sz + 2));
	memset(C, 0, *sz / comm_sz + 2);
	for (int i = 0; i < *sz / comm_sz + 1; i++) C[i] = 0;
	
	for (int i = 0; i < *sz / comm_sz; i++){
		C[i] += (partA[i] * q) % base;
		C[i + 1] += (partA[i] * q) / base + C[i] / base;
		C[i] %= base;
	}

	if (my_rank > 0 && my_rank < comm_sz - 1){
		int k;
		MPI_Recv(&k, 1, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int i = 0;
		C[i] += k;
		while(C[i] > 999){
			C[i] %= base;
			C[i + 1] += 1;
			i++; 
		}
		MPI_Send(C + (*sz / comm_sz), 1, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank == 0 && comm_sz > 1){
		MPI_Send(C + (*sz / comm_sz), 1, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank == comm_sz - 1 && comm_sz > 1){
		int k;
		MPI_Recv(&k, 1, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int i = 0;
		C[i] += k;
		while(C[i] > 999){
			C[i] %= base;
			C[i + 1] += 1;
			i++; 
		}
	}
	int *ans = NULL;
	ans = malloc(sizeof(int) * (*sz + 1));
	MPI_Gather(C, *sz / comm_sz, MPI_INT, ans, *sz / comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	if(my_rank == comm_sz -1){
		MPI_Send(C + (*sz / comm_sz), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	if(my_rank == 0){
		MPI_Recv(ans+(*sz), 1, MPI_INT, comm_sz - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if( ans[ *sz ] > 0) *sz = *sz + 1;
		else ans = realloc(ans, sizeof(int) * ( *sz ));
		
		
		return ans;
	}
	
	return C;
}

int* minus(int *A,int *B, int n, int m, int comm_sz, int my_rank){
	int *a = NULL;
	int *b = NULL;

	a = malloc(sizeof(int) * n);
	b = malloc(sizeof(int) * m);

	for (int i = 0; i < n; i++) a[i] = A[i];
	for (int i = 0; i < m; i++) b[i] = B[i];	
	
	int N = n, M = m;

	while(N % comm_sz > 0)N++;
	while(M < N)M++;

	a = realloc(a, sizeof(int) * N);
	b = realloc(b, sizeof(int) * M);
	
	memset(a+n, N-n, 0);
	memset(b+m, M-m, 0);
	
	int *partA = NULL;
	int *partB = NULL;

	partA = malloc(sizeof(int) * (N/comm_sz));
	partB = malloc(sizeof(int) * (M/comm_sz));
	
	MPI_Scatter(a, N/comm_sz, MPI_INT, partA, N/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(b, M/comm_sz, MPI_INT, partB, M/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	
	for (int i = 0; i < N/comm_sz; i++){
		partA[i] -= partB[i];
		if (partA[i] < 0){
			partA[i] += base;
			partA[i+1]--;
		}
	}
	if(my_rank+1<comm_sz){
		MPI_Send(partA + N / comm_sz, 1, MPI_INT, my_rank + 1,0, MPI_COMM_WORLD);
	}
	if(my_rank-1>=0){
		int ost;
		MPI_Recv(&ost, 1, MPI_INT, my_rank - 1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		partA[0] -= ost;
		int i = 0;
		while(partA[i] < 0){
			partA[i + 1]--;
			partA[i] += base;
			i++;
		}
	}
	MPI_Gather(partA, N/comm_sz, MPI_INT, a, N/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(partB, M/comm_sz, MPI_INT, b, M/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
		
	for (int i = 0; i < n; i++) A[i] = a[i];
}
void division(int* A, int* B, int n, int m, int my_rank, int comm_sz){
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (my_rank != 0){
		A = malloc(sizeof(int) * (n+5));
		B = malloc(sizeof(int) * (m+5));
		
	}

	int *a = NULL;
	int *b = NULL;

	a = malloc(sizeof(int) * (n+5));
	b = malloc(sizeof(int) * (m+5));

	for (int i = 0; i < n; i++) a[i] = A[i];
	for (int i = 0; i < m; i++) b[i] = B[i];	
	
	int N = m + 1;

	int *C = NULL;
	int *DIV = NULL;
	C = malloc(sizeof(int) * (n - m));
	if (my_rank == 0)a[++n] = 0;
	else n++;
	for (int i = n - m; i >= 0; i--){
		int ql = 0, qr = base, q;
		int t;
		for (ql, qr; qr - ql > 1;){
			q = (ql + qr) / 2;
			t = m;
			DIV = multiply(b, &t, q, comm_sz, my_rank);
			if(my_rank == 0){
				if ( compare(a + i, DIV, N, t) < 0){
					qr = q;
				} else {
					ql = q;
				}
			}
			MPI_Bcast(&ql, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&qr, 1, MPI_INT, 0, MPI_COMM_WORLD);
				
		}
		q = ql;
		t = m;
		if(my_rank == 0)C[i] = q;
		if (q == 0) continue;
		
		DIV = multiply(b, &t, q, comm_sz, my_rank);
			
		minus(a + i, DIV, N, t, comm_sz, my_rank);
	}
	if(my_rank == 0){
		for (int i = n - m; i>=0; i--) printf("%d ", C[i]);
		printf("\n");
	}
}

void subtract(int* A, int* B, int n, int m, int my_rank, int comm_sz){
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (my_rank != 0){
		A = malloc(sizeof(int) * (n+5));
		B = malloc(sizeof(int) * (m+5));
		
	}
	int *a = NULL;
	int *b = NULL;
	
	a = malloc(sizeof(int) * n);
	b = malloc(sizeof(int) * m);

	for (int i = 0; i < n; i++) a[i] = A[i];
	for (int i = 0; i < m; i++) b[i] = B[i];	
	
	int N = n, M = m;
			
	while(N % comm_sz > 0 && N < M)N++;
	if (M < N) M = N;

	a = realloc(a, sizeof(int) * N);
	b = realloc(b, sizeof(int) * M);
	
	memset(a+n, N-n, 0);
	memset(b+m, M-m, 0);
	
	int *partA = NULL;
	int *partB = NULL;

	partA = malloc(sizeof(int) * (N/comm_sz));
	partB = malloc(sizeof(int) * (M/comm_sz));
	
	MPI_Scatter(a, N/comm_sz, MPI_INT, partA, N/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(b, M/comm_sz, MPI_INT, partB, M/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < N/comm_sz; i++){
		partA[i] -= partB[i];
		if (partA[i] < 0){
			partA[i] += base;
			partA[i+1]--;
		}
	}
	if(my_rank+1<comm_sz){
		MPI_Send(partA + N / comm_sz, 1, MPI_INT, my_rank + 1,0, MPI_COMM_WORLD);
	}
	if(my_rank-1>=0){
		int ost;
		MPI_Recv(&ost, 1, MPI_INT, my_rank - 1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		partA[0] -= ost;
		int i = 0;
		while(partA[i] < 0){
			partA[i + 1]--;
			partA[i] += base;
			i++;
		}
	}
	MPI_Gather(partA, N/comm_sz, MPI_INT, a, N/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(partB, M/comm_sz, MPI_INT, b, M/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
		
	if(my_rank == 0){
		for (int i = N; i >= 0; i--) printf("%d ", a[i]);
		printf("\n");
	}
}

void addition(int* A, int* B, int n, int m, int my_rank, int comm_sz){
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (my_rank != 0){
		A = malloc(sizeof(int) * (n+5));
		B = malloc(sizeof(int) * (m+5));
		
	}
	int *a = NULL;
	int *b = NULL;

	a = malloc(sizeof(int) * n);
	b = malloc(sizeof(int) * m);

	for (int i = 0; i < n; i++) a[i] = A[i];
	for (int i = 0; i < m; i++) b[i] = B[i];	
	
	int N = n, M = m;

	while(N % comm_sz > 0 && N < M)N++;
	if (M < N) M = N;

	a = realloc(a, sizeof(int) * N);
	b = realloc(b, sizeof(int) * M);
	
	memset(a+n, N-n, 0);
	memset(b+m, M-m, 0);
	
	int *partA = NULL;
	int *partB = NULL;

	partA = malloc(sizeof(int) * (N/comm_sz + 10));
	partB = malloc(sizeof(int) * (M/comm_sz + 10));
	
	MPI_Scatter(a, N/comm_sz, MPI_INT, partA, N/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(b, M/comm_sz, MPI_INT, partB, M/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	
	for (int i = 0; i < N/comm_sz; i++){
		partA[i] += partB[i];
		if (partA[i] > base - 1){
			partA[i] -= base;
			partA[i + 1]++;
		}
	}
	if(my_rank+1<comm_sz){
		MPI_Send(partA + N / comm_sz, 1, MPI_INT, my_rank + 1,0, MPI_COMM_WORLD);
	}
	if(my_rank-1>=0){
		int ost;
		MPI_Recv(&ost, 1, MPI_INT, my_rank-1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		partA[0]+=ost;
		int i=0;
		while(partA[i] > base - 1){
			partA[i+1]++;
			partA[i] %= base;
			i++;
		}
	}
	MPI_Gather(partA, N/comm_sz, MPI_INT, a, N/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(partB, M/comm_sz, MPI_INT, b, M/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
		
	if(my_rank == 0){
		for (int i = N; i >= 0; i--) printf("%d ", a[i]);
		printf("\n");
	}
}

void multiplication(int* A, int* B, int n, int m, int my_rank, int comm_sz){
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (my_rank != 0){
		A = malloc(sizeof(int) * (n+5));
		B = malloc(sizeof(int) * (m+5));
		
	}
	int *a = NULL;
	int *b = NULL;
	int *ans = NULL;
	int *c = NULL;
	int *partA =NULL;

	a = malloc(sizeof(int) * n);
	b = malloc(sizeof(int) * m);

	for (int i = 0; i < n; i++) a[i] = A[i];
	for (int i = 0; i < m; i++) b[i] = B[i];	
	

	MPI_Bcast(b, m, MPI_INT, 0, MPI_COMM_WORLD);
	
	partA = malloc((n/comm_sz+1)*sizeof(int));

	MPI_Scatter(a, n/comm_sz, MPI_INT, partA, n/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	
	c = malloc(sizeof(int) * (n/comm_sz + m + 1));

	for (int i = 0; i < n/comm_sz + m; i++ )c[i] = 0;
	
	for (int i = 0; i < n/comm_sz; i++){
		for (int j = 0; j < m; j++){
			c[i + j] += (partA[i] * b[j]) % base;
			c[i + j + 1] += c[i+j] / base + (partA[i] * b[j]) / base;
			c[i + j] %= base;			
		}
	}
	if ( my_rank == 0 ){
		ans = malloc(sizeof(int) * (n + m + 1));
		
		for(int i = 0; i < n + m; i++) ans[i] = 0;
		
		MPI_Gather(c, n/comm_sz, MPI_INT, ans, n/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
		
		int *remainder = NULL;
		remainder = malloc(sizeof(int) * (m+1));
		int count = 0;		

		for (int i = n/comm_sz; i < n/comm_sz + m; i++){
			ans[i] += c[i];
			ans[i + 1] = ans[i + 1] + ans[i] / base;
			ans[i] = ans[i] % base;
		}
		
		for (int i = 1; i < comm_sz; i++){
			MPI_Recv(remainder, m, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int j = n/comm_sz * (i + 1); j < n/comm_sz * (i + 1) + m; j++){
				ans[j] += remainder[j - n/comm_sz * (i + 1)];
				ans[j + 1] = ans[j + 1] + ans[j] / base;
				ans[j] = ans[j] % base;		
			}	
		}
		
		for (int i = n + m - 1; i>=0; i--){
			if (ans[i] < 10)printf("00%d", ans[i]);
			else if(ans[i] < 100)printf("0%d", ans[i]);
			else printf("%d", ans[i]);
		}
		
		free(ans);
		free(remainder);
	}
	else {
		MPI_Gather(c, n/comm_sz, MPI_INT, ans, n/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Send(c + n/comm_sz, m - 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}
/*
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
*/