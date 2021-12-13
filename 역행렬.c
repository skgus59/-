#include "20211582.h"

int main(void){
	FILE *fp;//파일입출력
	int N;
	fp=fopen("matrix.txt","r");//파일오류점검
	if (fp==NULL){
		printf("input file open error!\n");
		return 1;
	}
	fscanf(fp, "%d", &N);//N입력받기
	double **matrix, **cofactormatrix, ** inversematrix;//동적할당
	matrix=(double**)malloc(N*sizeof(double*));
	for (int i=0;i<N;i++) *(matrix+i)=(double*)malloc(N*sizeof(double));
	cofactormatrix=(double**)malloc(N*sizeof(double*));
	for (int i=0;i<N;i++) *(cofactormatrix+i)=(double*)malloc(N*sizeof(double));
	inversematrix=(double**)malloc(N*sizeof(double*));
	for (int i=0;i<N;i++) *(inversematrix+i)=(double*)malloc(N*sizeof(double));

	if (N<=0) exit(1);//N이 0이하일때 제외
	else{
		for (int i=0;i<N;i++){
			for (int j=0;j<N;j++)
				fscanf(fp,"%lf", *(matrix+i)+j);//N만큼 NxN 배열 입력받기
		}
	}
	printMatrix(matrix,N);//nxn matrix print
	printf("\n");
	double A;
	A=detA(matrix,N);//double det값 detA 함수를 통해 구하기
	if (A==0.000000) exit(1);//detA가 0이면 종료
	else{
		cofactormatrix = cofactorMatrix(matrix,N);//cofactormatrix구하기
	}
	printMatrix(cofactormatrix,N);
	printf("\n");
	inversematrix = inverseMatrix(cofactormatrix,A,N);//inversematrix구하기
	printMatrix(inversematrix,N);
	for (int i=0;i<N;i++){
		free(*(inversematrix+i));//inversematrix free
	}
	free(inversematrix);
	fclose(fp);//파일 닫기
	return 0;
}

double detA(double **A, int n){
	int M=n-1;
	int K = 0;
	int plinus=1;
	double result=0;
	double **Matrix;
	Matrix=(double**)malloc(sizeof(double*)*M); //C1j를 받아올 배열
	for (int i=0; i<M; i++){
		*(Matrix+i)=(double*)malloc(sizeof(double)*M);
	}
	if (M==1){//n이 2
		result=(*(*(A+0)+0)*(*(*(A+1)+1))-*(*(A+0)+1)*(*(*(A+1)+0)));//n=2일때는 따로 계산
	}
	else if (M==0){//n이 1일때는 자기자신
		result=*(*(A+0)+0);
	}
	else{
		for (int j=0;j<n;j++){ //1j를 (1부터 n까지)
			for (int p = 0; p < M; p++) {
				K = 0;
				for (int q = 0; q < M; q++) {//세로1 가로 j에 위치한 값들은 제외하고 Matrix에 넣어야하기 때문에
					 if (q==j) K++;//q랑 j랑 같을 때는 K값을 추가해줌
					 *(*(Matrix + p) + q) = *(*(A + p + 1) + K);// 1j줄은 제외하고 나머지 MxM 배열을 Matrix에 담기
					 K++;
				}
			}
			if (j % 2 == 1) plinus = -1;//홀수인경우 -1
			else plinus = 1;//짝수인경우 +1
			result += plinus * (*(*(A + 0) + j)) * detA(Matrix, M);//a1j와 C1j를 곱해서result에 담기
		}
	}
	return result;//결과값 return
	for (int j=0;j<M;j++){
		free(*(Matrix+j));
		*(Matrix+j)=NULL;//Matrix free
	}
	free(Matrix);
	Matrix=NULL;
}

void printMatrix(double **A, int n){
	for (int p=0;p<n;p++){  //세로 p
		for (int q=0;q<n;q++)//가로q
			printf("%7.2f",*(*(A+p)+q));//소수점 세자리 출력
		printf("\n");
	}
}

double** cofactorMatrix(double** A, int n) {
	int M=n-1;
	int K, T=0;
	int plinus=1;
	double** cofactor;//cofactor 동적할당
	cofactor = (double**)malloc(sizeof(double*) * n);
	double** mmatrix;//mmatrix 동적할당
	mmatrix= (double**)malloc(sizeof(double*)*M);
	for (int x=0;x<n;x++) *(mmatrix +x)=(double*)malloc(sizeof(double)*n);
	for (int x = 0; x < n; x++) *(cofactor + x) = (double*)malloc(sizeof(double) * n);
	if (n == 2) {//n이 2일때a는따로 계산
		*(*(cofactor + 0) + 0) = *(*(A + 1) + 1); //cofactor 구하기
		*(*(cofactor + 0) + 1) = *(*(A + 1) + 0);
		*(*(cofactor + 1) + 0) = *(*(A + 0) + 1);
		*(*(cofactor + 1) + 1) = *(*(A + 0) + 0);
	}
	else if (n == 1){//n이 1일때 따로 계산
		*(*(cofactor + 0) + 0) = *(*(A + 0) + 0);
	}
	else {
		for (int i = 0; i < n; i++) { //C1? ~Cn?까지, 세로
			for (int j = 0; j < n; j++) {//C?1~C?n까지, 가로
				K = 0;//collum
				for (int p = 0; p <M; p++) {//mmatrix의 세로
					if (K == i) K++;
					T = 0;
					for (int q = 0; q <M; q++) {//mmatrix의 가로
						if (q == j) T++;
						*(*(mmatrix+p)+q) = *(*(A+K)+T);
						T++;
					}
					K++;
				}
				if ((i+j)%2==1) plinus=-1;//홀수일때는 -1
				else plinus=1;//짝수일때는 +1
				*(*(cofactor+i)+j)=plinus*detA(mmatrix, M);
			}
		}
	}
	for(int p=0;p<n;p++){
		for (int q=0;q<n;q++){
			*(*(A+p)+q)=*(*(cofactor+p)+q);//cofactor의 값을 A에 다시 넣기
		}
	}
	return A;//A return
	for (int x=0;x<n;x++){
		free(*(mmatrix+x));//mmatrix free
		*(mmatrix+x)=NULL;
	}
	free(mmatrix);
	mmatrix=NULL;
	for (int x = 0; x < n; x++){//cofactor free
		free(*(cofactor + x));
		*(cofactor+x)=NULL;
	}
	free(cofactor);
	cofactor=NULL;
}
double** inverseMatrix(double** A, double det, int n) {
	double** inverse;
	inverse = (double**)malloc(sizeof(double*) * n);//동적할당
	for (int x = 0; x < n; x++) *(inverse + x) = (double*)malloc(sizeof(double) * n);
	for (int p = 0; p < n; p++) {//세로
		for (int q = 0; q < n; q++){//가로
			*(*(inverse + p) + q) = (*(*(A + q) + p))/det;//Cij를 Cji로바꾸기 그리고 det값으로 나누기
		}
	}
	for (int p=0;p<n;p++){//세로
		for (int q=0;q<n;q++){//가로
			*(*(A+p)+q)=*(*(inverse+p)+q);//inverse값을 a에 다시 복사
		}
	}
	return A;
	for (int x = 0; x < n; x++){//inverse free
		free(*(inverse + x));
		*(inverse + x) = NULL;
	}
	free(inverse);
	inverse = NULL;
}