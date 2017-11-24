#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
//#define N 256//矩阵的维度

void creat_matrix(int *,int);//创建一个随机矩阵
void normal_matrix_mult(int *,int *,int *,int);//普通的暴力矩阵乘法 O(n^3)
void strassen_matrix_mult(int *,int *, int *,int);//strassen算法 O(n^lg7)
void print_matrix(int *,int);//输出矩阵
void matrix_Addition(int *, int *, int *,int);//矩阵加法
void matrix_Subtraction(int *, int *, int *,int);//矩阵减法

void creat_matrix(int *matrix, int dim){//创建一个随机矩阵
    int i;
    int temp;
    for(i=0;i<dim*dim;i++){
        temp=rand()%100;
        matrix[i]=temp;
    }
}

void normal_matrix_mult(int *matrix_A, int *matrix_B, int *ans, int dim){//普通的暴力矩阵乘法 O(n^3)
    int i, j, k, temp;

    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            temp = 0;
            for(k = 0; k < dim; k++){
                temp += matrix_A[i*dim+k] * matrix_B[k*dim+j];
            }
            ans[i*dim+j] = temp;
        }
    }
}



void print_matrix(int *matrix,int dim){//输出矩阵
    int i, j;
    for (i = 0; i <dim; i++){
        for (j = 0; j <dim; j++)
            printf("%d ",matrix[i*dim+j]);
            printf("\n");
    }
    printf("\n");
}


void matrix_Addition(int *matrix_A, int *matrix_B, int *ans,int dim){//矩阵加法
    int i, j;
    for (i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            ans[i*dim+j] = matrix_A[i*dim+j]+matrix_B[i*dim+j];
        }
    }
}


void matrix_Subtraction(int *matrix_A, int *matrix_B, int *ans,int dim){//矩阵减法
    int i, j;
    for (i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            ans[i*dim+j] = matrix_A[i*dim+j]-matrix_B[i*dim+j];
        }
    }
}


void strassen_matrix_mult(int *A, int *B, int *C,int dim){//strassen算法 O(n^lg7)
    int i, j;
    int new_i;
    int new_j;
    int n = dim/2;

    if(dim<=32){//理论上应该按照下面算法导论的写法，递归到dim==1，但是实际上小矩阵的时候用普通暴力算法快得多，否则小矩阵都要递归的话，递归
            //到最后会做大量无用功，经验告诉我们这个小矩阵的分界点大概在32，这种操作不会改变时间复杂度
     normal_matrix_mult(A, B, C,dim);
        return;
    }

/*
    if (dim==1){
        C[0]=A[0]*B[0];
        return;
    }
*/

    int *A11 = malloc(n*n*sizeof(int));//矩阵A分块4份
    int *A12 = malloc(n*n*sizeof(int));
    int *A21 = malloc(n*n*sizeof(int));
    int *A22 = malloc(n*n*sizeof(int));

    int *B11 = malloc(n*n*sizeof(int));//矩阵B分块4份
    int *B12 = malloc(n*n*sizeof(int));
    int *B21 = malloc(n*n*sizeof(int));
    int *B22 = malloc(n*n*sizeof(int));

    int *P1 = malloc(n*n*sizeof(int));
    int *P2 = malloc(n*n*sizeof(int));
    int *P3 = malloc(n*n*sizeof(int));
    int *P4 = malloc(n*n*sizeof(int));
    int *P5 = malloc(n*n*sizeof(int));
    int *P6 = malloc(n*n*sizeof(int));
    int *P7 = malloc(n*n*sizeof(int));

    int *C11 = malloc(n*n*sizeof(int));//矩阵C分块4份
    int *C12 = malloc(n*n*sizeof(int));
    int *C21 = malloc(n*n*sizeof(int));
    int *C22 = malloc(n*n*sizeof(int));

    int *S1 = malloc(n*n*sizeof(int));//中间量
    int *S2 = malloc(n*n*sizeof(int));

//给分块矩阵赋值
    for (new_i = 0, i = 0; new_i < n; new_i++, i++){
        for (new_j = 0, j = 0; new_j < n; new_j++, j++){
            A11[i*n+j] = A[new_i*dim+new_j];
            B11[i*n+j] = B[new_i*dim+new_j];
        }
        for (new_j= n, j = 0; new_j < dim; new_j++, j++){
            A12[i*n+j] = A[new_i*dim+new_j];
            B12[i*n+j] = B[new_i*dim+new_j];
        }
    }

    for (new_i = n, i = 0; new_i < dim; new_i++, i++){
        for (new_j = 0, j = 0; new_j < n; new_j++, j++){
            A21[i*n+j] = A[new_i*dim+new_j];
            B21[i*n+j] = B[new_i*dim+new_j];
        }
        for (new_j= n, j = 0; new_j< dim; new_j++, j++){
            A22[i*n+j] = A[new_i*dim+new_j];
            B22[i*n+j] = B[new_i*dim+new_j];
        }
    }

    //P1=(A11+A22)(B11+B22)
    matrix_Addition(A11, A22, S1, n);
    matrix_Addition(B11, B22, S2, n);
    strassen_matrix_mult(S1, S2, P1, n);

    //P2=(A21+A22)B11
    matrix_Addition(A21, A22, S1, n);
    strassen_matrix_mult(S1, B11, P2, n);

    //P3=A11(B12-B22)
    matrix_Subtraction(B12, B22, S2, n);
    strassen_matrix_mult(A11, S2, P3, n);

    //P4=A22(B21-B11)
    matrix_Subtraction(B21, B11, S2, n);
    strassen_matrix_mult(A22, S2, P4, n);

    //P5=(A11+A12)B22
    matrix_Addition(A11, A12, S1, n);
    strassen_matrix_mult(S1, B22, P5, n);

    //P6=(A21-A11)(B11+B12)
    matrix_Subtraction(A21, A11, S1, n);
    matrix_Addition(B11, B12, S2, n);
    strassen_matrix_mult(S1, S2, P6, n);

    //P7=(A12-A22)(B21+B22)
    matrix_Subtraction(A12, A22, S1, n);
    matrix_Addition(B21, B22, S2, n);
    strassen_matrix_mult(S1, S2, P7, n);

    //C11=(P1+P4-P5)
    matrix_Addition(P1, P4, C11, n);
    matrix_Subtraction(C11, P5, C11, n);
    matrix_Addition(C11, P7, C11, n);

    //C12=P3+P5
    matrix_Addition(P3, P5, C12, n);

    //C21=P2+P4
    matrix_Addition(P2, P4, C21, n);

    //C22=P1-P2+P3+P6
    matrix_Subtraction(P1, P2, C22, n);
    matrix_Addition(C22, P3, C22, n);
    matrix_Addition(C22, P6, C22, n);

//分块矩阵更新原矩阵
    for (new_i = 0, i = 0; new_i< n; new_i++, i++){
        for (new_j = 0, j = 0; new_j< n; new_j++, j++)
            C[new_i*dim+new_j] = C11[i*n+j];
        for (new_j = n, j = 0; new_j < dim; new_j++, j++)
            C[new_i*dim+new_j] = C12[i*n+j];
    }
    for (new_i = n, i = 0; new_i < dim; new_i++, i++){
        for (new_j = 0, j = 0; new_j < n;new_j++, j++)
            C[new_i*dim+new_j] = C21[i*n+j];
        for (new_j = n, j = 0;new_j < dim; new_j++, j++)
            C[new_i*dim+new_j] = C22[i*n+j];
    }

    free(S1);
    free(S2);
    free(A11);
    free(A12);
    free(A21);
    free(A22);
    free(B11);
    free(B12);
    free(B21);
    free(B22);
    free(P1);
    free(P2);
    free(P3);
    free(P4);
    free(P5);
    free(P6);
    free(P7);
    free(C11);
    free(C12);
    free(C21);
    free(C22);


}



int main()
{
    int *matrix_A;
    int *matrix_B;
    int *ans;
    int m=10;
    int N=16;
    int p;
    int temp=N;
    double *normal_ans;
    double *strassen_ans;
    normal_ans=malloc(m*sizeof(double));
    strassen_ans=malloc(m*sizeof(double));

    for(p=4;p<m;p++){
    N*=2;
    clock_t timer;//计时器
    srand(time(NULL));//为生成随机矩阵播种

    matrix_A=malloc(N*N*sizeof(int));
    matrix_B=malloc(N*N*sizeof(int));
    ans=malloc(N*N*sizeof(int));

    creat_matrix(matrix_A,N);
    creat_matrix(matrix_B,N);//初始化两个待乘矩阵

    /*
    对于普通算法的计时：
    */
    timer=clock();
    normal_matrix_mult(matrix_A,matrix_B,ans,N);
    timer=clock()-timer;
    normal_ans[p]=((float)timer)/CLOCKS_PER_SEC;
   // printf("%f\n", ((float)timer)/CLOCKS_PER_SEC);
   // print_matrix(ans);

    /*
    对于strassen算法的计时：
    */
    timer=clock();
    strassen_matrix_mult(matrix_A,matrix_B,ans,N);
    timer=clock()-timer;
    strassen_ans[p]=((float)timer)/CLOCKS_PER_SEC;
   // printf("%f\n", ((float)timer)/CLOCKS_PER_SEC);
  //  print_matrix(ans);

    free(matrix_A);
    free(matrix_B);
    free(ans);
    }
    for(p=4;p<m;p++){
        printf("%lf ",normal_ans[p]);
    }
    printf("\n");

    for(p=4;p<m;p++){
        printf("%lf ",strassen_ans[p]);
    }
        printf("\n");
    for(p=4;p<m;p++){
            temp*=2;
        printf("%d ",temp);
    }
    return 0;
}
