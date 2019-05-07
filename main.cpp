#include <stdio.h>
#include <stdlib.h>
#include <mm_malloc.h>
#include <time.h>

int pow(int a,int b){
    long int p=1;
    if (b==1) {
        return a;
    }
    while (b!=0){
        p*=a;
        b-=1;
    }
    return p;
}

int ** matrix_add(int n,int **p,int **q){
    int i,j;
    int **a;
    a=(int**)malloc(sizeof(int*)*n);
    for (i=0;i<n;i++)
    {
        *(a+i) = (int*)malloc(sizeof(int)*n);
    }
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            a[i][j]=p[i][j]+q[i][j];
        }
    }
    return a;
}

int ** matrix_sub(int n,int **p,int **q){
    int i,j;
    int **a;
    a=(int**)malloc(sizeof(int*)*n);
    for (i=0;i<n;i++)
    {
        *(a+i) = (int*)malloc(sizeof(int)*n);
    }
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            a[i][j]=p[i][j]-q[i][j];
        }
    }
    return a;
}

int ** matrix_mul(int n,int **p,int **q) {
    int i, j, k, sum;
    int **a = (int **)malloc(sizeof(int*)* n);
    for (i = 0; i < n; i++)
    {
        *(a + i) = (int *)malloc(sizeof(int) * n);
    }
    for (k = 0; k < n; k++) {
        for (i = 0; i < n; i++) {
            sum=0;
            for (j = 0; j < n; j++) {
                sum += p[k][j] * q[j][i];
            }
            a[k][i]=sum;
        }
    }
    return a;
}


int **native(int l,int **p, int **q) {
    int **c;
    int i,j,k,sum;

    c=(int**)malloc(sizeof(int*)*l);
    for (i=0;i<l;i++){
        *(c+i)=(int*)malloc(sizeof(int)*l);
    }
    for (k=0;k<l;k++) {
        for (i=0;i<l;i++) {
            sum = 0;
            for (j=0;j<l;j++) {
                sum = sum + p[i][j] * q[j][k];
            }
            c[i][k] =sum;
        }
    }

    return c;
}

int **strassen(int n,int **p,int **q) {
    int i, j;
    int **c;
    int **a11, **a12, **a21, **a22, **b11, **b12, **b21, **b22;
    int **a11_a22,**a21_a22,**a11_a12,**a21_a11,**a12_a22,**b11_b22,**b12_b22,**b21_b11,**b11_b12,**b21_b22;
    int **m1,**m2,**m3,**m4,**m5,**m6,**m7;
    int **c11,**c12,**c21,**c22,**c111,**c112,**c221,**c222;


    if (n == 2) {
        c = matrix_mul(n,p, q);
        return c;
    } else {
        int m=n/2;
        a11 = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
            *(a11 + i) = (int *) malloc(sizeof(int) * m);
        }
        a12 = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
            *(a12 + i) = (int *) malloc(sizeof(int) * m);
        }
        a21 = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
            *(a21 + i) = (int *) malloc(sizeof(int) * m);
        }
        a22 = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
            *(a22 + i) = (int *) malloc(sizeof(int) * m);
        }
        b11 = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
            *(b11 + i) = (int *) malloc(sizeof(int) * m);
        }
        b12 = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
            *(b12 + i) = (int *) malloc(sizeof(int) * m);
        }
        b21 = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
            *(b21 + i) = (int *) malloc(sizeof(int) * m);
        }
        b22 = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
            *(b22 + i) = (int *) malloc(sizeof(int) * m);
        }


        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                a11[i][j] = p[i][j];
                a12[i][j] = p[i][j + m];
                a21[i][j] = p[i + m][j];
                a22[i][j] = p[i + m][j + m];
                b11[i][j] = q[i][j];
                b12[i][j] = q[i][j + m];
                b21[i][j] = q[i + m][j];
                b22[i][j] = q[i + m][j + m];

            }
        }

        a11_a22=matrix_add(m,a11,a22);
        a21_a22=matrix_add(m,a21,a22);
        a11_a12=matrix_add(m,a11,a12);
        a21_a11=matrix_sub(m,a21,a11);
        a12_a22=matrix_sub(m,a12,a22);
        b11_b22=matrix_add(m,b11,b22);
        b12_b22=matrix_sub(m,b12,b22);
        b21_b11=matrix_sub(m,b21,b11);
        b11_b12=matrix_add(m,b11,b12);
        b21_b22=matrix_add(m,b21,b22);

        m1=strassen(m,a11_a22,b11_b22);
        m2=strassen(m,a21_a22,b11);
        m3=strassen(m,a11,b12_b22);
        m4=strassen(m,a22,b21_b11);
        m5=strassen(m,a11_a12,b22);
        m6=strassen(m,a21_a11,b11_b12);
        m7=strassen(m,a12_a22,b21_b22);

        c111=matrix_add(m,m1,m4);
        c112=matrix_sub(m,m7,m5);
        c221=matrix_sub(m,m1,m2);
        c222=matrix_add(m,m3,m6);
        c11=matrix_add(m,c111,c112);
        c12=matrix_add(m,m3,m5);
        c21=matrix_add(m,m2,m4);
        c22=matrix_add(m,c221,c222);

        int **c = (int **)malloc(sizeof(int*) * n);
        for (i = 0; i < n; i++)
        {
            *(c + i) = (int *)malloc(sizeof(int)*n);
        }

        for (i=0;i<m;i++){
            for(j=0;j<m;j++){
                c[i][j]=c11[i][j];
                c[i][j+m]=c12[i][j];
                c[i+m][j]=c21[i][j];
                c[i+m][j+m]=c22[i][j];
            }
        }

        return c;
    }
}


int main() {
    int l,i,j,n,a;
    int **p,**q,**c;
    int nat_start,nat_end,nat_running_time,str_start,str_end,str_running_time;
    FILE *fp=fopen("input.txt","r+");
    if (fp==NULL)
    {
        printf("No File Found!\n");
        return -1;
    }
    fscanf(fp,"%d",&l);

    for (a=0;a<9;a++){
        if (l<pow(2,a)){
            break;
        }
    }
    n=pow(2,a);
    p=(int**)malloc(sizeof(int*)*n);
    for (i=0;i<n;i++)
    {
        *(p+i) = (int*)malloc(sizeof(int)*n);
    }

    for (i=0;i<n;i++) {
        for (j = 0; j < n; j++) {
            if (i<l && j<l) {
                fscanf(fp, "%d", *(p + i) + j);
            }
            else{
                *(*(p+i)+j)=0;
            }
        }
    }
    q=(int**)malloc(sizeof(int*)*n);
    for (i=0;i<n;i++)
    {
        *(q+i) = (int*)malloc(sizeof(int)*n);
    }
    for (i=0;i<n;i++) {
        for (j = 0; j < n; j++) {
            if (i<l && j<l) {
                fscanf(fp, "%d", *(q + i) + j);
            }
            else{
                *(*(q+i)+j)=0;
            }
        }
    }

    fclose(fp);
    nat_start=clock();
    c=native(l,p,q);
    nat_end=clock();
    nat_running_time=nat_end-nat_start;
    printf("%d ",nat_running_time);
    fp=fopen("output_p1_m1.txt","w+");
    for (i=0;i<l;i++) {
        for (j = 0; j < l; j++) {
            if (j == l - 1) {
                fprintf(fp, "%d \n", c[i][j]);
            }
            else
                fprintf(fp, "%d ", c[i][j]);
        }
    }
    fclose(fp);
    free(c);

    str_start=clock();
    c=strassen(n,p,q);
    str_end=clock();
    str_running_time=str_end-str_start;
    printf("%d",str_running_time);
    fp=fopen("output_p1_m2.txt","w+");
    for (i=0;i<l;i++) {
        for (j = 0; j < l; j++) {
            if (j == l - 1) {
                fprintf(fp, "%d \n", c[i][j]);
            }
            else
                fprintf(fp, "%d ", c[i][j]);
        }
    }
    fclose(fp);
    free(c);

    fp=fopen("output_p2.txt","w+");
    fprintf(fp,"%d ",nat_running_time);
    fprintf(fp,"%d ",str_running_time);
    fclose(fp);
   return 0;

}

