#include <stdio.h>
#include <complex.h>
#include "NUFFT2.c"
// #include "mat.h"

# define K2 16900
#define len_fn 459
#define len_kloc 16256
int main()

{

    FILE *fp1r;
    FILE *fp1i;
    FILE *fp2;
    FILE *fp3r;
    FILE *fp3i;
    double sk_r[len_kloc];
    double sk_i[len_kloc];
    double KFFT_r[K2];
    double KFFT_i[K2];
    int K=130;
    int N=128;
    int J=6;
    int Ofactor=151;
    int len=len_kloc;

    double fn[len_fn];
    // double complex kloc1[len_kloc];
    double kloc_r[len_kloc];
    double kloc_i[len_kloc];
    int i, j;

    // Open the file
    fp1r = fopen("KFFTimag_C_r.txt", "r");
    fp1i=fopen("KFFTimag_C_i.txt","r");
    if (fp1r == NULL) {
        printf("Error opening file.\n");
        return 1;
    }


    //Read the KFFT_imag_C from the file
    for (i = 0; i < K2; i++) {
        printf("%d",i);
        fscanf(fp1r, "%lE", &KFFT_r[i]);
        fscanf(fp1i, "%lE", &KFFT_i[i]);
    }

    // Close the file
    fclose(fp1r);
    fclose(fp1i);

    // Read fn
    fp2=fopen("fn.txt","r");

    for (i=0;i<len_fn;i++) {
        fscanf(fp2,"%lf", &fn[i]);
    }

    fclose(fp2);

    //Read kloc
    fp3r=fopen("kloc1_r.txt","r");
    fp3i=fopen("kloc1_i.txt","r");
    for (i=0;i<len_fn;i++) {
        fscanf(fp3r, "%lf", &kloc_r[i]);
        fscanf(fp3i, "%lf", &kloc_i[i]);
    }

    fclose(fp3r);
    fclose(fp3i);


    // //Display the matrix
    // printf("Matrix read from file:\n");
    // for (i = 0; i < K2; i++) {
    //     printf("%e", KFFT_r[i]);
    //     printf("\n");
    // }

    gridlut(KFFT_r,KFFT_i,K,N,J,kloc_r,kloc_i,sk_r,sk_i,fn,Ofactor,len);



return 0;

}
