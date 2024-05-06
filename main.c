#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<math.h>
#include<time.h>
#include <C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h>
#include <C:\Program Files\JetBrains\CLion 2023.2.2\bin\mingw\lib\gcc\x86_64-w64-mingw32\13.1.0\include\omp.h>



//SYSTEM SI
int N_x=12000;
int N_z=9000;
double x_begin=0,x_end=4000.0, z_begin=0,z_end=50000.0;
double n_0 = 1.00028;
double pol=1; // Polarization type: 1 for 'Horz.' or 0 for 'Vert.'

// source parameters
double source_height = 10.0;
double gamma_horiz=5*3.14/180; //elv
double gamma_rastvor=0.5*3.14/180; //bw
double a_0 = 1.2e-6;//2.4e-6;
double source_frequency = 3.e9;



//standard
void standard_refraction(complex double* refractive_index, double current_x,int l){

    refractive_index [l] =  pow ((1 + (315*cexp(-1.36 *1.e-5*current_x))*1.e-6),2); //it is square refractive index n^2

}


//ducting
void duct_refraction(complex double* refractive_index, double current_x,int l){

    refractive_index [l] = 1 - a_0*current_x; //it is square refractive index n^2

}


//exponential function of refractive index
void exponential_refraction(complex double* refractive_index, double current_x,int l){

    refractive_index [l] = pow ((1 + (315*cexp(-current_x/7350))*1.e-6),2);  //it is square refractive index n^2
}


//linear function of refractive index
void linear_refraction(complex double* refractive_index, double current_x,int l){

    refractive_index [l] = pow((1 + (315-(-0.0392)*(current_x-7350))*1.e-6),2); //it is square refractive index n^2

}


void tridiag_matrix_algorithm(complex double* array_A, complex double* array_B, complex double* array_C, complex double* array_D,complex double* array_u) {
    complex double *alpha = malloc(N_x * sizeof(complex double));
    complex double *beta = malloc(N_x * sizeof(complex double));

    alpha[0] = -array_A[0]/(array_B[0]);
    beta[0] = array_D[0]/(array_B[0]);



    for (int i=1; i<N_x; i++) {

        alpha[i] = -array_A[i]/(array_B[i]+array_C[i]*alpha[i-1]);
        beta[i] = (array_D[i]-array_C[i]*beta[i-1])/(array_B[i]+array_C[i]*alpha[i-1]);

    }

    array_u[N_x-1]=beta[N_x-1];

    for (int i = N_x - 1; i >= 0; i--) {
        array_u[i] = alpha[i] * array_u[i+1] + beta[i];

    }

    free(alpha);
    free(beta);
}


int main() {
    clock_t start = clock();
    double dx = (x_end - x_begin) / N_x, dz = (z_end - z_begin) / N_z;
    complex double k_0 = 2.0*M_PI*source_frequency/3.e8;
    complex double B = -1.0 /(dx*dx) + (2.0*I*k_0)/dz + k_0*k_0*(n_0*n_0-1.0)/2;
    complex double A = 1.0 /(2.0*dx*dx), C = 1.0 /(2.0*dx*dx);

    FILE *file = NULL;
    complex double **array_u = malloc(N_z * sizeof(complex double*)); //the number of lines is multiplied by sizeof
    complex double **refractive_index = malloc(N_z * sizeof(complex double*));
    complex double *array_A = malloc(N_x * sizeof(complex double));
    complex double *array_B = malloc(N_x * sizeof(complex double));
    complex double *array_C = malloc(N_x * sizeof(complex double));
    complex double *array_D = malloc(N_x * sizeof(complex double));

    for (int i = 0; i < N_z; ++i) {

        array_u[i] = malloc(N_x * sizeof(complex double)); //the number of columns is multiplied by sizeof
        refractive_index[i]= malloc(N_x * sizeof(complex double));


    }

    double maximum_u=0;

// set initial values
    double source_omega = sqrt(2.0*log(2))/(k_0*sin(gamma_rastvor/2));
    double gamma_horiz_1 = 5*3.14/180;
    double gamma_horiz_2 = 10*3.14/180;
    double gamma_horiz_3 = 45*3.14/180;
    double gamma_horiz_4 = 15*3.14/180;
    double gamma_horiz_5 = 20*3.14/180;
    int counter = 0;



    for (int j = 1; j < N_x-1; j++) {
        double current_x_1 = x_begin + j*dx;
        array_u[0][j] = cexp(I*k_0*current_x_1*csin(gamma_horiz)-pow(current_x_1-source_height,2)/pow(source_omega,2)) + pow(-1,pol)*cexp(I*k_0*(-current_x_1)*csin(gamma_horiz)-pow(-current_x_1-source_height,2)/pow(source_omega,2));
//                + cexp(I*k_0*current_x_1*csin(gamma_horiz_1)-pow(current_x_1-source_height,2)/pow(source_omega,2)) + cpow(-1,pol)*cexp(I*k_0*(-current_x_1)*csin(gamma_horiz_1)-pow(-current_x_1-source_height,2)/pow(source_omega,2))
//               +cexp(I*k_0*current_x_1*csin(gamma_horiz_2)-pow(current_x_1-source_height,2)/pow(source_omega,2)) + cpow(-1,pol)*cexp(I*k_0*(-current_x_1)*csin(gamma_horiz_2)-pow(-current_x_1-source_height,2)/pow(source_omega,2))
//                +cexp(I*k_0*current_x_1*csin(gamma_horiz_3)-pow(current_x_1-source_height,2)/pow(source_omega,2)) + cpow(-1,pol)*cexp(I*k_0*(-current_x_1)*csin(gamma_horiz_3)-pow(-current_x_1-source_height,2)/pow(source_omega,2))
//                +cexp(I*k_0*current_x_1*csin(gamma_horiz_4)-pow(current_x_1-source_height,2)/pow(source_omega,2)) + cpow(-1,pol)*cexp(I*k_0*(-current_x_1)*csin(gamma_horiz_4)-pow(-current_x_1-source_height,2)/pow(source_omega,2))
//                +cexp(I*k_0*current_x_1*csin(gamma_horiz_5)-pow(current_x_1-source_height,2)/pow(source_omega,2)) + cpow(-1,pol)*cexp(I*k_0*(-current_x_1)*csin(gamma_horiz_5)-pow(-current_x_1-source_height,2)/pow(source_omega,2));
        //printf("%0.12f + i*(%0.12f)\n", creal(array_u[0][j]), cimag(array_u[0][j]));
        // if (cabs(array_u[0][j])<1)
        //maximum_u = fmax(cabs(array_u[0][j]),maximum_u);
        // printf("endof: %f\n",maximum_u);

    }

    for (int i = 0; i < N_z; i++) {
        array_u[i][0]=0;
        array_u[i][N_x-1]=0;

    }



    for (int k = 1; k < N_z; k++) {

        double current_z = z_begin + dz * k;

        for (int l = 0; l < N_x; l++) {

            double current_x = x_begin + dx * l;

            //recalculation A,B,C,D
            if (l == 0) {
                array_D[l] = 0;
                array_B[l] = B;
                array_A[l] = 0;
                array_C[l] = 0;
            } else if (l == N_x - 1) {
                array_D[l] = 0;
                array_B[l] = B;
                array_A[l] = 0;
                array_C[l] = 0;
            } else {
                // array_B[l]=B;
                array_A[l] = A;
                array_C[l] = C;
                standard_refraction(refractive_index[k], current_x, l);
                //exponential_refraction(refractive_index[k], current_x, l);
                //duct_refraction(refractive_index[k], current_x,l);
                //linear_refraction(refractive_index[k], current_x,l);
                //printf("%10.7e ",cabs(refractive_index[k][300]));
                array_B[l] = -1.0 / (dx * dx) + (2.0 * I * k_0) / dz + k_0 * k_0 * (refractive_index[k][l] - 1.0);
                array_D[l] = array_u[k - 1][l] * (2.0 * I * k_0 / dz + 1.0 / (dx * dx) -
                                                  k_0 * k_0 * (refractive_index[k][l] - 1.0) / 2) -
                             1.0 / (2.0 * dx * dx) * (array_u[k - 1][l + 1] + array_u[k - 1][l - 1]);
                if (cabs(array_B[l]) >= cabs(array_A[l]) + cabs(array_C[l]))
                    counter += 0;
                else
                    counter += 1;
            }

            //printf("ee: %f\n", cabs(array_D[l]));


        }


        tridiag_matrix_algorithm(array_A, array_B, array_C, array_D, array_u[k]);


        int h = round(0.85 * (N_x));
        //printf("%d\n",h);
        //Hanning window
        for (h; h < N_x; h++) {
            double current_x = x_begin + dx * h;
            array_u[k][h] *= csin(2 * M_PI * current_x / x_end) * csin(2 * M_PI * current_x / x_end);
        }

    }





    printf("hello \n");




    file = fopen("5_tilt.txt", "w+");

    if (file == NULL) {
        printf("FileIsNull\n");
        return 1;
    }

    for (int i = 0; i<N_z; i++) {
        for (int j = 0; j<N_x; j++){

            fprintf(file, "%1.10e ", cabs(array_u[i][j]));
        }

        fprintf(file,"\n");
    }


    fclose(file);



    for (int i = 0; i < N_z; i++) {
        free(array_u[i]);
        free(refractive_index[i]);
    }

    free(array_u);
    free(refractive_index);
    free(array_A);
    free(array_B);
    free(array_C);
    free(array_D);

    printf("Hello, World! : %d \n", counter);
    clock_t end = clock();
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Execution time: %f s\n", elapsed_time);
    return 0;
}



