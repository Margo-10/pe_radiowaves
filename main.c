#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<math.h>
#include<time.h>
//#include <C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h>
//#include <C:\Program Files\JetBrains\CLion 2023.2.2\bin\mingw\lib\gcc\x86_64-w64-mingw32\13.1.0\include\omp.h>



//SYSTEM SI
int N_z=6000;
int N_x=6750;
double x_begin=0,x_end=2700.0, z_begin=0,z_end=30000.0;
double n_0 = 1.00028;
double pol=1; // Polarization type: 1 for 'Horz.' or 0 for 'Vert.'

// source parameters
double source_height = 30.0;
double gamma_horiz=5*M_PI/180; //elv
double gamma_rastvor=2*M_PI/180; //bw
double a_0 = 1.2e-6;
double source_frequency = 3.e9;
//complex double eps_1 = 4.56+I*0.251;
double b=0.28;
double rho=2440.0;

//visibility
double V_0 = 6.5; //meters
double h_0 = 2.0;


double V = 100;

//for Libya and Sudan
double gamma_ = 1.07;
double C_ = 2.3*1.e-2;
double Humidity = 40;


//standard
void standard_refraction(complex double* refractive_index, double T, double* P, double* e, int l){

    refractive_index [l]  = cpow(1 + (77.6*P[l]/T - 5.6*e[l]/T + 3.75*1.e5*e[l]/(T*T))*1.e-6, 2);//77.6*P[l]/T[l] - 5.6*e[l]/T[l] + 3.75*1.e5*e[l]/(T[l]*T[l]); // //it is square refractive index n^2

}


//looyeng_model
void looyeng_model(complex double* eps, complex double* refractive_index, complex double eps_1, complex double* phi_1,  double current_x, int l){
//    if (current_x<2.0)
//        phi_1[l] = C_*pow(h_0/h_0,b)/(rho*pow(V_0,gamma_));
//    else
//        phi_1[l] = C_*pow(h_0/current_x,b)/(rho*pow(V_0,gamma_));


    phi_1[l] =  C_/(pow(V, gamma_)*rho);
    eps[l] = cpow((phi_1[l]*(cpow(eps_1,1.0/3.0)-cpow(refractive_index[l],1.0/3.0)) + cpow(refractive_index[l],1.0/3.0)), 3); //it is square refractive index n^2
}

void tridiag_matrix_algorithm(complex double* array_A, complex double* array_B, complex double* array_C, complex double* array_D,complex double* array_u){
    complex double *alpha = malloc(N_x * sizeof(complex double));
    complex double *beta = malloc(N_x * sizeof(complex double));

    alpha[0] = -array_A[0]/(array_B[0]);
    beta[0] = array_D[0]/(array_B[0]);

    //printf("hello \n");

    for (int i=1; i<N_x; i++) {

        alpha[i] = -array_A[i]/(array_B[i]+array_C[i]*alpha[i-1]);
        beta[i] = (array_D[i]-array_C[i]*beta[i-1])/(array_B[i]+array_C[i]*alpha[i-1]);

    }

    array_u[N_x-1]=beta[N_x-1];

    for (int i = N_x - 2; i > 0; i--) {
        array_u[i] = alpha[i] * array_u[i+1] + beta[i];

    }

    free(alpha);
    free(beta);
}




//ducting
//void duct_refraction(complex double* refractive_index, double current_x,int l){
//
//    refractive_index [l] = 1 - a_0*current_x; //it is square refractive index n^2
//
//}


//kharadly_and_jackson_model
//void kharadly_and_jackson_model(complex double* eps, complex double* refractive_index, complex double eps_1, double* phi_1,  double current_x, int l){ //double* phi_1,
//    if (l==0)
//        phi_1[l] = 0;
//    else
//        phi_1[l] = C_*pow(h_0/current_x,b)/(rho*pow(V_0,gamma_));
//
//    eps[l] = (2*refractive_index[l]*phi_1[l]*(eps_1-refractive_index[l]) + refractive_index[l]*eps_1 + 2*refractive_index[l]*refractive_index[l])/(eps_1+2*refractive_index[l]-phi_1[l]*(eps_1-refractive_index[l])); //it is square refractive index n^2
//
//}

//wagner
//void wagner_model(complex double* refractive_index, complex double* eps_0, complex double eps_1, complex double* phi_1, double current_x,int l){
//    phi_1[l] = C_*pow(h_0/current_x,b)/(rho*pow(V_0,gamma));
//    refractive_index [l] = 3*eps_0[l]*phi_1[l]*(eps_1-eps_0[l])/(eps_1+2*eps_0[l]) + eps_0[l]; //it is square refractive index n^2
//
//}


int main() {

    clock_t start = clock();
    double dx = (x_end - x_begin) / N_x, dz = (z_end - z_begin) / N_z;
    complex double k_0 = 2.0*M_PI*source_frequency/3.e8;
    complex double B = -1.0 /(dx*dx) + (2.0*I*k_0)/dz + k_0*k_0*(n_0*n_0-1.0)/2;
    complex double A = 1.0/(2.0*dx*dx), C = 1.0/(2.0*dx*dx);
    complex double eps_1 = (4.56 + 0.04*Humidity - 7.78*pow(Humidity,2)*1.e-4 + 5.56*pow(Humidity,3)*1.e-6) + I*(0.251+ 0.02*Humidity - 3.71*pow(Humidity,2)*1.e-4 + 2.76*pow(Humidity,3)*1.e-6);


    FILE *file = NULL;
    complex double **array_u = malloc(N_z * sizeof(complex double*)); //the number of lines is multiplied by sizeof
    complex double **refractive_index = malloc(N_z * sizeof(complex double*));
    complex double **eps = malloc(N_z * sizeof(complex double*));
    complex double *array_A = malloc(N_x * sizeof(complex double));
    complex double *array_B = malloc(N_x * sizeof(complex double));
    complex double *array_C = malloc(N_x * sizeof(complex double));
    complex double *array_D = malloc(N_x * sizeof(complex double));
    complex double *phi_1 = malloc(N_x * sizeof(complex double));
//    double *T = malloc(N_x * sizeof(double));
    double *P = malloc(N_x * sizeof(double));
    double *e = malloc(N_x * sizeof(double));


    for (int i = 0; i < N_z; ++i) {

        array_u[i] = malloc(N_x * sizeof(complex double)); //the number of columns is multiplied by sizeof
        refractive_index[i]= malloc(N_x * sizeof(complex double));
        eps[i] = malloc(N_x * sizeof(complex double));

    }

    double a =  6.1121;
    double b_e =  18.678;
    double c =  257.14;
    double d =  234.5;
    double t = 20; //Celsius
    double T = t + 273.15;
    double p_0 = 1013.25; //gPa
    double g = 9.81;
    double R = 8.31;
    double M = 0.02897; // kg/mol;
    double h_sea = 0;


    for (int j = 0; j < N_x; j++) {
        double current_x = x_begin + j*dx;
        P[j] = p_0*exp(-M*g*(current_x-h_sea)/(R*T));
        e[j] = Humidity / 100 * a * exp((b_e - t / d) * t / (t + c))*(1+1.e-4*(7.2+P[j]*(0.032 + 5.9*1.e-6*t*t)));
    }

    int h_max = (int)(0.75 * N_x);
    double maximum_u=0;
    double U_up = 0;
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
            if ((l == 0) || (l == N_x - 1)) {
                array_D[l] = 0;
                array_B[l] = B;
                array_A[l] = 0;
                array_C[l] = 0;
            }
            else {
                array_A[l] = A;
                array_C[l] = C;
                standard_refraction(refractive_index[k], T, P, e, l);

                //duct_refraction(refractive_index[k], current_x,l);
                looyeng_model(eps[k], refractive_index[k], eps_1, phi_1, current_x, l);
                //linear_refraction(refractive_index[k], current_x,l);
                //exponential_refraction(refractive_index[k], current_x, l);
                //kharadly_and_jackson_model(eps[k], refractive_index[k], eps_1, phi_1, current_x, l);
                //wagner_model(refractive_index[k],eps_0[k], eps_1, phi_1[k], current_x, l);
                //printf("%10.7e ",cabs(eps_0[k][300]));
//                array_B[l] = -1.0 / (dx * dx) + (2.0 * I * k_0) / dz + k_0 * k_0 * (refractive_index[k][l] - 1.0);
//                array_D[l] = array_u[k - 1][l] * (2.0 * I * k_0 / dz + 1.0 / (dx * dx) -
//                                                  k_0 * k_0 * (refractive_index[k][l] - 1.0) / 2) -
//                             1.0 / (2.0 * dx * dx) * (array_u[k - 1][l + 1] + array_u[k - 1][l - 1]);

                array_B[l] = -1.0 / (dx * dx) + (2.0 * I * k_0) / dz + k_0 * k_0 * (eps[k][l] - 1.0)/2;
                array_D[l] = array_u[k - 1][l] * (2.0 * I * k_0 / dz + 1.0 / (dx * dx) -  k_0 * k_0 * (eps[k][l] - 1.0) / 2) - 1.0 / (2.0 * dx * dx) * (array_u[k - 1][l + 1] + array_u[k - 1][l - 1]);
                if (cabs(array_B[l]) >= (cabs(array_A[l]) + cabs(array_C[l])))
                    counter += 0;
                else
                    counter += 1;

            }

        }

        tridiag_matrix_algorithm(array_A, array_B, array_C, array_D, array_u[k]);
//        printf("%.10f + i%.10f\n", creal(array_u[2][251]), cimag(array_u[2][251]));
//        printf("%10.7e ",cabs(array_u[2][251]));

        //Hanning window
        for (int ind = h_max; ind < N_x; ind++) {
            double current_xh = x_begin + dx * ind;
            array_u[k][ind] *= csin(2 * M_PI * current_xh / x_end) * csin(2 * M_PI * current_xh / x_end);
        }

    }


    double max_k = 0;
    for (int k = 0; k < N_x; k++) {
            double current_value = cabs(array_u[5999][9922]);
            if (current_value > U_up) {
                U_up = current_value;
                max_k = k;
        }
    }
    max_k = z_begin + dz * max_k;
    printf(": %f :%f  ",U_up, max_k); //Расстояние: %f grbg %fmax_k, maximum_u


//
//    file = fopen("field3ghz_satellite_5tilt_30vid_20hum.txt", "w+");
//
//    if (file == NULL) {
//        printf("FileIsNull\n");
//        return 1;
//    }
//
//    for (int i = 0; i<N_z; i++) {
//        for (int j = 0; j<N_x; j++){
//            double current_z = z_begin + dz * i;
//            fprintf(file, "%1.6e ", cabs(array_u[i][j]));
//            //fprintf(file, "%1.13e ", -20*log(cabs(array_u[i][j]))+ 20*log(4*M_PI) + 10*log(current_z)-30*log(3.e8/source_frequency));
//        }
//
//        fprintf(file,"\n");
//    }
//
//
//    fclose(file);


    for (int i = 0; i < N_z; i++) {
        free(array_u[i]);
        free(refractive_index[i]);
        free(eps[i]);
    }

    free(array_u);
    free(refractive_index);
    free(eps);
    free(array_A);
    free(array_B);
    free(array_C);
    free(array_D);
    free(phi_1);
//    free(T);
    free(P);
    free(e);

    printf("Convergence check : %d ", counter);
    clock_t end = clock();
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Execution time: %f ", elapsed_time);
    return 0;
}



