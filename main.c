#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include <ctype.h>


//declaration of parameters that will be read further from the file "parameters.txt"
int N_x;
int N_z;
double x_begin, x_end, z_begin, z_end;
double n_0;
double pol; // Polarization type: 1 for 'Horz.' or 0 for 'Vert.'
double source_height, gamma_elv, gamma_beamwidth;
double a_0, source_frequency;
double b,rho, V_0, h_0, V, gamma_, C_, H;
double a, b_e, c, d, t, p_0, g, R, M, h_sea;

//standard
void standard_refraction(complex double* refractive_index, double T, double* P, double* e, int l){
    refractive_index [l]  = cpow(1 + (77.6*P[l]/T - 5.6*e[l]/T + 3.75*1.e5*e[l]/(T*T))*1.e-6, 2);//77.6*P[l]/T[l] - 5.6*e[l]/T[l] + 3.75*1.e5*e[l]/(T[l]*T[l]); // //it is square refractive index n^2
}


//looyeng_model
void looyeng_model(complex double* eps, complex double* refractive_index, complex double eps_1, complex double* phi_1,  double current_x, int l){
    phi_1[l] =  C_/(pow(V, gamma_)*rho);
    eps[l] = cpow((phi_1[l]*(cpow(eps_1,1.0/3.0)-cpow(refractive_index[l],1.0/3.0)) + cpow(refractive_index[l],1.0/3.0)), 3); //it is square refractive index n^2
}

void tridiag_matrix_algorithm(complex double* array_A, complex double* array_B, complex double* array_C, complex double* array_D,complex double* array_u){
    complex double *alpha = malloc(N_x * sizeof(complex double));
    complex double *beta = malloc(N_x * sizeof(complex double));

    alpha[0] = -array_A[0]/(array_B[0]);
    beta[0] = array_D[0]/(array_B[0]);


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


//ducting model
//void duct_refraction(complex double* refractive_index, double current_x,int l){
//
//    refractive_index [l] = 1 - a_0*current_x; //it is square refractive index n^2
//
//}


//kharadly-Jackson model
//void kharadly_and_jackson_model(complex double* eps, complex double* refractive_index, complex double eps_1, double* phi_1,  double current_x, int l){ //double* phi_1,
//    if (l==0)
//        phi_1[l] = 0;
//    else
//        phi_1[l] = C_*pow(h_0/current_x,b)/(rho*pow(V_0,gamma_));
//
//    eps[l] = (2*refractive_index[l]*phi_1[l]*(eps_1-refractive_index[l]) + refractive_index[l]*eps_1 + 2*refractive_index[l]*refractive_index[l])/(eps_1+2*refractive_index[l]-phi_1[l]*(eps_1-refractive_index[l])); //it is square refractive index n^2
//
//}

//Wagner model
//void wagner_model(complex double* refractive_index, complex double* eps_0, complex double eps_1, complex double* phi_1, double current_x,int l){
//    phi_1[l] = C_*pow(h_0/current_x,b)/(rho*pow(V_0,gamma_));
//    refractive_index [l] = 3*eps_0[l]*phi_1[l]*(eps_1-eps_0[l])/(eps_1+2*eps_0[l]) + eps_0[l]; //it is square refractive index n^2
//
//}

void remove_comments(const char *input_file, const char *output_file) {
    FILE *infile = fopen(input_file, "r");
    FILE *outfile = fopen(output_file, "w");
    if (infile == NULL || outfile == NULL) {
        perror("Error opening file");
        return;
    }

    char line[256];
    while (fgets(line, sizeof(line), infile)) {
        char *start = strstr(line, "/*");
        char *end;

        while (start != NULL) {
            end = strstr(start, "*/");
            if (end != NULL) {
                // Shift the rest of the line to the left
                memmove(start, end + 2, strlen(end + 2) + 1);
            } else {
                // If there is no closing comment, just cut the line
                *start = '\0';
                break;
            }
            start = strstr(line, "/*"); // Looking for next comment
        }

        //Write the cleared string to the output file
        fprintf(outfile, "%s", line);
    }

    fclose(infile);
    fclose(outfile);
}

int main() {
    clock_t start = clock();
    const char *input_file = "parameters.txt";
    const char *output_file = "parameters_cleaned.txt";

    remove_comments(input_file, output_file);

    if ((output_file = fopen("parameters_cleaned.txt", "r")) == NULL)
        printf("Error!\n");
    else {

        fscanf(output_file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &N_x, &N_z,
               &x_begin, &x_end,
               &z_begin, &z_end,
               &n_0,
               &pol,
               &source_height,
               &gamma_elv,
               &gamma_beamwidth,
               &a_0,
               &source_frequency,
               &b,
               &rho,
               &V_0,
               &h_0,
               &V,
               &gamma_,
               &C_,
               &H,
               &a,
               &b_e,
               &c,
               &d,
               &t,
               &p_0,
               &g,
               &R,
               &M,
               &h_sea);
        fclose(output_file);
    }


    double dx = (x_end - x_begin) / N_x, dz = (z_end - z_begin) / N_z;
    complex double k_0 = 2.0*M_PI*source_frequency/3.e8;
    complex double B = -1.0 /(dx*dx) + (2.0*I*k_0)/dz + k_0*k_0*(n_0*n_0-1.0)/2;
    complex double A = 1.0/(2.0*dx*dx), C = 1.0/(2.0*dx*dx);
    complex double eps_1 = (4.56 + 0.04*H - 7.78*pow(H,2)*1.e-4 + 5.56*pow(H,3)*1.e-6) + I*(0.251+ 0.02*H - 3.71*pow(H,2)*1.e-4 + 2.76*pow(H,3)*1.e-6);
    complex double **array_u = malloc(N_z * sizeof(complex double*)); //the number of lines is multiplied by sizeof
    complex double **refractive_index = malloc(N_z * sizeof(complex double*));
    complex double **eps = malloc(N_z * sizeof(complex double*));
    complex double *array_A = malloc(N_x * sizeof(complex double));
    complex double *array_B = malloc(N_x * sizeof(complex double));
    complex double *array_C = malloc(N_x * sizeof(complex double));
    complex double *array_D = malloc(N_x * sizeof(complex double));
    complex double *phi_1 = malloc(N_x * sizeof(complex double));
    double *P = malloc(N_x * sizeof(double));
    double *e = malloc(N_x * sizeof(double));


    for (int i = 0; i < N_z; ++i) {

        array_u[i] = malloc(N_x * sizeof(complex double)); //the number of columns is multiplied by sizeof
        refractive_index[i]= malloc(N_x * sizeof(complex double));
        eps[i] = malloc(N_x * sizeof(complex double));

    }


    //calculation of water vapor pressure e and pressure at current altitude P
    double T = t + 273.15;
    double current_x_0;
    for (int j = 0; j < N_x; j++) {
        current_x_0 = x_begin + j*dx;
        P[j] = p_0*exp(-M*g*(current_x_0-h_sea)/(R*T));
        e[j] = H / 100 * a * exp((b_e - t / d) * t / (t + c))*(1+1.e-4*(7.2+P[j]*(0.032 + 5.9*1.e-6*t*t)));
    }

    int h_max = (int)(0.75 * N_x);

// set initial values
    double source_omega = sqrt(2.0*log(2))/(k_0*sin(gamma_beamwidth/2));
    int counter = 0;

    double current_x_1;
    for (int j = 1; j < N_x-1; j++) {

        current_x_1 = x_begin + j*dx;
        array_u[0][j] = cexp(I*k_0*current_x_1*csin(gamma_elv)-pow(current_x_1-source_height,2)/pow(source_omega,2)) + pow(-1,pol)*cexp(I*k_0*(-current_x_1)*csin(gamma_elv)-pow(-current_x_1-source_height,2)/pow(source_omega,2));
    }

    for (int i = 0; i < N_z; i++) {
        array_u[i][0]=0;
        array_u[i][N_x-1]=0;
    }

    double current_z;
    double current_x;

    for (int k = 1; k < N_z; k++) {

        current_z = z_begin + dz * k;

        for (int l = 0; l < N_x; l++) {


            current_x = x_begin + dx * l;

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
                looyeng_model(eps[k], refractive_index[k], eps_1, phi_1, current_x, l);
                //duct_refraction(refractive_index[k], current_x,l);
                //exponential_refraction(refractive_index[k], current_x, l);
                //kharadly_and_jackson_model(eps[k], refractive_index[k], eps_1, phi_1, current_x, l);
                //wagner_model(refractive_index[k],eps_0[k], eps_1, phi_1[k], current_x, l);
//                array_B[l] = -1.0 / (dx * dx) + (2.0 * I * k_0) / dz + k_0 * k_0 * (refractive_index[k][l] - 1.0)/2;
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

        //Hanning window
        double current_xh;
        for (int ind = h_max; ind < N_x; ind++) {
            current_xh = x_begin + dx * ind;
            array_u[k][ind] *= csin(2 * M_PI * current_xh / x_end) * csin(2 * M_PI * current_xh / x_end);
        }

    }


    FILE *file;
    file = fopen("output.bin", "wb");

    if (file == NULL) {
        printf("FileIsNull\n");
        return 1;
    }

    for (int i = 0; i < N_z; i++) {
        for (int j = 0; j < N_x; j++) {
            double abs_array_u = cabs(array_u[i][j]); // вычисление модуля
            fwrite(&abs_array_u, sizeof(double), 1, file); // запись модуля в файл
        }
    }

    fclose(file);


//
//    FILE *file;
//    file = fopen("output.txt", "w+");
//
//    if (file == NULL) {
//        printf("FileIsNull\n");
//        return 1;
//    }
//
//    for (int i = 0; i<N_z; i++) {
//        for (int j = 0; j<N_x; j++){
//            fprintf(file, "%1.6e ", cabs(array_u[i][j]));
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



