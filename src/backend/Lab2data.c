#include <stdio.h>
#include <math.h>
#include "Lab2data.h"
#define MAX(a, b) a > b ? a : b
#define DINT(a) (double)(int)a

double phi_func(double x) {
    if (x > 0.0) {
        if (x > 45.0) {
            return 1.0;
        } else {
            double temp = (((((x * 5.383e-6 + 4.88906e-5) * x + 3.80036e-5) * x + 
                  .0032776263) * x + .0211410061) * x + .049867347) * x + 1.;
            temp *= temp;
            temp *= temp;
            temp *= temp;
            return (1. - .5 / (temp * temp));                  
        }
    } else {
        if (x < -45.0) {
            return 0.0;
        } else {
            x = -x;
            return (0.50 / pow((1.0 + (0.0498673470 +
                    (0.0211410061 + (0.0032776263 + (0.0000380036 +
                    (0.0000488906 + 0.0000053830 * x) * x) * x) * x) * x) *
                    x), 16));
        }
    }
}

double phi_quantile(double beta) {
    double delta;
    if (beta >= 0.5) {
        delta = beta;
    } else {
        delta = 1.0 - beta;
    }
    double w = sqrt(-log((1.0 - delta) * (1.0 - delta)));
    double low = w - (2.515517 + (0.802853 + 0.010328 * w) * w) /
                 (1.0 + (1.432788 + (0.189269 + 0.001308 * w) * w) * w);
    double high = low + 4.5E-4;
    low -= 4.5E-4;
    double mid;
    do {
        mid = (low + high) / 2.0;
        if (phi_func(mid) > delta) {
            low = mid;
        } else {
            high = mid;
        }
    } while (fabs(mid - (low + high) / 2.0) > 1.0E-15);
    if (beta >= 0.5) {
        return mid;
    } else {
        return -mid;
    }
}

double student_t_quantile(double beta, int l) {
    double delta;
    double dfl;
    double high;
    double low;
    double mid;
    double test;
    double w;
    double x;
    double y;
    if (beta == 0.5) {
        return 0.0;
    } else if (l >= 7) {
        x = phi_quantile(beta);
        y = x * x;
        dfl = (double)l;
        return (x * (1.0 +
           ((-945.0 + (-3600.0 + (2880.0 + 23040.0 * dfl) * dfl) * dfl)
           +((-1920.0 + (4080.0 + (15360.0 + 23040.0 * dfl) * dfl) *
           dfl) + ((1482.0 + (4560.0 + 4800.0 * dfl) * dfl) + ((776.0 +
           720.0 * dfl) + 79.0 * y) * y) * y) * y)/(pow(dfl, 4) * 92160.0)));
    } else if (l == 1) {
        return tan((beta - 0.5) * 3.1415926535898);
    } else if (l == 2) {
        return ((2.0 * beta - 1.0) / sqrt(2.0 * beta * (1.0 - beta)));
    } else {
        if (beta >= 0.5) {
            delta = beta;
        } else {
            delta = 1.0 - beta;
        }
        low = phi_quantile(delta);
        high = (2.0 * delta - 1.0) / sqrt(2.0 * delta * (1.0 - delta));
        do {
            mid = (low + high) / 2.0;
            w = mid * mid;
            if (l == 3.0) {
                test = sqrt(3.0) * tan(3.1415926535898 *
                       (delta - 0.5)-(mid * sqrt(3.0))/(3.0 + w));
            } else if (l == 4.0) {
                test = (2.0 * delta - 1.0) * pow(sqrt(4.0 + w), 3) / (6.0 + w);
            } else if (l == 5.0) {
                test = sqrt(5.0) * tan(3.1415926535898 *
                       (delta - 0.5) - (mid * sqrt(5.0)) / (3.0 *
                       pow((5.0 + w), 2)) * (25.0 + 3.0 * w));     
            } else if (l == 6) {
                test = (2.0 * delta - 1.0) * pow(sqrt(6.0 + w), 5) /
                       (w * w + 15.0 * w + 67.5);
            }
        } while (fabs(mid - (low + high) / 2.0) > 1.0E-15);
        if (beta >= 0.5) {
            return mid;
        } else {
            return -mid;
        }
    }
}

void size(int t, int l_upper, int* b_1, int* b_2_sqrt, int* l_1,
          int* l_2_sqrt, int* t_prime) {
    int alpha = 0;
    int alpha_max;
    int i_max;
    int temp_1;
    int temp_2;
    int b[46] = {1,5,5,7,7,7,8,9,12,12,12,15,17,17,17,19,
                 21,21,22,22,26,27,29,31,32,33,36,36,36,
                 39,41,43,46,49,50,51,51,53,56,60,63,66,
                 67,70,84,85};
    int l[46] = {3,7,21,15,25,35,11,13,17,51,85,21,36,60,
                 84,27,25,35,31,93,37,57,41,66,45,47,51,
                 68,85,55,87,61,65,69,71,60,84,75,79,85,
                 89,93,95,99,85,96};
    i_max = 0;
    alpha_max = 0;
    temp_2 = 0;
    for (int i = 0; i < 46; i++) {
        if ((t >= b[i] * l[i]) && (l[i] <= l_upper)) {
            alpha = (int)(log(t / l[i] / b[i]) /
                          log(2.0));
            temp_1 = pow(2, alpha) * b[i] * l[i];
            if (temp_1 >= temp_2) {
                if ((temp_1 > temp_2) ||
                    (l[i] * b[i] > l[i_max] * b[i_max])) {
                    i_max = i;
                    alpha_max = alpha;
                    temp_2 = temp_1;
                }
            } 
        }
    }
    *b_1 = b[i_max];
    *b_2_sqrt = 3;
    if (*b_1 > 1) {
        *b_2_sqrt = (int)(sqrt(2.0) * *b_1 + 0.5);
    } 
    *l_1 = l[i_max];
    *l_2_sqrt = (int)(sqrt(2.0) * *l_1 + 0.5);
    *t_prime = temp_2;
}

void batch_updates(int j_dummy, int r_dummy, double s, double* w_vector,
                   double* s_vector, double* theta_vector, double* xi_vector,
                   int r_max) {
    //all r's must be changed to r - 1 in array indices
    int j = j_dummy;
    int r = r_dummy;
    double w;
    while (!(j % 2)) {
        w = s - s_vector[r - 1];
        xi_vector[r - 1] += w * w_vector[r - 1];
        w_vector[r - 1] = w;
        s_vector[r - 1] = s;
        j >>= 1;
        r += 2;
    }
    if (j == 1) {
        w = s;
        theta_vector[r - 1] = 0.5 * w * w;
        xi_vector[r - 1] = 0.0;
    } else {
        w = s - s_vector[r - 1];
        theta_vector[r - 1] += w * w_vector[r - 1];
    }
    w_vector[r - 1] = w;
    s_vector[r - 1] = s;
}

int write_tableau_file(double delta, int head, int s_num, int t, int t_prime,
                        double* beta_vector, double* ind_sum,
                        int* method_vector, int* row_vector,
                        int b_matrix[S_MAX][V_MAX], int l_matrix[S_MAX][V_MAX],
                        double p_matrix[S_MAX][V_MAX],
                        double v_matrix[S_MAX][V_MAX],
                        double x_bar_matrix[S_MAX][V_MAX], FILE* fout) {
    int b = 0;
    double h = 0.0;
    int l = 0;
    int lines = 0;
    double rel = 0.0;
    double p = 0.0;
    int row = 0;
    int s = 0;
    double v = 0;
    double x_bar = 0;
    if (head == 1) {
        perror("Not implemented");
        return 1;
    }
    fprintf(fout, "\n%s\n", "Final Tableau");
    for (int i = 0; i < 36; i++) {
        fprintf(fout, "%c", ' ');
    }
    fprintf(fout, "%s\n", "Mean Estimation");
    for (int i = 0; i < 36; i++) {
        fprintf(fout, "%c", ' ');
    }
    for (int i = 0; i < 15; i++) {
        fprintf(fout, "%c", '*');
    }
    fprintf(fout, "%c", 10);
    for (int i = 0; i < 36; i++) {
        fprintf(fout, "%c", ' ');
    }
    fprintf(fout, "(t = %9d)\n\n", t);
    for (int i = 0; i < 47; i++) {
        fprintf(fout, "%c", ' ');
    }
    fprintf(fout, "%4.1lf%%\n", (1.0 - delta) * 100.0);
    fprintf(fout, "%s", "            _       Standard Error      Confidence ");
    fprintf(fout, "%s\n", "Interval                    _");
    fprintf(fout, "%s", "Series      X      Sqrt[VAR.EST./t]      Lower");
    fprintf(fout, "%s\n\n", "       Upper      (Upper-Lower)/|X|");
    for (int series = 0; series < s_num; series++) {
        row = row_vector[series];
        b = b_matrix[series][row - 1];
        l = l_matrix[series][row - 1];
        s = l * b;
        if (s == t_prime) {
            s = t;
        }
        x_bar = ind_sum[series] / (double)s;
        v = v_matrix[series][row - 1];
        h = student_t_quantile(1.0 - delta / 2.0, l - 1) *
            sqrt(b * v / (double)s);
        rel = 0.0;
        if (x_bar) {
            rel = 2 * h / fabs(x_bar);
        }
        fprintf(fout, "%4d  %11.3e    ", series + 1, x_bar);
        fprintf(fout, " %11.3e    %11.3e  ", sqrt(b * v / s), x_bar - h);
        fprintf(fout, "%11.3e      %11.3e\n\n", x_bar + h, rel);
    }
    fprintf(fout, "%s\n", "   _");
    fprintf(fout, "%s\n", "   X is based on all t observations.");
    fprintf(fout, "   The Variance Estimator is based on first %6.2lf",
            l * b * 100.0 / t);
    fprintf(fout, "%s\n", "% of the t observations");
    for (int series = 0; series < s_num; series++) {
        row = row_vector[series];
        for (int i = 0; i < row; i++) {
            b = b_matrix[series][i];
            l = l_matrix[series][i];
            p = p_matrix[series][i];
            x_bar = x_bar_matrix[series][i];
            s = l * b;
            if (s == t_prime) {
                s = t;
                x_bar = ind_sum[series] / (double)s;
            }
            v = v_matrix[series][i];
            h = student_t_quantile(1.0 - delta / 2.0, l - 1) *
                sqrt(b * v / (double)s);
            fprintf(fout, "%3d %8d %8d %8d ", i + 1, b * l, l, b);
            fprintf(fout, "%11.3e %11.3e %11.3e ", x_bar, x_bar - h, x_bar + h);
            fprintf(fout, "%11.3e   %6.4f\n", sqrt((double)b * v), p);
        }
        l = l * b;
        b = 1;
        p = p_matrix[series][row];
        v = v_matrix[series][row];
        h = student_t_quantile(1.0 - delta / 2.0, s - 1) * sqrt(b * v / s);
        fprintf(fout, "    %8d %8d %8d ", t, t, b);
        fprintf(fout, "%11.3e %11.3e %11.3e ", x_bar, x_bar - h, x_bar + h);
        fprintf(fout, "%11.3e   %6.4lf\n", sqrt((double)b * v), p);
        if (beta_vector[series] < 0) {
            beta_vector[series] = 0;
        }
        fprintf(fout, "\n %5.2lf", beta_vector[series]);
        fprintf(fout, "%s", " significance leval for independence");
        fprintf(fout, "%s\n", " testing.");
        fprintf(fout, "  Review %3d", row);
        fprintf(fout, " used the first %6.2lf", l * b * 100.0 / t);
        fprintf(fout, "%s\n",
                "% of the t observations for variance estimation.");
    }
    return 0;
}

StatCode batch_means(FILE* fin, FILE* fout, int t, int s_num,
                     double* psi_vector, double delta, int rule, double beta,
                     int l_upper, int screen) {
    int b_1 = 0;
    int b_2_sqrt = 0;
    int head = 0;
    int l_1 = 0;
    int l_2_sqrt = 0;
    int t_prime = 0;
    int method_vector[S_MAX] = {0};
    double beta_vector[S_MAX] = {0.0};
    int row_vector[S_MAX] = {0};
    int b_matrix[S_MAX][V_MAX] = {{0}};
    double ind_sum[S_MAX] = {0.0};
    int l_matrix[S_MAX][V_MAX] = {{0}};
    double p_matrix[S_MAX][V_MAX] = {{0.0}};
    double rel[S_MAX] = {0.0};
    double v_matrix[S_MAX][V_MAX] = {{0.0}};
    double x_bar_matrix[S_MAX][V_MAX] = {{0.0}};
    int b_vector[S_MAX] = {0};
    int b_sqrt_vector[S_MAX] = {0};
    int j_vector[S_MAX] = {0};
    int j_sqrt_vector[S_MAX] = {0};
    int l_vector[S_MAX] = {0};
    int l_sqrt_vector[S_MAX] = {0};
    int r_sqrt_vector[S_MAX] = {0};
    int r_vector[S_MAX] = {0};
    int test_vector[S_MAX] = {0};
    double s_vector[S_MAX] = {0.0};
    double w_ind_vector[S_MAX] = {0.0};
    double y_ind_vector[S_MAX] = {0.0};
    double y_sqrt_vector[S_MAX] = {0.0};
    double y_vector[S_MAX] = {0.0};
    double z_ind_vector[S_MAX] = {0.0};
    //start
    double s_matrix[S_MAX][V_MAX] = {{0.0}};
    double theta_matrix[S_MAX][V_MAX] = {{0.0}};
    double w_matrix[S_MAX][V_MAX] = {{0.0}};
    double xi_matrix[S_MAX][V_MAX] = {{0.0}};
    //end change in index order due to row major in c vs column major in fortran
    double phi = 0;
    int b = 0;
    int b_sqrt = 0;
    int i = 0;
    int ii = 0;
    int j = 0;
    int j_sqrt = 0;
    int jj = 0;
    char kk = 0;
    int l = 0;
    int l_sqrt = 0;
    int old_row = 0;
    int r = 0;
    int r_max = 0;
    int row = 0;
    int r_sqrt = 0;
    int swap_int = 0;
    int ss = 0;
    int test = 0;
    double c = 0;
    double rho = 0;
    double s = 0;
    double swap_dbl = 0;
    double tau = 0;
    double x = 0;
    double y = 0;
    double y_sqrt = 0;
    head = 0;
    if (!t_prime) {
        for (int series = 0; series < s_num; series++) {
            method_vector[series] = rule;
            beta_vector[series] = beta;
        }
        size(t, l_upper, &b_1, &b_2_sqrt, &l_1, &l_2_sqrt, &t_prime);
    }
    do {
    if (fin) {
        for (int series = 0; series < s_num;
             fscanf(fin, "%lf", &psi_vector[series++]));
    }
    for (int series = 0; series < s_num; series++) {
        if (ss < s_num) {
            if (!ss) {
                if (s_num > S_MAX) {
                    perror("Maximum number of series exceeded");
                    return MAX_NUM_SERIES;
                }
                i = 1;
            }
            ss++;
            b = b_1;
            b_sqrt = b_2_sqrt;
            l = l_1;
            l_sqrt = l_2_sqrt;
            j = 0;
            y = 0.0;
            r = 1;
            j_sqrt = 0;
            y_sqrt = 0.0;
            r_sqrt = 2;
            r_max = MAX(2 * DINT(log((double)(t / b)) / log(2.0)) + 1,
                        2 * DINT(log((double)(t / b_sqrt)) / log(2.0)) + 2);
            if (r_max + 1 > V_MAX) {
                perror("Maximum vector size exceeded");
                return MAX_VECTOR_SIZE;
            }
            test = 1;
            s = 0;
            row_vector[series] = 0;
        } else {
            b = b_vector[series];
            l = l_vector[series];
            s = s_vector[series];
            b_sqrt = b_sqrt_vector[series];
            j = j_vector[series];
            j_sqrt = j_sqrt_vector[series];
            l_sqrt = l_sqrt_vector[series];
            r = r_vector[series];
            r_sqrt = r_sqrt_vector[series];
            test = test_vector[series];
            y = y_vector[series];
            y_sqrt = y_sqrt_vector[series];
        }
        
        x = psi_vector[series];
        if (i <= t_prime) {
            s += x;
        }
        ind_sum[series] += x;
        y_ind_vector[series] += x * x;
        if (i == 1) {
            z_ind_vector[series] = 0.5 * y_ind_vector[series];
        } else {
            z_ind_vector[series] += x * w_ind_vector[series];
        }
        w_ind_vector[series] = x;
        if (!(i % b_sqrt)) {
            j_sqrt++;
            batch_updates(j_sqrt, r_sqrt, s, w_matrix[series],
                          s_matrix[series], theta_matrix[series],
                          xi_matrix[series], r_max);
            y_sqrt += w_matrix[series][r_sqrt - 1] *
                      w_matrix[series][r_sqrt - 1];
            //changed from r_sqrt to r_sqrt - 1 b/c c array starts from zero
            //passing row vector(column vector in fortran)
        }
        if (!(i % b)) {
            j++;
            batch_updates(j, r, s, w_matrix[series],
                          s_matrix[series], theta_matrix[series],
                          xi_matrix[series], r_max);
            y += w_matrix[series][r - 1] * w_matrix[series][r - 1];
            //changed from r to r - 1 b/c c array starts from zero
            if (j == l) {
                //all row's are changed to row - 1
                old_row = row_vector[series];
                row_vector[series]++;
                row = row_vector[series];
                b_matrix[series][row - 1] = b;
                l_matrix[series][row - 1] = l;
                x_bar_matrix[series][row - 1] = s / ((double)(b * l));
                tau = (s * s) / ((double)l);
                rho = y - tau;
                v_matrix[series][row - 1] = rho / ((double)b * (double)b *
                                               ((double)(l - 1)));
                if (rho > 0.0) {
                    c = (-tau + theta_matrix[series][r - 1] +
                         xi_matrix[series][r - 1] + 0.5 *
                         w_matrix[series][r - 1] * w_matrix[series][r - 1]) /
                         rho;
                    p_matrix[series][row - 1] = 1.0 - phi_func(c /
                    sqrt((double)(l - 2) / ((double)l * (double)l - 1.0)));
                } else {
                    p_matrix[series][row - 1] = 0.0;
                }
                if (j % 2) {
                    y -= w_matrix[series][r - 1] * w_matrix[series][r - 1];
                }
                y += 2.0 * xi_matrix[series][r - 1];
                j >>= 1;
                r += 2;
                if (test == 1 &&
                    p_matrix[series][row - 1] <= beta_vector[series]) {
                    if (b > 1) {
                        b_sqrt <<= 1;
                        if (j_sqrt % 2) {
                            y_sqrt -= w_matrix[series][r_sqrt - 1] *
                                      w_matrix[series][r_sqrt - 1];
                        }
                        y_sqrt += 2.0 * xi_matrix[series][r_sqrt - 1];
                        j_sqrt >>= 1;
                        r_sqrt += 2;
                    }
                    b <<= 1;
                } else {
                    if (b == 1) {
                        b <<= 1;
                    } else {
                        swap_int = b << 1;
                        b = b_sqrt;
                        b_sqrt = swap_int;
                        swap_int = l << 1;
                        l = l_sqrt;
                        l_sqrt = swap_int;
                        swap_int = j;
                        j = j_sqrt;
                        j_sqrt = swap_int;
                        swap_dbl = y;
                        y = y_sqrt;
                        y_sqrt = swap_dbl;
                        swap_int = r;
                        r = r_sqrt;
                        r_sqrt = swap_int;
                    }
                    test = method_vector[series];
                }
            }
        }
        if (i > 2) {
            //puts("i > 2");
            tau = (ind_sum[series] * ind_sum[series]) / (double)i;
            rho = y_ind_vector[series] - tau;
            //changed from row_vector[series] + 1
            v_matrix[series][row_vector[series]] = rho / ((double)(i - 1));
            if (rho > 0) {
                c = (-tau + z_ind_vector[series] + 0.5 *
                     w_ind_vector[series] *
                     w_ind_vector[series]) / rho;
                //changed from row_vector[series] + 1
                p_matrix[series][row_vector[series]] = 1.0 -
                phi_func(c /
                sqrt((double)(i - 2) / ((double)i * (double)i - 1.0)));
            } else {
                //changed from row_vector[series] + 1
                p_matrix[series][row_vector[series]] = 0.0;
            }
        }
        if (i < t) {
            b_vector[series] = b;
            b_sqrt_vector[series] = b_sqrt;
            j_vector[series] = j;
            j_sqrt_vector[series] = j_sqrt;
            l_vector[series] = l;
            l_sqrt_vector[series] = l_sqrt;
            r_vector[series] = r;
            r_sqrt_vector[series] = r_sqrt;
            s_vector[series] = s;
            test_vector[series] = test;
            y_vector[series] = y;
            y_sqrt_vector[series] = y_sqrt;
        }
    }
    if (i == t) {
        //puts("writing to file");
        int stat = write_tableau_file(delta, head, s_num, i, t_prime,
                                      beta_vector, ind_sum, method_vector,
                                      row_vector, b_matrix,
                                      l_matrix, p_matrix,
                                      v_matrix,
                                      x_bar_matrix, fout);
        if (stat) {
            return IO_ERROR;
        }
    }
    i++;
    } while (fin && i <= t); 
    return NORMAL;
}
