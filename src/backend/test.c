#include <stdio.h>
#include <stdlib.h>
#include "Lab2data.h"

int main(int argc, char** argv) {
    int t, s_num, rule, l_upper, screen;
    double beta, delta;
    if (argc < 4) {
        perror("Insufficient # of args!");
        return 1;
    }
    FILE* par = fopen(argv[2], "r");
    fscanf(par, "%d\n%d\n%lf\n%d\n%lf\n%d\n%d\n", &t, &s_num, &delta, &rule,
                                                  &beta, &l_upper, &screen);
    fclose(par);
    double data[s_num];
    FILE* fin = fopen(argv[1], "r");
    FILE* fout = fopen(argv[3], "w");
    int stat = batch_means(fin, fout, t, s_num, data, delta, rule, beta, l_upper, screen);
    fclose(fin);
    fclose(fout);
    //printf("%d\n", stat);
    return stat;
}
