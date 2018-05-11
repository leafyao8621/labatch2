enum {
    S_MAX = 60,
    V_MAX = 60
};
typedef enum {
    NORMAL,
    MAX_NUM_SERIES,
    MAX_VECTOR_SIZE,
    IO_ERROR
} StatCode;
StatCode batch_means(FILE* fin, FILE* fout, int t, int s_num, double* psi_vector,
                 double delta, int rule, double beta, int l_upper, int screen);
