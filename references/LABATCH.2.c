

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* copy.c contains interactive version that lists X_BAR, C.V.(X_BAR), and     */
/*        Sqrt[L*W(L,B)] on screen for each series at each iteration.         */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                                                                            */
/*         L        A    BBBB     A    TTTTT   CCC   H   H        222         */
/*         L       A A   B   B   A A     T    C   C  H   H       2   2        */
/*         L      A   A  B   B  A   A    T    C      H   H          2         */
/*         L      AAAAA  BBBB   AAAAA    T    C      HHHHH         2          */
/*         L      A   A  B   B  A   A    T    C      H   H   *    2           */
/*         L      A   A  B   B  A   A    T    C   C  H   H  * *  2            */
/*         LLLLL  A   A  BBBB   A   A    T     CCC   H   H   *   22222        */
/*                                                                            */
/*                                                                            */
/*                     (Revision of LABATCH Version 1.0)                      */
/*                                                                            */
/*                             C Implementation                               */
/*                                                                            */
/*                               G. S. Fishman                                */
/*                                                                            */
/*                        Translated by C. Arguelles                          */
/*                                 from the                                   */
/*                              FORTRAN Version                               */
/*                                                                            */
/*                               October 1997                                 */
/*                                                                            */
/*                                                                            */
/*                      Department of Operations Research                     */
/*                        210 Smith Building  CB # 3180                       */
/*                 University of North Carolina at Chapel Hill                */
/*                         Chapel Hill, NC 27599-3180                         */
/*                               (919) 962-8401                               */
/*                        email: gfish@fish.or.unc.edu                        */
/*                          fax: (919) 962-0391                               */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* LABATCH.2 is described in:                                                 */
/*                                                                            */
/* Fishman, G.S. (1997). LABATCH.2: An enhanced implementation of the batch   */
/* means method, Technical Report No. 97/04,Operations Research Department,   */
/* University of North Carolina at Chapel Hill.                               */
/*                                                                            */
/* The original LABATCH analysis package is based on:                         */
/*                                                                            */
/* Yarberry, L.S. (1993). Incorporating a dynamic batch size selection        */
/* mechanism in a fixed-sample-size batch means procedure. Ph.D. Dissertation,*/
/* Department of Operations Research, University of North Carolina at         */
/* Chapel Hill.                                                               */
/*                                                                            */
/* Use of Version 1.0 is illustrated in:                                      */
/*                                                                            */
/* Fishman, G.S. (1996). Monte Carlo: Concepts, Algorithms, and Applications, */
/* Springer-Verlag, New York.                                                 */
/*                                                                            */
/* and                                                                        */
/*                                                                            */
/* Fishman, G.S. and L.S. Yarberry (1997). An implementation of the batch     */
/* means method, INFORMS Journal on Computing, 9,296-310.                     */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                                                                            */
/*   *******************************************************************      */
/*   *                                                                 *      */
/*   * To reduce the potential for error introduced by unauthorized    *      */
/*   *modifications,it is recommended that this source code be obtained*      */
/*   *from http://www.or.unc.edu/~gfish/labatch.2.html                 *      */
/*   *                                                                 *      */
/*   *******************************************************************      */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*   Batch Means                                                              */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*                                                                            */
enum {S_MAX = 60, V_MAX = 60};
int BATCH_MEANS(IN_UNIT, OUT_UNIT, T, S_NUM, PSI_VECTOR, DELTA, RULE,
                BETA, L_UPPER, SCREEN)
int    IN_UNIT;
int    OUT_UNIT;
int    T;
int    S_NUM;
double PSI_VECTOR[];
double DELTA;
int    RULE;
double BETA;
int    L_UPPER;
int    SCREEN;
{
     extern FILE *in_file, *out_file;
     int i__1;
     double d__1, d__2, d__3, d__4;

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                                                                            */
/*    Input Parameters (All are local to subroutine BATCH_MEANS.)             */
/*                                                                            */
/*       IN_UNIT      := 0  if data are to be repeatedly                      */
/*                          transferred from a main program                   */
/*                                                                            */
/*                    := 30 if data are to be read one vector                 */
/*                          at a time from input file fort.30                 */
/*                          (n.b., 30 is merely an example)                   */
/*                                                                            */
/*       OUT_UNIT     := designated unit for writing output                   */
/*                                                                            */
/*       T            := total number of observations                         */
/*                                                                            */
/*       S_NUM        := number of sample series                              */
/*                                                                            */
/*       PSI_VECTOR   := vector of S_NUM sample values                        */
/*                                                                            */
/*       DELTA        := desired significance level for                       */
/*                       confidence intervals  (.01 or .05 suggested)         */
/*                                                                            */
/*       RULE         := 1  if ABATCH Rule is used                            */
/*                    := 0  if LBATCH Rule is used                            */
/*                                                                            */
/*       BETA         := significance level for testing for independence      */
/*                       ( .10 suggested)                                     */
/*                    N.B. If BETA =  1, each review uses the FNB rule.        */
/*                                 = -1, each review uses the SQRT rule.       */
/*                                                                            */
/*       L_UPPER      := upper bound on initial number of batches             */
/*                       (3 <= L_UPPER <= 100)                                */
/*                                                                            */
/*       SCREEN       := 1 if interim review estimates are to be              */
/*                         displayed on the screen.                           */
/*                    := 0 otherwise.                                         */
/*                                                                            */
/*                                                                            */
/*    Array Dimensions                                                        */
/*                                                                            */
/*       S_MAX        := maximal allowable number of sample series            */
/*       V_MAX        := maximal allowable number of reviews  + 1             */
/*                       (one for the independent case)                       */
/*                                                                            */
/*                                                                            */
/*    Important Parameters and Variables                                      */
/*                                                                            */
/*       B_1              := initial batch size                               */
/*       B_2_SQRT         := minimal permissible batch size using             */
/*                           SQRT Rule                                        */
/*                                                                            */
/*       L_1              := initial number of batches                        */
/*       L_2_SQRT         := minimal permissible number of batches            */
/*                           using SQRT Rule                                  */
/*                                                                            */
/*       HEAD             := 1 to produce headings in the computed            */
/*                             tableau file (default is 1 in program)         */
/*                        := 0 to omit headings                               */
/*                                                                            */
/*       T_PRIME          := observations 1,...,T_PRIME are used              */
/*                           to compute B*W(L,B) on final review              */
/*                                                                            */
/*       Subroutine SIZE chooses (L_1,B_1) first to maximize T_PRIME          */
/*                       and then to maximize the number of reviews.          */
/*                                                                            */
/*----------------------------------------------------------------------------*/

     static int                B_1;
     static int                B_2_SQRT;
     static int                HEAD;
     static int                L_1;
     static int                L_2_SQRT;
     static int                T_PRIME;

     static double             BETA_VECTOR[S_MAX];
     static int                METHOD_VECTOR[S_MAX];

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Working Storage                                                         */
/*                                                                            */
/*    In the following definitions, the subscripts I and J are used           */
/*    to denote the value of a quantity for series I on review J.             */
/*                                                                            */
/*       ROW_VECTOR(I)    := number of completed reviews                      */
/*       B_MATRIX(I,J)    := batch size                                       */
/*       L_MATRIX(I,J)    := number of batches                                */
/*       P_MATRIX(I,J)    := p-value                                          */
/*       V_MATRIX(I,J)    := sample variance of a batch mean                  */
/*                           (W(L,B) in the output tableaus)                  */
/*       X_BAR_MATRIX(I,J):= sample mean                                      */
/*                                                                            */
/*    Statistics for the independent case of series I are maintained          */
/*    in column # ROW_VECTOR(I)+1 of B_MATRIX, L_MATRIX, P_MATRIX,            */
/*    V_MATRIX and X_BAR_MATRIX.                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

     static int                ROW_VECTOR[S_MAX];
     static int                B_MATRIX[S_MAX][V_MAX + 1] = {{0}};
     static double             IND_SUM[S_MAX] = {0.};
     static int                L_MATRIX[S_MAX][V_MAX + 1] = {{0}};
     static double             P_MATRIX[S_MAX][V_MAX + 1] = {{0.}};
     static double             REL[S_MAX] = {{0.}};
     static double             V_MATRIX[S_MAX][V_MAX + 1] = {{0.}};
     static double             X_BAR_MATRIX[S_MAX][V_MAX + 1] = {{0.}};

     int                       B;
     int                       B_SQRT;
     static int                I;
     int                       II;
     int                       J;
     int                       J_SQRT;
     static int                JJ;
     char                      KK;
     int                       L;
     int                       L_SQRT;
     static int                OLD_ROW;
     int                       R;
     static int                R_MAX;
     int                       ROW = 0;
     int                       R_SQRT;
     int                       SERIES;
     int                       SWAP_INT;
     static int                SS = 0;
     int                       TEST;
     double                    C; 
     double                    RHO; 
     double                    S; 
     double                    SWAP_DBL; 
     double                    TAU; 
     double                    X; 
     double                    Y; 
     double                    Y_SQRT; 
     static int                B_VECTOR[S_MAX] = {0};
     static int                B_SQRT_VECTOR[S_MAX] = {0};
     static int                J_VECTOR[S_MAX] = {0};
     static int                J_SQRT_VECTOR[S_MAX] = {0};
     static int                L_VECTOR[S_MAX] = {0};
     static int                L_SQRT_VECTOR[S_MAX] = {0};
     static int                R_SQRT_VECTOR[S_MAX] = {0};
     static int                R_VECTOR[S_MAX] = {0};
     static int                TEST_VECTOR[S_MAX] = {0};
     static double             S_VECTOR[S_MAX] = {0.};
     static double             W_IND_VECTOR[S_MAX] = {0.};
     static double             Y_IND_VECTOR[S_MAX] = {0.};
     static double             Y_SQRT_VECTOR[S_MAX] = {0.};
     static double             Y_VECTOR[S_MAX] = {0.};
     static double             Z_IND_VECTOR[S_MAX] = {0.};
     static double             *S_MATRIX[S_MAX];
     static double             *THETA_MATRIX[S_MAX];
     static double             *W_MATRIX[S_MAX];
     static double             *XI_MATRIX[S_MAX];

     extern double             PHI();
     extern int                SIZE();
     extern int                BATCH_UPDATES();
     extern int                WRITE_TABLEAU_FILE();

     if (SS == 0) {
        for (i__1 = 0; i__1 < S_MAX; i__1++) {
           S_MATRIX[i__1]     = (double *) calloc(V_MAX + 1, sizeof(double));
           THETA_MATRIX[i__1] = (double *) calloc(V_MAX + 1, sizeof(double));
           W_MATRIX[i__1]     = (double *) calloc(V_MAX + 1, sizeof(double));
           XI_MATRIX[i__1]    = (double *) calloc(V_MAX + 1, sizeof(double));
        }

/*----------------------------------------------------------------------------*/

        HEAD = 1;
        for (SERIES = 0; SERIES < S_NUM; SERIES++) {
           METHOD_VECTOR[SERIES] = RULE;
           BETA_VECTOR[SERIES]   = BETA;
        }
        SIZE(T, L_UPPER, &B_1, &B_2_SQRT, &L_1, &L_2_SQRT, &T_PRIME);
     }
L2:
     if (IN_UNIT > 0) {
        for (SERIES = 0; SERIES < S_NUM; SERIES++) {
           fscanf(in_file, "%lf", &PSI_VECTOR[SERIES]);
        }
     }
     for (SERIES = 0; SERIES < S_NUM; SERIES++) {
        if (SS < S_NUM) {
/*----------------------------------------------------------------------------*/
/*    Begin initialization.                                                   */
/*----------------------------------------------------------------------------*/
           if (SS == 0) {
              if (S_NUM > S_MAX) {
                 printf("Maximum number of series exceeded.\n");
                 exit(0);
              }
              I = 1;
           }
           SS++;
           B      = B_1;
           B_SQRT = B_2_SQRT;
           L      = L_1;
           L_SQRT = L_2_SQRT;
           J      = 0;
           Y      = 0.;
           R      = 1;
           J_SQRT = 0;
           Y_SQRT = 0.;
           R_SQRT = 2;
             d__3 = log((double) (T / B)) / log(2.);
             d__4 = log((double) (T / B_SQRT)) / log(2.);
             d__1 = ((int) d__3) * 2 + 1;
             d__2 = ((int) d__4) * 2 + 2;
           R_MAX  = (d__1 > d__2) ? d__1 : d__2;
           if (R_MAX + 1 > V_MAX) {
              printf("Maximum vector size exceeded.\n");
              exit(0);
           }
           TEST = 1;
           S    = 0.;
           ROW_VECTOR[SERIES] = 0;
/*----------------------------------------------------------------------------*/
/*    End initialization.                                                     */
/*----------------------------------------------------------------------------*/
        } else {
/*----------------------------------------------------------------------------*/
/*    Begin restoration.                                                      */
/*----------------------------------------------------------------------------*/
           B      = B_VECTOR[SERIES];
           L      = L_VECTOR[SERIES];
           S      = S_VECTOR[SERIES];
           B_SQRT = B_SQRT_VECTOR[SERIES];
           J      = J_VECTOR[SERIES];
           J_SQRT = J_SQRT_VECTOR[SERIES];
           L_SQRT = L_SQRT_VECTOR[SERIES];
           R      = R_VECTOR[SERIES];
           R_SQRT = R_SQRT_VECTOR[SERIES];
           TEST   = TEST_VECTOR[SERIES];
           Y      = Y_VECTOR[SERIES];
           Y_SQRT = Y_SQRT_VECTOR[SERIES];
/*----------------------------------------------------------------------------*/
/*    End restoration.                                                        */
/*----------------------------------------------------------------------------*/
        }
        X = PSI_VECTOR[SERIES];
        if (I <= T_PRIME) {
           S = S + X;
        }
        IND_SUM[SERIES] += X;
/*----------------------------------------------------------------------------*/
/*    Update summary data for independent case.                               */
/*----------------------------------------------------------------------------*/
        Y_IND_VECTOR[SERIES] += X * X;
        if (I == 1) {
           Z_IND_VECTOR[SERIES] = Y_IND_VECTOR[SERIES] * .5;
        } else {
           Z_IND_VECTOR[SERIES] += X * W_IND_VECTOR[SERIES];
        }
        W_IND_VECTOR[SERIES] = X;
        if (I % B_SQRT == 0) {
/*----------------------------------------------------------------------------*/
/*    If a batch of size B_SQRT is complete:                                  */
/*----------------------------------------------------------------------------*/
           J_SQRT++;
           BATCH_UPDATES(J_SQRT, R_SQRT, S, W_MATRIX[SERIES],
                         S_MATRIX[SERIES], THETA_MATRIX[SERIES],
                         XI_MATRIX[SERIES], R_MAX);
           Y_SQRT += W_MATRIX[SERIES][R_SQRT] *
                     W_MATRIX[SERIES][R_SQRT];
        }
        if (I % B == 0) {
/*----------------------------------------------------------------------------*/
/*    If a batch of size B is complete                                        */
/*----------------------------------------------------------------------------*/
           J++;
           BATCH_UPDATES(J, R, S, W_MATRIX[SERIES],
                         S_MATRIX[SERIES], THETA_MATRIX[SERIES],
                         XI_MATRIX[SERIES], R_MAX);
           Y += W_MATRIX[SERIES][R] * W_MATRIX[SERIES][R];
           if (J == L) {
/*----------------------------------------------------------------------------*/
/*    If the current review is over, compute statistics.                      */
/*----------------------------------------------------------------------------*/
              OLD_ROW                   = ROW_VECTOR[SERIES];
              ROW_VECTOR[SERIES]++;
              ROW                       = ROW_VECTOR[SERIES];
              B_MATRIX[SERIES][ROW]     = B;
              L_MATRIX[SERIES][ROW]     = L;
              X_BAR_MATRIX[SERIES][ROW] = S / (double) (B * L);
              TAU                       = (S * S) / (double) L;
              RHO                       = Y - TAU;
              V_MATRIX[SERIES][ROW]     = RHO / ((double) B *
                               (double) B * (double) (L - 1));
              if (RHO > 0.) {
                 C = (-TAU + THETA_MATRIX[SERIES][R] +
                     XI_MATRIX[SERIES][R] + W_MATRIX[SERIES][R] * .5 *
                     W_MATRIX[SERIES][R]) / RHO;
                 d__1 = C / sqrt((double) (L - 2) / ((double) L *
                        (double) L - 1.));
                 P_MATRIX[SERIES][ROW] = 1. - PHI(d__1);
              } else {
                 P_MATRIX[SERIES][ROW] = 0.; 
              }
/*----------------------------------------------------------------------------*/
/*    End compute statistics.                                                 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*    If J is odd:                                                            */
/*----------------------------------------------------------------------------*/
              if (J % 2 == 1) {
                 Y -= W_MATRIX[SERIES][R] * W_MATRIX[SERIES][R];
              }
              Y += XI_MATRIX[SERIES][R] * 2.;
              J /= 2;
              R += 2;
              if (TEST == 1 && P_MATRIX[SERIES][ROW] <= 
                  BETA_VECTOR[SERIES]) {
                 if (B > 1) {
                    B_SQRT *= 2;
/*----------------------------------------------------------------------------*/
/*    If J_SQRT is odd:                                                       */
/*----------------------------------------------------------------------------*/
                    if (J_SQRT % 2 == 1) {
                       Y_SQRT -= W_MATRIX[SERIES][R_SQRT] *
                                 W_MATRIX[SERIES][R_SQRT];
                    }
                    Y_SQRT += XI_MATRIX[SERIES][R_SQRT] * 2.;
                    J_SQRT /= 2;
                    R_SQRT += 2;
                 }
/*----------------------------------------------------------------------------*/
/*    FNB Rule                                                                */
/*----------------------------------------------------------------------------*/
                 B *= 2;
              } else {
                 if (B == 1) {
/*----------------------------------------------------------------------------*/
/*    FNB Rule                                                                */
/*----------------------------------------------------------------------------*/
                    B *= 2;
                 } else {
                    SWAP_INT = 2 * B;
/*----------------------------------------------------------------------------*/
/*    SQRT Rule                                                               */
/*----------------------------------------------------------------------------*/
                    B        = B_SQRT;
                    B_SQRT   = SWAP_INT;
                    SWAP_INT = 2 * L;
/*----------------------------------------------------------------------------*/
/*    SQRT Rule                                                               */
/*----------------------------------------------------------------------------*/
                    L        = L_SQRT;
                    L_SQRT   = SWAP_INT;
                    SWAP_INT = J;
                    J        = J_SQRT;
                    J_SQRT   = SWAP_INT;
                    SWAP_DBL = Y;
                    Y        = Y_SQRT;
                    Y_SQRT   = SWAP_DBL;
                    SWAP_INT = R;
                    R        = R_SQRT;
                    R_SQRT   = SWAP_INT;
                 }
                 TEST = METHOD_VECTOR[SERIES];
              }
/*----------------------------------------------------------------------------*/
/*    Compute statistics for independent case.                                */
/*----------------------------------------------------------------------------*/
              TAU =  (IND_SUM[SERIES] * IND_SUM[SERIES]) / (double) I;
              RHO = Y_IND_VECTOR[SERIES] - TAU;
              V_MATRIX[SERIES][ROW_VECTOR[SERIES] + 1] = RHO / (double) (I - 1);
              if (RHO > 0.) {
                 C = (-TAU + Z_IND_VECTOR[SERIES] + W_IND_VECTOR[SERIES] * .5 *
                     W_IND_VECTOR[SERIES]) / RHO;
                 d__1 = C / sqrt((double) (I - 2) / ((double) I *
                        (double) I - 1.));
                 P_MATRIX[SERIES][ROW_VECTOR[SERIES] + 1] = 1. - PHI(d__1);
              } else {
                 P_MATRIX[SERIES][ROW_VECTOR[SERIES] + 1] = 0.;
              }
           }
        }
/*----------------------------------------------------------------------------*/
/*    Begin storage.                                                          */
/*----------------------------------------------------------------------------*/
        if (I < T) {                                        
           B_VECTOR[SERIES]      = B;
           B_SQRT_VECTOR[SERIES] = B_SQRT;
           J_VECTOR[SERIES]      = J;
           J_SQRT_VECTOR[SERIES] = J_SQRT;
           L_VECTOR[SERIES]      = L;
           L_SQRT_VECTOR[SERIES] = L_SQRT;
           R_VECTOR[SERIES]      = R;
           R_SQRT_VECTOR[SERIES] = R_SQRT;
           S_VECTOR[SERIES]      = S;
           TEST_VECTOR[SERIES]   = TEST;
           Y_VECTOR[SERIES]      = Y;
           Y_SQRT_VECTOR[SERIES] = Y_SQRT;
        }
/*----------------------------------------------------------------------------*/
/*    End storage.                                                            */
/*----------------------------------------------------------------------------*/
L10:
        ;
     }
        if (I == T) {
                    WRITE_TABLEAU_FILE(DELTA, HEAD, S_NUM, T,
                                T_PRIME, BETA_VECTOR, IND_SUM, 
                                METHOD_VECTOR, ROW_VECTOR, B_MATRIX, 
                                L_MATRIX, P_MATRIX, V_MATRIX, 
                                X_BAR_MATRIX, OUT_UNIT);
        }             
/*---------------------------------------------------------------------------*/
/*   Display selected interim results on screen.                             */
/*---------------------------------------------------------------------------*/
if (SCREEN == 1) {
   if (ROW > OLD_ROW ) { 
      if (OLD_ROW < 1) { 
         printf("\n\n\n");
         printf("LABATCH.2 INTERIM REVIEW STATISTICAL ANALYSIS\n\n\n");
         JJ = 6;
         if (S_NUM < 6) {      
            JJ = S_NUM; }       
         printf("X_BAR, C.V.(X_BAR), and Sqrt[B*W(L,B)]");
         printf(" for Series 1 Through %2d\n", JJ);
         printf("**************************************");
         printf("************************ \n\n\n");
         printf("         No. of  \n");
         printf("Review    Obs.                 ");
         for (II = 1; II<= JJ; II++) {       
            printf("     %1d     ", II); }  
         printf("\n\n\n");
      }
      printf("%3d %10d", ROW, I);     
      printf(" X_BAR         ");
      for (SERIES = 0; SERIES < JJ ; SERIES++) {                       
         d__1 = X_BAR_MATRIX[SERIES][ROW];
         printf("%12.3e", d__1);}                                   
      printf("\n               C.V.(X_BAR)   ");
      for (SERIES = 0; SERIES < JJ ; SERIES++) {                       
         REL[SERIES] = 0.;
         if (X_BAR_MATRIX[SERIES][ROW] != 0.) {
            REL[SERIES] = sqrt(V_MATRIX[SERIES][ROW]/
                              L_MATRIX[SERIES][ROW])/
                              fabs(X_BAR_MATRIX[SERIES][ROW]);
         }}
      for (SERIES = 0; SERIES < JJ ; SERIES++) {                       
         d__1 = REL[SERIES];
         printf("%12.3e", d__1);}                                   
      printf("\n               Sqrt[B*W(L,B)]");
      for (SERIES = 0; SERIES < JJ ; SERIES++) {                       
         d__1 =  sqrt(B_MATRIX[SERIES][ROW] * V_MATRIX[SERIES][ROW]); 
         printf("%12.3e", d__1);}                                   
      printf("\n");
L15:  printf("continue[y/n]? ");
L16:  if ((KK = getchar()) == '\n') {
         goto L15; } 
      else {
/*----------------------------------------------------------------------*/
/*    If continue = n, write tableau file and terminate execution.      */
/*----------------------------------------------------------------------*/
         if (KK == 'n') {
            WRITE_TABLEAU_FILE(DELTA, HEAD, S_NUM, I,
                            T_PRIME, BETA_VECTOR, IND_SUM, 
                            METHOD_VECTOR, ROW_VECTOR, B_MATRIX, 
                            L_MATRIX, P_MATRIX, V_MATRIX, 
                            X_BAR_MATRIX, OUT_UNIT);
            exit(0);   
         }
         else {
            if (KK != 'y') {
               goto L16;   
            }
            else {
               KK = getchar(); 
            }
         }
       }
/*----------------------------------------------------------------------*/
/*    End screen display for this review.                               */
/*----------------------------------------------------------------------*/
      OLD_ROW = ROW;
   } 
}     
     I++; 
     if (IN_UNIT > 0) {
        if (I <= T) {
           goto L2;
        }
     }
     return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*    Batch Updates                                                           */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int BATCH_UPDATES(J_DUMMY, R_DUMMY, S, W_VECTOR, 
                  S_VECTOR, THETA_VECTOR, XI_VECTOR, R_MAX)
int    J_DUMMY;
int    R_DUMMY;
double S;
double W_VECTOR[];
double S_VECTOR[];
double THETA_VECTOR[];
double XI_VECTOR[];
int    R_MAX;
{
     int                       J;
     int                       R;
     double                    W;

     J = J_DUMMY;
     R = R_DUMMY;
/*----------------------------------------------------------------------------*/
/*    While J is even:                                                        */
/*----------------------------------------------------------------------------*/
     while (J % 2 == 0) {
        W             = S - S_VECTOR[R];
        XI_VECTOR[R] += W * W_VECTOR[R];
        W_VECTOR[R]   = W;
        S_VECTOR[R]   = S;
        J            /= 2;
        R            += 2;
     }
     if (J == 1) {
        W               = S;
        THETA_VECTOR[R] = W * .5 * W;
        XI_VECTOR[R]    = 0.;
     } else {
        W                = S - S_VECTOR[R];
        THETA_VECTOR[R] += W * W_VECTOR[R];
     }
     W_VECTOR[R] = W;
     S_VECTOR[R] = S;
     return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*    Write tableau file.                                                     */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int WRITE_TABLEAU_FILE(DELTA, HEAD, S_NUM, T, T_PRIME, 
        BETA_VECTOR, IND_SUM, METHOD_VECTOR, ROW_VECTOR, B_MATRIX, 
        L_MATRIX, P_MATRIX, V_MATRIX, X_BAR_MATRIX, OUT_UNIT)
double DELTA;
int    HEAD;
int    S_NUM;
int    T;
int    T_PRIME;

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Array Dimensions                                                        */
/*                                                                            */
/*       S_MAX       := the maximum number of sample series                   */
/*       V_MAX       := the maximum number of reviews + 1                     */
/*                           (one for the independent case)                   */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    Input Parameters                                                        */
/*                                                                            */
/*       DELTA            := desired significance level for                   */
/*                           confidence intervals                             */
/*       HEAD             := 1 to produce headings in the computed            */
/*                           tableau file                                     */
/*                        := 0 to omit headings                               */
/*       S_NUM            := number of sample series                          */
/*       T                := number oof observations                          */
/*       T_PRIME          := number of observations used for variance         */
/*                           analysis                                         */
/*       BETA_VECTOR      := vector of significance levels for                */
/*                           testing independence for each Sample             */
/*       IND_SUM          := sum of the values of the observations            */
/*       METHOD_VECTOR    := vector determining which rule to use             */
/*                           (set to 1 for ABATCH, 0 for LBATCH)              */
/*       OUT_UNIT         := designated unit for writing output               */
/*                                                                            */
/*    Working Storage                                                         */
/*                                                                            */
/*    In the following definitions, the subscripts I and J are used           */
/*    to denote the value of a quantity for series I on review J.             */
/*                                                                            */
/*       ROW_VECTOR(I)    := number of completed reviews                      */
/*       B_MATRIX(I,J)    := batch size                                       */
/*       L_MATRIX(I,J)    := number of batches                                */
/*       P_MATRIX(I,J)    := p-value                                          */
/*       V_MATRIX(I,J)    := sample variance of a batch mean                  */
/*                           (W(L,B) in the output tableaus)                  */
/*       X_BAR_MATRIX(I,J):= sample mean                                      */
/*                                                                            */
/*    statistics for the independent case of series I are maintained          */
/*    in column # ROW_VECTOR(I)+1 of B_MATRIX, L_MATRIX, P_MATRIX,            */
/*    V_MATRIX and X_BAR_MATRIX.                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

double BETA_VECTOR[];
double IND_SUM[];
int    METHOD_VECTOR[];
int    ROW_VECTOR[];
int    B_MATRIX[][V_MAX + 1];
int    L_MATRIX[][V_MAX + 1];
double P_MATRIX[][V_MAX + 1];
double V_MATRIX[][V_MAX + 1];
double X_BAR_MATRIX[][V_MAX + 1];
int    OUT_UNIT;
{

     extern FILE *out_file;
     int         i__1;
     double      d__1, d__2;

     int                       B;
     int                       L;
     int                       LINES;
     double                    REL;
     int                       ROW;
     int                       S; 
     int                       SERIES;
     double                    H;
     double                    P;
     double                    V;
     double                    X_BAR;

     extern double             STUDENT_T_QUANTILE();

     if (HEAD == 1) {
        fprintf(out_file, "\n\n");
        fprintf(out_file, "*---------------------------------------");
        fprintf(out_file, "---------------------------------------*\n");
        fprintf(out_file, "|---------------------------------------");
        fprintf(out_file, "---------------------------------------|\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|          L        A    BBBB     A    T");
        fprintf(out_file, "TTTT   CCC   H   H        222          |\n");
        fprintf(out_file, "|          L       A A   B   B   A A    ");
        fprintf(out_file, " T    C   C  H   H       2   2         |\n");
        fprintf(out_file, "|          L      A   A  B   B  A   A   ");
        fprintf(out_file, " T    C      H   H          2          |\n");
        fprintf(out_file, "|          L      AAAAA  BBBB   AAAAA   ");
        fprintf(out_file, " T    C      HHHHH         2           |\n");
        fprintf(out_file, "|          L      A   A  B   B  A   A   ");
        fprintf(out_file, " T    C      H   H   *    2            |\n");
        fprintf(out_file, "|          L      A   A  B   B  A   A   ");
        fprintf(out_file, " T    C   C  H   H  * *  2             |\n");
        fprintf(out_file, "|          LLLLL  A   A  BBBB   A   A   ");
        fprintf(out_file, " T     CCC   H   H   *   22222         |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                      (Revision of LABA");
        fprintf(out_file, "TCH Version 1.0)                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                              C Impleme");
        fprintf(out_file, "ntation                                |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                G. S. F");
        fprintf(out_file, "ishman                                 |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                         Translated by ");
        fprintf(out_file, "C. Arguelles                           |\n");      
        fprintf(out_file, "|                                  from ");
        fprintf(out_file, "the                                    |\n");      
        fprintf(out_file, "|                               FORTRAN ");
        fprintf(out_file, "Version                                |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                October");
        fprintf(out_file, " 1997                                  |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                       Department of Op");
        fprintf(out_file, "erations Research                      |\n");
        fprintf(out_file, "|                         210 Smith Buil");
        fprintf(out_file, "ding  CB # 3180                        |\n");
        fprintf(out_file, "|                  University of North C");
        fprintf(out_file, "arolina at Chapel Hill                 |\n");
        fprintf(out_file, "|                          Chapel Hill, ");
        fprintf(out_file, "NC 27599-3180                          |\n");
        fprintf(out_file, "|                                (919) 9");
        fprintf(out_file, "62-8401                                |\n");
        fprintf(out_file, "|                         email: gfish@f");
        fprintf(out_file, "ish.or.unc.edu                         |\n");
        fprintf(out_file, "|                           fax: (919) 9");
        fprintf(out_file, "62-0391                                |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|---------------------------------------");
        fprintf(out_file, "---------------------------------------|\n");
        fprintf(out_file, "|---------------------------------------");
        fprintf(out_file, "---------------------------------------|\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "| LABATCH.2 is described in:            ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "| Fishman, G.S. (1997). LABATCH.2: An en");
        fprintf(out_file, "hanced implementation of the batch     |\n");
        fprintf(out_file, "| means method, Technical Report No. 97/");
        fprintf(out_file, "04, Operations Research Department,    |\n");
        fprintf(out_file, "| University of North Carolina at Chapel");  
        fprintf(out_file, " at Chapel Hill.                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");  
        fprintf(out_file, "| The original LABATCH analysis package ");
        fprintf(out_file, "is based on:                           |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "| Yarberry, L.S. (1993). Incorporating a");
        fprintf(out_file, " dynamic batch size selection mechanism|\n");
        fprintf(out_file, "| in a fixed-sample-size batch means pro");
        fprintf(out_file, "cedure. Ph.D. Dissertation, Department |\n");
        fprintf(out_file, "| of Operations Research, University of ");
        fprintf(out_file, "North Carolina at Chapel Hill.         |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "| Use of Version 1.0 is illustrated in: ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "| Fishman, G.S. (1996). Monte Carlo: Con");
        fprintf(out_file, "cepts, Algorithms, and Applications,   |\n");
        fprintf(out_file, "| Springer-Verlag, New York.            ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "| and                                   ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "| Fishman, G.S. and L.S. Yarberry (1997)");
        fprintf(out_file, ". An implementation of the batch means |\n");
        fprintf(out_file, "| method, INFORMS Journal on Computing, ");
        fprintf(out_file, "9, 296-310.                            |\n");
        fprintf(out_file, "|                                       ");
        fprintf(out_file, "                                       |\n");
        fprintf(out_file, "|---------------------------------------");
        fprintf(out_file, "---------------------------------------|\n");
        fprintf(out_file, "*---------------------------------------");
        fprintf(out_file, "---------------------------------------*\n\n");
     }

     fprintf(out_file, "\n\n");
     fprintf(out_file, "Final Tableau                            \n");
     fprintf(out_file, "                                    Mean ");
     fprintf(out_file, "Estimation                               \n");
     fprintf(out_file, "                                    *****");
     fprintf(out_file, "**********                               \n");
     fprintf(out_file, "                                    (t =%9d", T);
     fprintf(out_file, " )                                       \n\n");
     fprintf(out_file, "                                         ");
     fprintf(out_file, "     %4.1f%%                             \n",
            (1. - DELTA)*100.);
     fprintf(out_file, "            _       Standard Error      C");
     fprintf(out_file, "onfidence Interval                    _  \n");
     fprintf(out_file, "Series      X      Sqrt[B*W(L,B)/t]     ");
     fprintf(out_file, "Lower        Upper      (Upper-Lower)/|X|\n\n");

     for (SERIES = 0; SERIES < S_NUM; SERIES++) {
        ROW   = ROW_VECTOR[SERIES];
        B     = B_MATRIX[SERIES][ROW];
        L     = L_MATRIX[SERIES][ROW];
        S     = L * B;
        if (S == T_PRIME) {
           S = T; }

        X_BAR = IND_SUM[SERIES] / (double) S;
        V     = V_MATRIX[SERIES][ROW];
        d__1  = 1. - DELTA / 2.;
        i__1  = L - 1;
        H     = STUDENT_T_QUANTILE(d__1, i__1) * sqrt(B * V / S);
        REL   = 0.;
        if (X_BAR != 0.) {
           REL = 2*H/fabs(X_BAR); }
        fprintf(out_file, "%4d  %11.3e    ", SERIES + 1, X_BAR);
        fprintf(out_file, " %11.3e    %11.3e  ", sqrt(B * V / S), X_BAR - H);
        d__1 = REL;           
        fprintf(out_file, "%11.3e      %11.3e\n\n", X_BAR + H, d__1);
     }
     fprintf(out_file, "   _                                    ");
     fprintf(out_file, "                                        \n");
     fprintf(out_file, "   X(t) is based on all t observations. \n");
     fprintf(out_file, "   W(L,B) is based on first %6.2f", L * B * 100. / T);
     fprintf(out_file, "%% of the t observations.                \n");
     if (HEAD == 1) {
        for (LINES = 0; LINES < 66 - 13 - 2 * S_NUM; LINES++) {
           fprintf(out_file, "\n");
        }
     }
     for (SERIES = 0; SERIES < S_NUM; SERIES++) {
        if (HEAD == 1) {
           if (METHOD_VECTOR[SERIES] == 1) {
              fprintf(out_file, "Interim Review Tableau          \n");
              fprintf(out_file, "                       ");
              fprintf(out_file, "ABATCH Data Analysis for Series %d\n", SERIES+1);
           } else {
              fprintf(out_file, "Interim Review Tableau           \n");
              fprintf(out_file, "                       ");
              fprintf(out_file, "LBATCH Data Analysis for Series %d\n", SERIES+1);
           }
           fprintf(out_file, "\n");
           fprintf(out_file, "                                          ");
           fprintf(out_file, "           %4.1f%%                        \n",
                  (1. - DELTA)*100.);
           fprintf(out_file, "                                     _    ");
           fprintf(out_file, "   Confidence Interval                    \n");
           fprintf(out_file, "Review    L*B      L         B       X    ");
           fprintf(out_file, "    Lower       Upper Sqrt[B*W(L,B)]  p-va");
           fprintf(out_file, "lue                                       \n");
        }
        for (ROW = 1; ROW <= ROW_VECTOR[SERIES]; ROW++) {
           B     = B_MATRIX[SERIES][ROW];
           L     = L_MATRIX[SERIES][ROW];
           P     = P_MATRIX[SERIES][ROW];
           X_BAR = X_BAR_MATRIX[SERIES][ROW];
           S     = L*B;
           if (S == T_PRIME) {
                S   = T;
              X_BAR = IND_SUM[SERIES]/S;
           }
           V     = V_MATRIX[SERIES][ROW];
           d__1  = 1. - DELTA / 2.;
           i__1  = L - 1;
           H     = STUDENT_T_QUANTILE(d__1, i__1) * sqrt(B * V / S);
           fprintf(out_file, "%3d %8d %8d %8d ", ROW, B * L, L, B);
           fprintf(out_file, "%11.3e %11.3e %11.3e ", X_BAR, X_BAR-H, X_BAR+H);
           fprintf(out_file, "%11.3e   %6.4f\n", sqrt((double) B * V), P); 
        }
        if (HEAD == 1) {
           fprintf(out_file, "\n  If data are independent:\n\n");
        }
        L = L*B;
        B = 1;
        P = P_MATRIX[SERIES][ROW_VECTOR[SERIES] + 1];
        V = V_MATRIX[SERIES][ROW_VECTOR[SERIES] + 1];
    d__1  = 1. - DELTA / 2.;
    i__1  = S - 1;
        H = STUDENT_T_QUANTILE(d__1, i__1) * sqrt(B * V / S);
        fprintf(out_file, "    %8d %8d %8d ", T, T, B);        
        fprintf(out_file, "%11.3e %11.3e %11.3e ", X_BAR, X_BAR-H, X_BAR+H);
        fprintf(out_file, "%11.3e   %6.4f\n", sqrt((double) B * V), P); 
        if(BETA_VECTOR[SERIES] < 0)   {
           BETA_VECTOR[SERIES] = 0;
        }
        fprintf(out_file, "\n %5.2f", BETA_VECTOR[SERIES]);
        fprintf(out_file, " significance level for independence");
        fprintf(out_file, " testing.\n");
        fprintf(out_file, "  Review %3d", ROW - 1);
        fprintf(out_file, " used the first %6.2f", L * B * 100. / T);
        fprintf(out_file, "%% of the t observations for W(L,B).   \n");
        if (HEAD == 1) {
           for (LINES = 0; LINES < 66 - 14 - ROW_VECTOR[SERIES]; LINES++) {
              fprintf(out_file, "\n");
           }
        }
     }
     return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*   PHI function                                                             */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double PHI(X)
double X;
{
     double d__1;
/*----------------------------------------------------------------------------*/
/* This function computes the value of the normal c.d.f. at the               */
/* point X.                                                                   */
/*                                                                            */
/* Reference: Abramowitz, M. and I.A. Stegun (1964). Handbook of              */
/* Mathematical Functions with Formulas, Graphs and Mathematical              */
/* Tables, U.S. Government Printing Office, Washington, D.C.                  */
/*----------------------------------------------------------------------------*/

     double RC;
     double Y;

     if (X > 0.) {
        if (X > 45.) {
           RC = 1.;
        } else {
           Y = X;
           d__1 = (((((Y * 5.383e-6 + 4.88906e-5) * Y + 3.80036e-5) * Y + 
                  .0032776263) * Y + .0211410061) * Y + .049867347) * Y + 1.;
           d__1 *= d__1;
           d__1 *= d__1;
           d__1 *= d__1;
           RC = 1. - .5 / (d__1 * d__1);
        }
    } else {
        if (X < -45.) {
           RC = 0.;
        } else {
           Y = -X;
           d__1 = (((((Y * 5.383e-6 + 4.88906e-5) * Y + 3.80036e-5) * Y + 
                  .0032776263) * Y + .0211410061) * Y + .049867347) * Y + 1.;
           d__1 *= d__1;
           d__1 *= d__1;
           d__1 *= d__1;
           RC = .5 / (d__1 * d__1);
        }
    }
    return RC;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*   PHI quantile function                                                    */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double PHI_QUANTILE(BETA_DUMMY)
double BETA_DUMMY;
{
     double d__1;
/*----------------------------------------------------------------------------*/
/* This function computes the quantile of order BETA from the normal          */
/* distribution.                                                              */
/*                                                                            */
/* Reference: Abramowitz, M. and I.A. Stegun (1964). Handbook of              */
/* Mathematical Functions with Formulas, Graphs and Mathematical              */
/* Tables, U.S. Government Printing Office, Washington, D.C.                  */
/*----------------------------------------------------------------------------*/

     double                    BETA;
     double                    DELTA;
     double                    HIGH;
     double                    LOW;
     double                    MID;
     double                    RC;
     double                    W;

     extern double             PHI();

     BETA = BETA_DUMMY;
     if (BETA >= .5) {
        DELTA = BETA;
     } else {
        DELTA = 1. - BETA;
     }
     W = sqrt(-log((1. - DELTA) * (1. - DELTA)));
     LOW = W - ((W * .010328 + .802853) * W + 2.515517) / (((W * .001308 + 
           .189269) * W + 1.432788) * W + 1.);
     HIGH = LOW + .00045;
     LOW -= .00045;
     while (fabs(MID - (LOW + HIGH) / 2.) > .000000000000001) {
        MID = (LOW + HIGH) / 2.;
        if (PHI(MID) < DELTA) {
           LOW = MID;
        } else {
           HIGH = MID;
        }
     }
     if (BETA >= .5) {
        RC = MID;
     } else {
        RC = -MID;
     }
     return RC;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*   Student t quantile function                                              */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double STUDENT_T_QUANTILE(BETA,L)
double BETA;
int    L;
{
     double d__1, d__2;
/*----------------------------------------------------------------------------*/
/* This function computes the quantile of order BETA from the                 */
/* Student t distribution with L degrees of freedom.                          */
/*                                                                            */
/* Reference: Abramowitz, M. and I.A. Stegun (1964). Handbook of              */
/* Mathematical Functions with Formulas, Graphs and Mathematical              */
/* Tables, U.S. Government Printing Office, Washington, D.C.                  */
/*----------------------------------------------------------------------------*/

     double                    DELTA;
     double                    DFL;
     double                    HIGH;
     double                    LOW;
     double                    MID;
     double                    RC;
     double                    TEST;
     double                    W; 
     double                    X; 
     double                    Y;

     extern double             PHI_QUANTILE();

     if (BETA == .5) {
        RC = 0.;
     } else {
        if (L >= 7) {
           X = PHI_QUANTILE(BETA);
           Y = X * X;
           DFL = (double) L;
           d__1 = DFL;
           d__1 *= d__1;
           RC = X * ((((DFL * 23040. + 2880.) * DFL - 3600.) * DFL - 945. + 
                  (((DFL * 23040. + 15360.) * DFL + 4080.) * DFL - 1920. + ((
                  DFL * 4800. + 4560.) * DFL + 1482. + (DFL * 720. + 776. + Y * 
                  79.) * Y) * Y) * Y) / (d__1 * d__1 * 92160.) + 1.);
        } else {
           switch (L) {
           case 1:
              RC = tan((BETA - .5) * 3.1415926535898);
           case 2:
              RC = (BETA * 2. - 1.) / sqrt(BETA * 2. * (1. - BETA));
           default:
              if (BETA >= .5) {
                 DELTA = BETA;
              } else {
                 DELTA = 1. - BETA;
              }
              LOW = PHI_QUANTILE(DELTA);
              HIGH = (DELTA * 2. - 1.) / sqrt(DELTA * 2. * (1. - DELTA));
              while (fabs(MID - (LOW + HIGH) / 2.) > .000000000000001) {
                 MID = (LOW + HIGH) / 2.;
                 W = MID * MID;
                 switch (L) {
                 case 3:
                    TEST = sqrt(3.) * tan((DELTA - .5) * 3.1415926535898 -
                           MID * sqrt(3.) / (W + 3.));
                 case 4:
                    d__1 = sqrt(W + 4.);
                    TEST = (DELTA * 2. - 1.) * d__1 * d__1 * d__1 / (W + 6.);
                 case 5:
                    d__1 = W + 5.;
                    TEST = sqrt(5.) * tan((DELTA - .5) * 3.1415926535898 - MID
                           * sqrt(5.) / (d__1 * d__1 * 3.) * (W * 3. + 25.));
                 case 6:
                    d__1 = sqrt(W + 6.);
                    d__2 = d__1;
                    d__1 *= d__1;
                    TEST = (DELTA * 2. - 1.) * (d__2 * (d__1 * d__1)) / (W * W +
                           W * 15. + 67.5);
                 }
                 if (TEST > MID) {
                    LOW = MID;
                 } else {
                    HIGH = MID;
                 }
              }
              if (BETA >= .5) {
                 RC = MID;
              } else {
                 RC = -MID;
              }
           }
        }
     }
     return RC;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/* subroutine SIZE finds L_1 and B_1 first to maximize T_PRIME                */
/*           (<=100) and then to maximize the number of reviews.              */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int SIZE(T, L_UPPER, B_1, B_2_SQRT, L_1, L_2_SQRT, T_PRIME)
int T;
int L_UPPER;
int *B_1;
int *B_2_SQRT;
int *L_1;
int *L_2_SQRT;
int *T_PRIME;
{
     int               ALPHA;
     int               ALPHA_MAX;
     int               I;
     static int        I_MAX;
     int               TEMP_1;
     int               TEMP_2;

     static int B[] = {0,1,5,5,7,7,7,8,9,12,12,12,15,17,17,17,19,21,21,22,
         22,26,27,29,31,32,33,36,36,36,39,41,43,46,49,50,51,51,53,56,60,63,
         66,67,70,84,85};
     static int L[] = {0,3,7,21,15,25,35,11,13,17,51,85,21,36,60,84,27,25,
         35,31,93,37,57,41,66,45,47,51,68,85,55,87,61,65,69,71,60,84,75,79,
         85,89,93,95,99,85,96};

     I_MAX = 1;
     ALPHA_MAX = 0;
     TEMP_2 = 0;

     for (I = 1; I <= 46; I++) {
        if (T >= B[I] * L[I] && L[I] <= L_UPPER) {
           ALPHA = (int) (log((double) T / (L[I] * B[I])) / log(2.));
           TEMP_1 = pow(2., (double) ALPHA) * B[I] * L[I];
           if (TEMP_1 >= TEMP_2) {
              if (TEMP_1 > TEMP_2 || L[I] * B[I] < L[I_MAX] * B[I_MAX]) {
                 I_MAX = I;
                 ALPHA_MAX = ALPHA;
                 TEMP_2 = TEMP_1;
              }
           }
        }
     }

     *B_1 = B[I_MAX];
     *B_2_SQRT = 3;
     if (*B_1 > 1) {
        *B_2_SQRT = (int) (sqrt(2.) * *B_1 + .5);
     }

     *L_1 = L[I_MAX];
     *L_2_SQRT = (int) (sqrt(2.) * *L_1 + .5);
     *T_PRIME = TEMP_2;
     return 0;
}
