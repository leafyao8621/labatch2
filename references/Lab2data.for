C     Last change:  CA   22 Mar 2002   10:29 am




C copy.for contains interactive version that lists X_BAR, C.V.(X_BAR), and 
C          Sqrt[VAR.EST.] on screen for each series at each iteration.


C*------------------------------------------------------------------------------*
C|------------------------------------------------------------------------------|
C|                                                                              |
C|                                                                              |
C|                                                                              |
C|          L        A    BBBB     A    TTTTT   CCC   H   H        222          |
C|          L       A A   B   B   A A     T    C   C  H   H       2   2         |
C|          L      A   A  B   B  A   A    T    C      H   H          2          |
C|          L      AAAAA  BBBB   AAAAA    T    C      HHHHH         2           |
C|          L      A   A  B   B  A   A    T    C      H   H   *    2            |
C|          L      A   A  B   B  A   A    T    C   C  H   H  * *  2             |
C|          LLLLL  A   A  BBBB   A   A    T     CCC   H   H   *   22222         |
C|                                                                              |
C|                                                                              |
C|                                                                              |
C|                      (Revision of LABATCH Version 1.0)                       |
C|                                                                              |
C|                           FORTRAN Implementation                             |
C|                                                                              |
C|                                G. S. Fishman                                 |
C|                                                                              |
C|                                October 1997                                  |
C|                                                                              |
C|                                                                              |
C|                       Department of Operations Research                      |
C|                         210 Smith Building  CB # 3180                        |
C|                  University of North Carolina at Chapel Hill                 |
C|                          Chapel Hill, NC 27599-3180                          |
C|                                (919) 962-8401                                |
C|                         email: gfish@fish.or.unc.edu                         |
C|                           fax: (919) 962-0391                                |
C|                                                                              |
C|------------------------------------------------------------------------------|
C|------------------------------------------------------------------------------|
C|                                                                              |
C| LABATCH.2 is described in:                                                   |
C|                                                                              |
C| Fishman, G.S. (1997). LABATCH.2: An enhanced implementation of the batch     |
C| means method, Technical Report No. 97/04,Operations Research Department,     |
C| University of North Carolina at Chapel Hill.                                 |
C|                                                                              |
C| The original LABATCH analysis package is based on:                           |
C|                                                                              |
C| Yarberry, L.S. (1993). Incorporating a dynamic batch size selection mechanism|
C| in a fixed-sample-size batch means procedure. Ph.D. Dissertation, Department |
C| of Operations Research, University of North Carolina at Chapel Hill.         |
C|                                                                              |
C|                                                                              |
C| Use of Version 1.0 is illustrated in:                                        |
C|                                                                              |
C| Fishman, G.S. (1996). Monte Carlo: Concepts, Algorithms, and Applications,   |
C| Springer-Verlag, New York.                                                   |
C|                                                                              |
C| and                                                                          |
C|                                                                              |
C| Fishman, G.S. and L.S. Yarberry (1997). An implementation of the batch means |
C| method, INFORMS Journal on Computing, 9, 296-300.                            |
C|                                                                              |
C|------------------------------------------------------------------------------|
C*------------------------------------------------------------------------------*
C
          integer                 IN_UNIT
          integer                 L_UPPER
          integer                 OUT_UNIT
          integer                 RULE
          integer                 SCREEN
          integer                 T
          integer                 S_NUM
          double precision        BETA
          double precision        DATA(60)             
          double precision        DELTA
          double precision        PSI
C------------------------------------------------------------------------------
C     Read values for input parameters.  Definitions are the same as 
C          in subroutine BATCH_MEANS.
C-----------------------------------------------------------------------------
      open (unit = 11, file = 'lab2par.dat')
      read (unit=11, fmt=*) T,S_NUM,DELTA,RULE,BETA,
     @                      L_UPPER,SCREEN
	IN_UNIT = 30
	OUT_UNIT = 50
      close unit = 11
      open (unit = IN_UNIT, file = 'lab2in.dat')
      open (unit = OUT_UNIT, file = 'lab2out.dat')
      call BATCH_MEANS(IN_UNIT,OUT_UNIT,T,S_NUM,DATA,DELTA,RULE,
     @                 BETA,L_UPPER,SCREEN)
      close unit = OUT_UNIT
      close unit = IN_UNIT
      end
C
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C
C    *******************************************************************
C    *                                                                 *
C    * To reduce the potential for error introduced by unauthorized    *
C    *modifications,it is recommended that this source code be obtained*
C    *from http://www.or.unc.edu/~gfish/labatch.2.html                 *
C    *                                                                 *
C    *******************************************************************
C
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C     Batch Means
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C
      subroutine BATCH_MEANS(IN_UNIT,OUT_UNIT,T,S_NUM,PSI_VECTOR,
     @                       DELTA,RULE,BETA,L_UPPER,SCREEN)


C
C-------------------------------------------------------------------
C
C
C
C-------------------------------------------------------------------
C
C
C
C-------------------------------------------------------------------
C
C 
C     Input Parameters (All are local to subroutine BATCH_MEANS.)
C 
C        IN_UNIT      := 0  if data are to be repeatedly
C                           transferred from a main program
C
C                     := 30 if data are to be read one vector
C                           at a time from input file fort.30
C                           (n.b., 30 is merely an example)
C
C        OUT_UNIT     := designated unit for writing output
C
C        T            := total number of observations                
C 
C        S_NUM        := number of sample series
C 
C 
C        PSI_VECTOR   := vector of S_NUM sample values
C 
C        DELTA        := desired significance level for
C                        confidence intervals  (.01 or .05 suggested)
C 
C        RULE         := 1  if ABATCH Rule is used            
C                     := 0  if LBATCH Rule is used
C 
C        BETA         := significance level for testing for independence
C                        ( .10 suggested)
C                     N.B. If BETA =  1.0, each review uses the FNB rule.
C                                  = -1.0, each review uses the SQRT rule.
C 
C        L_UPPER      := upper bound on initial number of batches
C                        (3 <= L_UPPER <= 100)
C
C        SCREEN       := 1 if interim review estimates are to be
C                          displayed on the screen
C                     := 0 otherwise
C 
C
C     Array Dimensions
C
C        S_MAX        := maximal allowable number of sample series
C        V_MAX        := maximal allowable number of reviews  + 1
C                        (one for the independent case)
C
C-----------------------------------------------------------------------------
C
      integer                      S_MAX
      integer                      V_MAX
      parameter (S_MAX=60, V_MAX=60)
C
C-----------------------------------------------------------------------------
C
C     Important Parameters and Variables
C 
C        B_1              := initial batch size
C        B_2_SQRT         := minimal permissible batch size using
C                            SQRT Rule
C 
C        L_1              := initial number of batches
C        L_2_SQRT         := minimal permissible number of batches
C                            using SQRT Rule
C
C        HEAD             := 1 to produce headings in the computed
C                              tableau file (default is 1 in program)
C                         := 0 to omit headings
C 
C        T_PRIME          := observations 1,...,T_PRIME are used
C                            to compute B*W(L,B) on final review
C 
C        Subroutine SIZE chooses (L_1,B_1) first to maximize T_PRIME
C                        and then to maximize the number of reviews.
C
C-------------------------------------------------------------------
C 
      integer                      B_1
      integer                      B_2_SQRT
      double precision             BETA
      double precision             DELTA
      integer                      HEAD
      integer                      IN_UNIT
      integer                      L_1
      integer                      L_2_SQRT
      integer                      L_UPPER
      integer                      OUT_UNIT
      integer                      RULE
      integer                      SCREEN
      integer                      T
      integer                      T_PRIME
      integer                      S_NUM
C
      integer                      METHOD_VECTOR(S_MAX)
      double precision             BETA_VECTOR(S_MAX)
      double precision             PSI_VECTOR(S_MAX)
C
C-------------------------------------------------------------------
C
C     Working Storage
C
C     In the following definitions, the subscripts I and J are used
C     to denote the value of a quantity for series I on review J.   
C
C        ROW_VECTOR(I)    := number of completed reviews     
C        B_MATRIX(I,J)    := batch size
C        L_MATRIX(I,J)    := number of batches
C        P_MATRIX(I,J)    := p-value
C        V_MATRIX(I,J)    := sample variance of a batch mean         
C                            (W(L,B) in the output tableaus)
C        X_BAR_MATRIX(I,J):= sample mean 
C
C     Statistics for the independent case of series I are maintained
C     in column # ROW_VECTOR(I)+1 of B_MATRIX, L_MATRIX, P_MATRIX,
C     V_MATRIX and X_BAR_MATRIX.
C
C-------------------------------------------------------------------
C
C
      integer                      ROW_VECTOR(S_MAX)
      integer                      B_MATRIX(S_MAX,V_MAX)
      double precision             IND_SUM(S_MAX)
      integer                      L_MATRIX(S_MAX,V_MAX)
      double precision             P_MATRIX(S_MAX,V_MAX)
      double precision             REL(S_MAX)
      double precision             V_MATRIX(S_MAX,V_MAX)
      double precision             X_BAR_MATRIX(S_MAX,V_MAX)
C
      double precision             PHI
      integer                      B
      integer                      B_SQRT
      integer                      I
      integer                      II
      integer                      J
      integer                      J_SQRT
      integer                      JJ
      character                    KK
      integer                      L
      integer                      L_SQRT
      integer                      OLD_ROW        
      integer                      R
      integer                      R_MAX
      integer                      ROW
      integer                      R_SQRT
      integer                      SERIES
      integer                      SWAP_INT
      integer                      SS
      integer                      TEST
      double precision             C
      double precision             RHO
      double precision             S
      double precision             SWAP_DBL
      double precision             TAU
      double precision             X
      double precision             Y
      double precision             Y_SQRT
      integer                      B_VECTOR(S_MAX)
      integer                      B_SQRT_VECTOR(S_MAX)
      integer                      J_VECTOR(S_MAX)
      integer                      J_SQRT_VECTOR(S_MAX)
      integer                      L_VECTOR(S_MAX)
      integer                      L_SQRT_VECTOR(S_MAX)
      integer                      R_SQRT_VECTOR(S_MAX)
      integer                      R_VECTOR(S_MAX)
      integer                      TEST_VECTOR(S_MAX)
      double precision             S_VECTOR(S_MAX)
      double precision             W_IND_VECTOR(S_MAX)
      double precision             Y_IND_VECTOR(S_MAX)
      double precision             Y_SQRT_VECTOR(S_MAX)
      double precision             Y_VECTOR(S_MAX)
      double precision             Z_IND_VECTOR(S_MAX)
      double precision             S_MATRIX(V_MAX,S_MAX)
      double precision             THETA_MATRIX(V_MAX,S_MAX)
      double precision             W_MATRIX(V_MAX,S_MAX)
      double precision             XI_MATRIX(V_MAX,S_MAX)
C
C-------------------------------------------------------------------
C
      HEAD = 0
      if (T_PRIME .eq. 0) then
         do 1  SERIES = 1,S_NUM
            METHOD_VECTOR(SERIES) = RULE
 1          BETA_VECTOR(SERIES)   = BETA
         call SIZE(T,L_UPPER,B_1,B_2_SQRT,L_1,L_2_SQRT,T_PRIME)
      endif
 2    continue
      if (IN_UNIT .gt. 0) then
C-------------------------------------------------------------------
C     The following line has been modified for reading ARENA output.
C     The file reads the variable TNOW.
C-------------------------------------------------------------------
         read(unit=IN_UNIT,fmt=*) (PSI_VECTOR(II), II=1,S_NUM)
C	   read(unit=IN_UNIT,fmt=*) TNOW,(PSI_VECTOR(II), II=1,S_NUM)
      endif
      do 10 SERIES = 1,S_NUM
         if (SS .lt. S_NUM) then
C-------------------------------------------------------------------
C     Begin initialization.
C-------------------------------------------------------------------
            if (SS .eq. 0) then
               if (S_NUM .gt. S_MAX) then
                  print *,'Maximum number of series exceeded.'
                  stop
               endif
               I = 1
            endif
            SS     = SS + 1
            B      = B_1
            B_SQRT = B_2_SQRT
            L      = L_1
            L_SQRT = L_2_SQRT
            J      = 0
            Y      = 0.0D0
            R      = 1
            J_SQRT = 0
            Y_SQRT = 0.0D0
            R_SQRT = 2
            R_MAX = max(
     @                  2*dint(dlog(dfloat(T/B))/dlog(2.0D0))+1,
     @                  2*dint(dlog(dfloat(T/B_SQRT))/dlog(2.0D0))+2)
            if (R_MAX+1 .gt. V_MAX) then
               print *, 'Maximum vector size exceeded.'
               stop
            endif
            TEST = 1
            S    = 0.0D0
            ROW_VECTOR(SERIES) = 0
C-----------------------------------------------------------------
C     End initialization.
C-------------------------------------------------------------------
         else
C-------------------------------------------------------------------
C     Begin restoration.
C-------------------------------------------------------------------
            B      = B_VECTOR(SERIES)
            L      = L_VECTOR(SERIES)
            S      = S_VECTOR(SERIES)
            B_SQRT = B_SQRT_VECTOR(SERIES)
            J      = J_VECTOR(SERIES)
            J_SQRT = J_SQRT_VECTOR(SERIES)
            L_SQRT = L_SQRT_VECTOR(SERIES)
            R      = R_VECTOR(SERIES)
            R_SQRT = R_SQRT_VECTOR(SERIES)
            TEST   = TEST_VECTOR(SERIES)
            Y      = Y_VECTOR(SERIES)
            Y_SQRT = Y_SQRT_VECTOR(SERIES)
C-----------------------------------------------------------------
C     End restoration.
C-------------------------------------------------------------------
      endif
      X = PSI_VECTOR(SERIES)
      if (I .le. T_PRIME) then
         S = S + X
      endif
      IND_SUM(SERIES) = IND_SUM(SERIES) + X
C-------------------------------------------------------------------
C     Update summary data for independent case.
C-------------------------------------------------------------------
         Y_IND_VECTOR(SERIES) = Y_IND_VECTOR(SERIES) + X * X
         if (I .eq. 1) then
               Z_IND_VECTOR(SERIES) = 0.5D0 * Y_IND_VECTOR(SERIES)
            else
               Z_IND_VECTOR(SERIES) = Z_IND_VECTOR(SERIES) + X *
     @            W_IND_VECTOR(SERIES)
         endif
         W_IND_VECTOR(SERIES) = X
         if (mod(I,B_SQRT) .eq. 0) then

C-------------------------------------------------------------------
C     If a batch of size B_SQRT is complete:
C-------------------------------------------------------------------
            J_SQRT = J_SQRT + 1
            call BATCH_UPDATES(J_SQRT, R_SQRT, S,
     @           W_MATRIX(1,SERIES), S_MATRIX(1,SERIES),
     @           THETA_MATRIX(1,SERIES), XI_MATRIX(1,SERIES), R_MAX)
            Y_SQRT = Y_SQRT + W_MATRIX(R_SQRT,SERIES) *
     @            W_MATRIX(R_SQRT,SERIES)
         endif
         if (mod(I,B) .eq. 0) then
C-------------------------------------------------------------------
C     If a batch of size B is complete
C-------------------------------------------------------------------
            J = J + 1
            call BATCH_UPDATES(J, R, S, W_MATRIX(1,SERIES),
     @           S_MATRIX(1,SERIES), THETA_MATRIX(1,SERIES),
     @           XI_MATRIX(1,SERIES), R_MAX)
            Y = Y + W_MATRIX(R,SERIES) * W_MATRIX(R,SERIES)
            if (J .eq. L) then
C-------------------------------------------------------------------
C     If the current review is over, compute statistics.
C-------------------------------------------------------------------
               OLD_ROW                  = ROW_VECTOR(SERIES)
               ROW_VECTOR(SERIES)       = ROW_VECTOR(SERIES) + 1
               ROW                      = ROW_VECTOR(SERIES)
               B_MATRIX(SERIES,ROW)     = B
               L_MATRIX(SERIES,ROW)     = L
               X_BAR_MATRIX(SERIES,ROW) = S / dfloat(B*L)
               TAU                      = (S * S) / dfloat(L)
               RHO                      = Y - TAU
               V_MATRIX(SERIES,ROW) = RHO / (dfloat(B) * dfloat(B) *
     @            dfloat(L-1))
               if (RHO .gt. 0.0D0) then
                     C = (-TAU + THETA_MATRIX(R,SERIES) +
     @                   XI_MATRIX(R,SERIES) + 0.5D0 *
     @                   W_MATRIX(R,SERIES) * W_MATRIX(R,SERIES)) /
     @                   RHO
                     P_MATRIX(SERIES,ROW) = 1.0D0 - 
     @               PHI(C / dsqrt(dfloat(L-2) / (dfloat(L) *
     @               dfloat(L) - 1.0D0)))
                  else
                     P_MATRIX(SERIES,ROW) = 0.0D0
               endif
C-------------------------------------------------------------------
C     End compute statistics.
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C     If J is odd:
C-------------------------------------------------------------------
               if (mod(J,2) .eq. 1) then
                  Y = Y - W_MATRIX(R,SERIES) * W_MATRIX(R,SERIES)
               endif
               Y = Y + 2.0D0 * XI_MATRIX(R,SERIES)
               J = int(J/2)
               R = R + 2
               if (TEST .eq . 1 .and. P_MATRIX(SERIES,ROW) .le.
     @            BETA_VECTOR(SERIES)) then
                  if (B .gt. 1) then
                     B_SQRT = 2 * B_SQRT
C-------------------------------------------------------------------
C     If J_SQRT is odd:
C-------------------------------------------------------------------
                     if (mod(J_SQRT,2) .eq. 1) then
                        Y_SQRT = Y_SQRT - W_MATRIX(R_SQRT,SERIES)
     @                          * W_MATRIX(R_SQRT,SERIES)
                     endif
                     Y_SQRT = Y_SQRT + 2.0D0 *
     @                        XI_MATRIX(R_SQRT,SERIES)
                     J_SQRT = int(J_SQRT / 2)
                     R_SQRT = R_SQRT + 2
                  endif
C-------------------------------------------------------------------
C     FNB Rule
C-------------------------------------------------------------------
                  B = 2 * B
               else
                  if (B .eq. 1) then
C-------------------------------------------------------------------
C     FNB Rule
C-------------------------------------------------------------------
                        B = 2 * B
                     else
                        SWAP_INT = 2 * B
C-------------------------------------------------------------------
C     SQRT Rule
C-------------------------------------------------------------------
                        B        = B_SQRT
                        B_SQRT   = SWAP_INT
                        SWAP_INT = 2 * L
C----------------    -----------------------------------------------
C     SQRT Rule
C-------------------------------------------------------------------
                        L        = L_SQRT
                        L_SQRT   = SWAP_INT
                        SWAP_INT = J
                        J        = J_SQRT
                        J_SQRT   = SWAP_INT
                        SWAP_DBL = Y
                        Y        = Y_SQRT
                        Y_SQRT   = SWAP_DBL
                        SWAP_INT = R
                        R        = R_SQRT
                        R_SQRT   = SWAP_INT
                  endif
                  TEST = METHOD_VECTOR(SERIES)
               endif
            endif
         endif
C-------------------------------------------------------------------
C     Compute statistics for independent case.
C-------------------------------------------------------------------
         if (I .gt. 2) then
            TAU = (IND_SUM(SERIES)*IND_SUM(SERIES)) / dfloat(I)
            RHO = Y_IND_VECTOR(SERIES) - TAU
            V_MATRIX(SERIES,ROW_VECTOR(SERIES)+1) =
     @         RHO/dfloat(I-1)
            if (RHO .gt. 0.0D0) then
               C = (-TAU + Z_IND_VECTOR(SERIES) + 0.5D0 *
     @             W_IND_VECTOR(SERIES) *
     @             W_IND_VECTOR(SERIES)) / RHO
                   P_MATRIX(SERIES,ROW_VECTOR(SERIES)+1) =
     @             1.0D0 - PHI(C/dsqrt(dfloat(I-2)/(dfloat(I)*
     @             dfloat(I)-1.0D0)))
            else
               P_MATRIX(SERIES,ROW_VECTOR(SERIES)+1)=0.0D0
            endif
         endif
C-------------------------------------------------------------------
C     Begin storage.
C-------------------------------------------------------------------
         if( I .lt. T) then
            B_VECTOR(SERIES)      = B
            B_SQRT_VECTOR(SERIES) = B_SQRT
            J_VECTOR(SERIES)      = J
            J_SQRT_VECTOR(SERIES) = J_SQRT
            L_VECTOR(SERIES)      = L
            L_SQRT_VECTOR(SERIES) = L_SQRT
            R_VECTOR(SERIES)      = R
            R_SQRT_VECTOR(SERIES) = R_SQRT
            S_VECTOR(SERIES)      = S
            TEST_VECTOR(SERIES)   = TEST
            Y_VECTOR(SERIES)      = Y
            Y_SQRT_VECTOR(SERIES) = Y_SQRT
         endif
C-------------------------------------------------------------------
C     End storage.
C-------------------------------------------------------------------
10    continue
      if (I .eq. T) then
               call WRITE_TABLEAU_FILE
     @              (DELTA,HEAD,S_NUM,I,T_PRIME,BETA_VECTOR,
     @              IND_SUM,METHOD_VECTOR,ROW_VECTOR,
     @              B_MATRIX,L_MATRIX,P_MATRIX,V_MATRIX,
     @              X_BAR_MATRIX,OUT_UNIT)
      endif
C-------------------------------------------------------------------
C     Display selected interim results on screen.
C-------------------------------------------------------------------
      if (SCREEN .eq. 1) then
         if (ROW .gt. OLD_ROW) then            
            if (OLD_ROW .lt. 1) then
               write (unit=*,fmt=11)
 11            format(///,1X,'LABATCH.2 INTERIM REVIEW STATISTICAL ',
     @                'ANALYSIS                             ',///)
               JJ = min(6,S_NUM)
               write(unit=*,fmt=12) JJ,(II, II=1,JJ)
 12            format(1X,'X_BAR, C.V.(X_BAR), and Sqrt[VAR.EST.'
     @                '] for Series 1 Through ',I2/,
     @                1X,'****************************************',
     @                '**********************',///,
     @                1X,'        No. of',/,
     @                1X,'Review    Obs.                  '6(I1,10X))
               write(unit=*,fmt=122)
 122           format(/)
            endif 
            write(unit=*,fmt=13) ROW,I,(X_BAR_MATRIX(SERIES,ROW),
     @                           SERIES=1,JJ)
            do 125 SERIES = 1,JJ
               REL(SERIES) = 0d0
               if (X_BAR_MATRIX(SERIES,ROW) .ne. 0d0) then
                  REL(SERIES) = dsqrt(V_MATRIX(SERIES,ROW)/
     @                               L_MATRIX(SERIES,ROW))/
     @                               dabs(X_BAR_MATRIX(SERIES,ROW))
               endif
 125        continue
            write(unit=*,fmt=131) (REL(SERIES),SERIES=1,JJ)    
            write(unit=*,fmt=132) (sqrt(B_MATRIX(SERIES,ROW)*
     @                            V_MATRIX(SERIES,ROW)),
     @                            SERIES=1,JJ)
 13         format(I3,I10,1X,'X_BAR         ',(6D11.4)/)
 131        format(14X,'C.V.(X_BAR)   ',(6D11.4)/)
 132        format(14X,'Sqrt[VAR.EST.]',(6D11.4)/)
 14         write (unit=*,fmt=15)
 15         format(1X,'continue[y/n]? ', $ )
            read(unit=*,fmt=*) KK
            if ((KK .ne. 'y').and.(KK .ne. 'n')) then
               go to 14
            endif
C------------------------------------------------------------------
C     If continue = n, write tableau file and terminate execution.
C------------------------------------------------------------------ 
            if (KK .eq. 'n') then
               call WRITE_TABLEAU_FILE
     @              (DELTA,HEAD,S_NUM,I,T_PRIME,BETA_VECTOR,
     @              IND_SUM,METHOD_VECTOR,ROW_VECTOR,
     @              B_MATRIX,L_MATRIX,P_MATRIX,V_MATRIX,
     @              X_BAR_MATRIX,OUT_UNIT)
            endif
            if (KK .ne. 'y') then
               stop
            endif
C-------------------------------------------------------------------
C     End screen display for this review.
C-------------------------------------------------------------------
            OLD_ROW = ROW
         endif
      endif   
      I = I + 1
      if (IN_UNIT .gt. 0) then
         if (I .le. T) then
            go to 2 
         endif
      endif
      return
      end
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C     Batch Updates
C-------------------------------------------------------------------
C-------------------------------------------------------------------
      subroutine BATCH_UPDATES(J_DUMMY, R_DUMMY, S, W_VECTOR,
     @   S_VECTOR, THETA_VECTOR, XI_VECTOR, R_MAX)
C
      integer                      J_DUMMY
      integer                      J
      integer                      R_DUMMY
      integer                      R
      integer                      R_MAX
      double precision             S
      double precision             W
      double precision             W_VECTOR(R_MAX)
      double precision             S_VECTOR(R_MAX)
      double precision             THETA_VECTOR(R_MAX)
      double precision             XI_VECTOR(R_MAX)
C
C
      J = J_DUMMY
      R = R_DUMMY
C-------------------------------------------------------------------
C     While J is even:
C-------------------------------------------------------------------
20    if (mod(J,2) .eq. 0) then
         W            = S - S_VECTOR(R)
         XI_VECTOR(R) = XI_VECTOR(R) + W * W_VECTOR(R)
         W_VECTOR(R)  = W
         S_VECTOR(R)  = S
         J            = J / 2
         R            = R + 2
         go to 20
      endif
      if (J .eq. 1) then
            W               = S
            THETA_VECTOR(R) = 0.5 * W * W
            XI_VECTOR(R)    = 0.0
         else
            W               = S - S_VECTOR(R)
            THETA_VECTOR(R) = THETA_VECTOR(R) + W * W_VECTOR(R)
      endif
      W_VECTOR(R) = W
      S_VECTOR(R) = S
      return
      end
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C     Write tableau file.
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C
      subroutine WRITE_TABLEAU_FILE(DELTA,HEAD,S_NUM,T,T_PRIME,
     @       BETA_VECTOR,IND_SUM,METHOD_VECTOR,ROW_VECTOR,B_MATRIX,
     @       L_MATRIX,P_MATRIX,V_MATRIX,X_BAR_MATRIX,OUT_UNIT)
C
C-------------------------------------------------------------------
C
C     Array Dimensions
C
C        S_MAX       := the maximum number of sample series
C        V_MAX       := the maximum number of reviews + 1     
C                            (one for the independent case)
C
C-------------------------------------------------------------------
C
      integer                      S_MAX
      integer                      V_MAX
C
      parameter (S_MAX=60, V_MAX=60)
C
C-------------------------------------------------------------------
C
C     Input Parameters
C
C        DELTA            := desired significance level for
C                            confidence intervals
C        HEAD             := 1 to produce headings in the computed
C                            tableau file
C                         := 0 to omit headings
C        S_NUM            := number of sample series
C        T                := number of observations
C        T_PRIME          := number of observations used for variance
C                            analysis
C        BETA_VECTOR      := vector of significance levels for
C                            testing independence for each Sample
C        IND_SUM          := sum of the values of the observations
C        METHOD_VECTOR    := vector determining which rule to use
C                            (set to 1 for ABATCH, 0 for LBATCH)
C        OUT_UNIT         := designated unit for writing output
C
C     Working Storage
C
C     In the following definitions, the subscripts I and J are used
C     to denote the value of a quantity for series I on review J.   
C
C        ROW_VECTOR(I)    := number of completed reviews     
C        B_MATRIX(I,J)    := batch size
C        L_MATRIX(I,J)    := number of batches
C        P_MATRIX(I,J)    := p-value
C        V_MATRIX(I,J)    := sample variance of a batch mean          
C                            (W(L,B) in the output tableaus)
C        X_BAR_MATRIX(I,J):= sample mean 
C
C     statistics for the independent case of series I are maintained
C     in column # ROW_VECTOR(I)+1 of B_MATRIX, L_MATRIX, P_MATRIX,
C     V_MATRIX and X_BAR_MATRIX.
C
C-------------------------------------------------------------------
C
      double precision             BETA_VECTOR(S_MAX)
      double precision             IND_SUM(S_MAX)
      integer                      METHOD_VECTOR(S_MAX)
      integer                      ROW_VECTOR(S_MAX)
      integer                      B_MATRIX(S_MAX,V_MAX)
      integer                      L_MATRIX(S_MAX,V_MAX)
      double precision             P_MATRIX(S_MAX,V_MAX)
      double precision             V_MATRIX(S_MAX,V_MAX)
      double precision             X_BAR_MATRIX(S_MAX,V_MAX)
C
      integer                      B
      double precision             DELTA
      double precision             H
      integer                      HEAD
      integer                      L
      integer                      LINES
      integer                      OUT_UNIT
      double precision             REL
      double precision             P
      integer                      ROW
      integer                      S
      integer                      SERIES
      integer                      S_NUM
      integer                      T
      integer                      T_PRIME
      double precision             V
      double precision             X_BAR
C   
      double precision             STUDENT_T_QUANTILE
      if ( HEAD.eq. 1) then
         write (unit=OUT_UNIT,fmt=100)
100      format(/,/,
     @'*---------------------------------------',
     @'---------------------------------------*',/,
     @'|---------------------------------------',
     @'---------------------------------------|',/,
     @'|                                       ',
     @'                                       |',/,
     @'|          L        A    BBBB     A    T',
     @'TTTT   CCC   H   H        222          |',/,
     @'|          L       A A   B   B   A A    ',
     @' T    C   C  H   H       2   2         |',/,
     @'|          L      A   A  B   B  A   A   ',
     @' T    C      H   H          2          |')
         write(unit=OUT_UNIT,fmt=101)
 101     format(
     @'|          L      AAAAA  BBBB   AAAAA   ',
     @' T    C      HHHHH         2           |',/,
     @'|          L      A   A  B   B  A   A   ',
     @' T    C      H   H   *    2            |',/,
     @'|          L      A   A  B   B  A   A   ',
     @' T    C   C  H   H  * *  2             |',/,
     @'|          LLLLL  A   A  BBBB   A   A   ',
     @' T     CCC   H   H   *   22222         |',/,
     @'|                                       ',
     @'                                       |')
         write(unit=OUT_UNIT,fmt=102)
 102     format(
     @'|                                       ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                      (Revision of LABA',
     @'TCH Version 1.0)                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                           FORTRAN Impl',
     @'ementation                             |',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                               G. S. Fi',
     @'shman                                  |',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                               October,',
     @' 1997                                  |',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |')
         write(unit=OUT_UNIT,fmt=103)
 103     format(
     @'|                       Department of Op',
     @'erations Research                      |',/,
     @'|                         210 Smith Buil',
     @'ding  CB # 3180                        |',/,
     @'|                  University of North C',
     @'arolina at Chapel Hill                 |',/,
     @'|                          Chapel Hill, ',
     @'NC 27599-3180                          |',/,
     @'|                                (919) 9',
     @'62-8401                                |',/,
     @'|                         email: gfish@f',
     @'ish.or.unc.edu                         |',/,
     @'|                           fax: (919) 9',
     @'62-0391                                |',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |')
         write(unit=OUT_UNIT,fmt=104)
 104     format(
     @'|---------------------------------------',
     @'---------------------------------------|',/,
     @'|---------------------------------------',
     @'---------------------------------------|',/,
     @'|                                       ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'| LABATCH.2 is described in:            ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'| Fishman, G.S. (1997). LABATCH.2: An en',
     @'hanced implementation of the batch     |',/,
     @'| means method, Technical Report No. 97/',
     @'04, Operations Research Department,    |',/,
     @'| University of North Carolina at Chapel',
     @' Hill.                                 |',/,
     @'|                                       ',
     @'                                       |',/,  
     @'| The original LABATCH analysis package ',
     @'is based on:                           |',/,
     @'|                                       ',
     @'                                       |',/,
     @'| Yarberry, L.S. (1993). Incorporating a',
     @' dynamic batch size selection mechanism|',/,
     @'| in a fixed-sample-size batch means pro',
     @'cedure. Ph.D. Dissertation, Department |',/,
     @'| of Operations Research, University of ',
     @'North Carolina at Chapel Hill.         |')
         write(unit=OUT_UNIT,fmt=105)
 105     format(
     @'|                                       ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'| Use of Version 1.0 is illustrated in: ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'| Fishman, G.S. (1996). Monte Carlo: Con',
     @'cepts, Algorithms, and Applications,   |',/,
     @'| Springer-Verlag, New York.            ',
     @'                                       |')
         write(unit=OUT_UNIT,fmt=106)
 106     format(
     @'|                                       ',
     @'                                       |',/,
     @'| and                                   ',
     @'                                       |',/,
     @'|                                       ',
     @'                                       |',/,
     @'| Fishman, G.S. and L.S. Yarberry (1997)',
     @'. An implementation of the batch means |',/,
     @'| method, INFORMS Journal on Computing, ',
     @'9, 296-310.                            |',/,
     @'|                                       ',
     @'                                       |')
         write(unit=OUT_UNIT,fmt=107)
 107     format(
     @'|---------------------------------------',
     @'---------------------------------------|',/,
     @'*---------------------------------------',
     @'---------------------------------------*',/,/)
      endif
         write (unit=OUT_UNIT,fmt=110) T,(1.0-DELTA)*100.0
 110     format(/,
     @       'Final Tableau                            ',/,
     @       '                                    Mean ',
     @       'Estimation                               ',/,
     @       '                                    *****',
     @       '**********                               ',/,
     @       '                                    (t =',I9,
     @       ' )                                       ',/,/,
     @       '                                         ',
     @       '      ',F4.1, '%                         ',/,
     @       '            _       Standard Error      C',
     @       'onfidence Interval                    _  ',/, 
     @       'Series      X      Sqrt[VAR.EST./t]      ',
     @       'Lower       Upper      (Upper-Lower)/|X| ',/)
      do 28 SERIES = 1, S_NUM
          ROW = ROW_VECTOR(SERIES)
            B = B_MATRIX(SERIES,ROW)
            L = L_MATRIX(SERIES,ROW)
            S = L*B
            if (S .eq. T_PRIME) then
               S = T
            endif
        X_BAR = IND_SUM(SERIES)/dfloat(S)
            V = V_MATRIX(SERIES,ROW)
            H = STUDENT_T_QUANTILE(1.0D0-DELTA/2.0D0, L-1) *
     @         dsqrt(B*V/dfloat(S))
         REL = 0D0
         if (X_BAR .ne. 0D0) then
            REL = 2*H/dabs(X_BAR)
         endif 
 28      write (unit=OUT_UNIT,fmt=140) SERIES,X_BAR,sqrt(B*V/S),
     @         X_BAR-H,X_BAR+H,REL
 140     format (I4,2X,D11.4,4X,D11.4,4X,D11.4,2X,D11.4,6X,D11.4,/)
      write (unit=OUT_UNIT,fmt=145) (100.0 * L*B/T)
 145  format(
     @       '   _                                    ',    
     @       '                                        ',/,
     @       '   X is based on all t observations.    ',/,
     @       '   The Variance Estimator is based on first',F7.2,
     @       '% of the t observations.                ')
      if (HEAD .eq. 1) then
         do 29 LINES = 1,66-13-2*S_NUM                   
            write(unit=OUT_UNIT, fmt=146)
 146        format()
 29      continue
      endif
      do 30 SERIES = 1, S_NUM
         if ( HEAD.eq. 1) then
            if (METHOD_VECTOR(SERIES) .eq. 1) then
               write(unit=OUT_UNIT, fmt=150) SERIES
150            format(
     @               'Interim Review Tableau           ',/,
     @               23X,
     @               'ABATCH Data Analysis for Series ',I2)
            else
               write(unit=OUT_UNIT, fmt=200) SERIES
200           format(
     @               'Interim Review Tableau          ',/,
     @               23X,
     @               'LBATCH Data Analysis for Series ',I2)
         endif
         write (unit=OUT_UNIT,fmt=360) (1. - DELTA)*100.0
 360     format (/,
     @       '                                         ',
     @       '          ',F4.1, '%                     ',/,
     @       '                                    _    ',
     @       '   Confidence Interval                   ',/,
     @       'Review   B*M        B        M      X    ',
     @       '    Lower       Upper Sqrt[VAR.EST.] p-va',
     @       'lue                                       ',/)
         endif
         do 40 ROW = 1,ROW_VECTOR(SERIES)
            B     = B_MATRIX(SERIES,ROW)
            L     = L_MATRIX(SERIES,ROW)
            P     = P_MATRIX(SERIES,ROW)
            X_BAR = X_BAR_MATRIX(SERIES,ROW) 
            S     = L*B
            if (S .eq. T_PRIME) then
               S     = T
               X_BAR = IND_SUM(SERIES)/dfloat(S)
            endif
            V     = V_MATRIX(SERIES,ROW)
            H     = STUDENT_T_QUANTILE(1.0D0-DELTA/2.0D0, L-1) *
     @              dsqrt(B*V/dfloat(S))
            write (unit=OUT_UNIT,fmt=425) ROW,B*L, L,B,X_BAR,
     @             X_BAR-H,X_BAR+H,dsqrt(dfloat(B)*V),sngl(P)
425         format(I3,1X,I8,1X,I8,1X,I8,1X,3D11.4,1X,D11.4,2X,F6.4)
40       continue
         if ( HEAD.eq. 1) then
            write (unit=OUT_UNIT,fmt=450)
450         format (/,'  If data are independent:',/)
         endif
         L = B * L
         B = 1
         P = P_MATRIX(SERIES,ROW_VECTOR(SERIES)+1)
         V = V_MATRIX(SERIES,ROW_VECTOR(SERIES)+1)
         H = STUDENT_T_QUANTILE(1.0D0-DELTA/2.0D0, S-1) *
     @              dsqrt(B*V/dfloat(S))
         write (unit=OUT_UNIT,fmt=500) T,T,B,X_BAR,X_BAR-H,X_BAR+H,
     @         dsqrt(dfloat(B)*V),sngl(P)
500      format (4X,I8,1X,I8,1X,I8,1X,3D11.4,1X,D11.4,2X,F6.4)
         write (unit=OUT_UNIT,fmt=525) max(0.0D0,BETA_VECTOR(SERIES))
525      format (/,1X,F5.2,' significance level for independence',
     @              ' testing.')
         write (unit=OUT_UNIT,fmt=535) ROW-1,(100.0*L*B/T)
535      format ('  Review',I3,
     @              ' used the first',F7.2,
     @              '% of the t observations for variance estimation.')    
         if (HEAD .eq. 1) then
            do 50 LINES = 1,66-15-ROW_VECTOR(SERIES)      
               write(unit=OUT_UNIT, fmt=550)
550            format()
50          continue
         endif
30    continue
      return
      end
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C     PHI function
C-------------------------------------------------------------------
C-------------------------------------------------------------------
      double precision function PHI(X)
C
C-------------------------------------------------------------------
C  This function computes the value of the normal c.d.f. at the
C  point X.
C
C  Reference: Abramowitz, M. and I.A. Stegun (1964). Handbook of
C  Mathematical Functions with Formulas, Graphs and Mathematical
C  Tables, U.S. Government Printing Office, Washington, D.C.
C-------------------------------------------------------------------
C
      double precision             X
C
      double precision             Y
C     
      if (X .gt. 0.0D0) then
            if (X .gt. 45.0D0) then
                  PHI = 1.0D0
               else
                  Y   = X
                  PHI = 1.0D0-0.50/(1.0D0+(0.0498673470D0+
     @               (0.0211410061D0+(0.0032776263D0+(0.0000380036D0
     @               +(0.0000488906D0+0.0000053830D0*Y)*Y)*Y)*Y)*Y)
     @               *Y)**16
            endif
         else
            if (X .lt. -45.0D0) then
                  PHI = 0.0D0
               else
                  Y   = -X
                  PHI = 0.50/(1.0D0+(0.0498673470D0+(0.0211410061D0+
     @               (0.0032776263D0+(0.0000380036D0+(0.0000488906D0
     @               +0.0000053830D0*Y)*Y)*Y)*Y)*Y)*Y)**16
            endif
      endif
      return
      end
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C     PHI quantile function
C-------------------------------------------------------------------
C-------------------------------------------------------------------
      double precision function PHI_QUANTILE(BETA)
C
C-------------------------------------------------------------------
C  This function computes the quantile of order BETA from the normal
C  distribution.
C
C  Reference: Abramowitz, M. and I.A. Stegun (1964). Handbook of
C  Mathematical Functions with Formulas, Graphs and Mathematical
C  Tables, U.S. Government Printing Office, Washington, D.C.
C-------------------------------------------------------------------
C
      double precision             BETA
C
      double precision             PHI
C
      double precision             DELTA
      double precision             HIGH
      double precision             LOW
      double precision             MID
      double precision             W
C
      if (BETA .ge. 0.5D0) then
            DELTA = BETA
         else
            DELTA = 1.0D0-BETA
      endif
      W   = dsqrt(-dlog((1.0D0-DELTA)*(1.0D0-DELTA)))
      LOW = W-(2.515517D0+(0.802853D0+0.010328D0*W)*W)/
     @   (1.0D0+(1.432788D0+(0.189269D0+0.001308D0*W)*W)*W)
      HIGH = LOW + 4.5D-4
      LOW  = LOW - 4.5D-4
 60      MID = (LOW + HIGH)/2.0D0
         if (PHI(MID) .lt. DELTA) then
               LOW = MID
            else
               HIGH = MID
         endif
         if (dabs(MID-(LOW+HIGH)/2.0D0) .gt. 1.0D-15) go to 60
      if (BETA .ge. 0.5D0) then
            PHI_QUANTILE = MID
         else
            PHI_QUANTILE = -MID
      endif
      return
      end
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C     Student t quantile function
C-------------------------------------------------------------------
C-------------------------------------------------------------------
      double precision function STUDENT_T_QUANTILE(BETA,L)
C
C-------------------------------------------------------------------
C  This function computes the quantile of order BETA from the
C  Student t distribution with L degrees of freedom.
C
C  Reference: Abramowitz, M. and I.A. Stegun (1964). Handbook of
C  Mathematical Functions with Formulas, Graphs and Mathematical
C  Tables, U.S. Government Printing Office, Washington, D.C.
C-------------------------------------------------------------------
C
      double precision             BETA
      integer                      L
C
      double precision             PHI_QUANTILE
C
      double precision             DELTA
      double precision             DFL
      double precision             HIGH
      double precision             LOW
      double precision             MID
      double precision             TEST
      double precision             W
      double precision             X
      double precision             Y
C
      if (BETA .eq. 0.5D0) then
            STUDENT_T_QUANTILE = 0.0D0
         elseif (L .ge. 7) then
            X = PHI_QUANTILE(BETA)
            Y = X*X
            DFL = dfloat(L)
            STUDENT_T_QUANTILE = X*(1.0D0+
     @      ((-945.0D0+(-3600.0D0+(2880.0D0+23040.0D0*DFL)*DFL)*DFL)
     @      +((-1920.0D0+(4080.0D0+(15360.0D0+23040.0D0*DFL)*DFL)*
     @      DFL)+((1482.0D0+(4560.0D0+4800.0D0*DFL)*DFL)+((776.0D0
     @      +720.0D0*DFL)+79.0D0*Y)*Y)*Y)*Y)/(DFL**4*92160.0D0))
         elseif (L .eq. 1) then
            STUDENT_T_QUANTILE = dtan((BETA-0.5D0)*
     @         3.1415926535898D0)
         elseif (L .eq. 2) then
            STUDENT_T_QUANTILE = (2.0D0*BETA-1.0D0)/
     @         dsqrt(2.0D0*BETA*(1.0D0-BETA))
         else
            if (BETA .ge. 0.5D0) then
                  DELTA = BETA
               else
                  DELTA = 1.0D0-BETA
            endif
            LOW  = PHI_QUANTILE(DELTA)
            HIGH = (2.0D0*DELTA-1.0D0)/
     @         dsqrt(2.0D0*DELTA*(1.0D0-DELTA))
 70         MID  = (LOW + HIGH)/2.0D0
               W = MID * MID
               if (L .eq. 3) then
                     TEST = dsqrt(3.0D0)*dtan(3.1415926535898D0*
     @                  (DELTA-0.5D0)-(MID*dsqrt(3.0D0))/(3.0D0+W))
                  elseif (L .eq. 4) then
                     TEST = (2.0D0*DELTA-1.0D0)*
     @                  dsqrt(4.0D0+W)**3/(6.0D0+W)
                  elseif (L .eq. 5) then
                     TEST = dsqrt(5.0D0)*dtan(3.1415926535898D0*
     @                  (DELTA-0.5D0)-(MID*dsqrt(5.0D0))/(3.0D0*
     @                  (5.0D0+W)**2)*(25.0D0+3.0D0*W))
                  elseif (L .eq. 6) then
                     TEST = (2.0D0*DELTA-1.0D0)*dsqrt(6.0D0+W)**5/
     @                  (W*W+15.0D0*w+67.5D0)
               endif
               if (TEST .gt. MID) then
                     LOW = MID
                  else
                     HIGH = MID                  
               endif
               if (dabs(MID-(LOW+HIGH)/2.0D0) .gt. 1.0D-15) go to 70
            if (BETA .ge. 0.5D0) then
                  STUDENT_T_QUANTILE = MID
               else
                  STUDENT_T_QUANTILE = -MID
            endif
      endif
      return
      end
C-----------------------------------------------------------------
C-----------------------------------------------------------------
C subroutine SIZE finds L_1 and B_1 first to maximize T_PRIME
C            (<=100) and then to maximize the number of reviews.
C-----------------------------------------------------------------
C-----------------------------------------------------------------
C
      subroutine SIZE(T,L_UPPER,B_1,B_2_SQRT,L_1,L_2_SQRT,T_PRIME)
         integer        ALPHA
         integer        ALPHA_MAX
         integer        B_1
         integer        B_2_SQRT
         integer        I
         integer        I_MAX
         integer        L_1
         integer        L_2_SQRT
         integer        L_UPPER
         integer        T
         integer        T_PRIME
         integer        TEMP_1
         integer        TEMP_2
      integer B(46) /1,5,5,7,7,7,8,9,12,12,12,15,17,17,17,19,
     @               21,21,22,22,26,27,29,31,32,33,36,36,36,
     @               39,41,43,46,49,50,51,51,53,56,60,63,66,
     @               67,70,84,85/
      integer L(46) /3,7,21,15,25,35,11,13,17,51,85,21,36,60,
     @               84,27,25,35,31,93,37,57,41,66,45,47,51,
     @               68,85,55,87,61,65,69,71,60,84,75,79,85,
     @               89,93,95,99,85,96/
C 
      I_MAX     = 1
      ALPHA_MAX = 0
      TEMP_2    = 0
      do 2100 I = 1, 46
         if ((T .ge. B(I)*L(I)) .and. (L(I) .le. L_UPPER)) then
            ALPHA  = idint(dlog(dfloat(T/ L(I)/B(I)))/dlog(2.0D0))
            TEMP_1 = 2**ALPHA*B(I)* L(I) 
            if (TEMP_1 .ge. TEMP_2) then
               if ((TEMP_1 .gt. TEMP_2)
     @             .or. (L(I)*B(I) .lt. L(I_MAX)*B(I_MAX))) then
                  I_MAX     = I
                  ALPHA_MAX = ALPHA
                  TEMP_2    = TEMP_1
               endif
            endif
         endif
 2100 continue
      B_1 = B(I_MAX)
      B_2_SQRT = 3
      if (B_1 .gt. 1) then 
         B_2_SQRT = idint(dsqrt(2.0D0)*B_1 + .5)
      endif
      L_1 =  L(I_MAX)
      L_2_SQRT = idint(dsqrt(2.0D0)*L_1 + .5)
      T_PRIME = TEMP_2
      return
      end
