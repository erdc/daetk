C Work performed under the auspices of the U.S. Department of Energy
C by Lawrence Livermore National Laboratory under contract number 
C W-7405-Eng-48.
C
C Copyright 1998 the Regents of the University of California.
C All rights reserved.
C
C***BEGIN PROLOGUE  DWEBS
C***REFER TO  DDASPK
C***DATE WRITTEN   950914   (YYMMDD)
C***MODIFIED TO INCLUDE SENSITIVITY ANALYSIS 990115
C
C***AUTHORS  A. C. Hindmarsh, P. N. Brown
C            Lawrence Livermore National Laboratory
C            Livermore, CA 94551, USA
C
C            L. R. Petzold
C            University of California
C            Santa Barbara, CA  93106, USA
C
C            Shengtai Li
C            University of California
C            Santa Barbara, CA  93106, USA
C
C***DESCRIPTION
C
C-----------------------------------------------------------------------
C Example program for DDASPK.
C DAE system derived from ns-species interaction PDE in 2 dimensions.
C Sensitivity analysis with respect to parameters alpha and beta.
C
C This is the double precision version.
C-----------------------------------------------------------------------
C
C This program solves a DAE system that arises from a system
C of partial differential equations.  The PDE system is a food web
C population model, with predator-prey interaction and diffusion on
C the unit square in two dimensions.  The dependent variable vector is
C
C         1   2        ns
C   c = (c , c , ..., c  )
C
C and the PDEs are as follows..
C
C     i               i      i
C   dc /dt  =  d(i)*(c    + c   )  +  R (x,y,c)  (i=1,...,ns/2)
C                     xx     yy        i
C
C                     i      i
C   0       =  d(i)*(c    + c   )  +  R (x,y,c)  (i=(ns/2)+1,...,ns)
C                     xx     yy        i
C
C where
C                  i          ns         j
C   R (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
C    i                       j=1
C
C The number of species is ns = 2*np, with the first np being prey and
C the last np being predators.  The coefficients a(i,j), b(i), d(i) are
C
C   a(i,i) = -a  (all i)
C   a(i,j) = -g  (i .le. np, j .gt. np)
C   a(i,j) =  e  (i .gt. np, j .le. np)
C   b(i) =  b*(1 + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y))  (i .le. np)
C   b(i) = -b*(1 + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y))  (i .gt. np)
C   d(i) = dprey  (i .le. np)
C   d(i) = dpred  (i .gt. np)
C
C The various scalar parameters are set in subroutine setpar.
C
C The boundary conditions are.. normal derivative = 0.
C A polynomial in x and y is used to set the initial conditions.
C
C The PDEs are discretized by central differencing on a MX by MY mesh.
C
C The DAE system is solved by DDASPK with three different method options:
C (1) direct band method for the linear systems (internal Jacobian),
C (2) preconditioned Krylov method for the linear systems, without
C     block-grouping in the reaction-based factor, and
C (3) preconditioned Krylov method for the linear systems, with
C     block-grouping in the reaction-based factor.
C
C In the Krylov cases, the preconditioner is the product of two factors:
C (a) The spatial factor uses a fixed number of Gauss-Seidel iterations
C based on the diffusion terms only.
C (b) The reaction-based factor is a block-diagonal matrix based on
C the partial derivatives of the interaction terms R only.
C With block-grouping, only a subset of the ns by ns blocks are computed.
C An integer flag, JPRE, is set in the main program to specify whether
C the preconditioner is to use only one of the two factors or both,
C and in which order.
C
C The reaction-based preconditioner factor is set up and solved in
C seven subroutines -- 
C   DMSET2, DRBDJA, DRBDPS  in the case of no block-grouping, and
C   DGSET2, GSET1, DRBGJA, DRBGPS  in the case of block-grouping.
C These routines are provided separately for general use on problems
C arising from reaction-transport systems.
C
C Two output files are written.. one with the problem description and
C performance statistics on unit LOUT = 9, and one with solution 
C profiles at selected output times on unit LCOUT = 10.
C The solution file is written only in the case of the direct method.
C-----------------------------------------------------------------------
C Note.. in addition to the main program and subroutines given below,
C this program requires the BLAS routine DAXPY.
C-----------------------------------------------------------------------
C References
C [1] Peter N. Brown and Alan C. Hindmarsh,
C     Reduced Storage Matrix Methods in Stiff ODE Systems,
C     J. Appl. Math. & Comp., 31 (1989), pp. 40-91.
C [2] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
C     Using Krylov Methods in the Solution of Large-Scale Differential-
C     Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488.
C [3] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
C     Consistent Initial Condition Calculation for Differential-
C     Algebraic Systems, LLNL Report UCRL-JC-122175, August 1995;
C     submitted to SIAM J. Sci. Comp.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   SETPAR, DGSET2, CINIT, DDASPK, OUTWEB
C       
C***END PROLOGUE  DWEB
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RESWEB, DJACILU, DPSOLILU, H_RES
C
C Set output unit numbers for main output and tabulated solution.
      DATA LOUT/6/, LCOUT/10/
C Load parameters for sparse preconditioner.
      PARAMETER (LENPFAC = 10)
      PARAMETER (LENPLUFAC = 2)
      PARAMETER (IPREMETH = 2)  ! =1 means ILUT preconditioner used
                                ! =2 means ILUTP preconditioner used
      PARAMETER (LFILILUT = 10)
      PARAMETER (IREORDER = 1)
      PARAMETER (ISRNORM = 1)
      PARAMETER (NORMTYPE = 2)
      PARAMETER (JACOUT = 0)
      PARAMETER (JSCALCOL = 1)
      PARAMETER (TOLILUT = 0.001)
      PARAMETER (PERMTOL = 0.01)
      PARAMETER (NPARA = 2)
      PARAMETER (MAXS = 2, MAXNY = 800, MAXN = (NPARA+1)*MAXNY)
      PARAMETER (NRPD = 2+MAXNY)
      PARAMETER (LENWP = 3*LENPFAC*MAXNY + LENPLUFAC*MAXNY
     1                   + ISRNORM*MAXNY + 2*(MAXNY+1) )
      PARAMETER (LENIWP = 5*(MAXNY+1) + 4*LENPFAC*MAXNY
     1                    + 2*LENPLUFAC*MAXNY + IREORDER*2*MAXNY
     2                    + (IPREMETH-1)*2*MAXNY + 3*MAXNY +NRPD)
      PARAMETER (LRW = 67301, LIW = 46447)
C 
      DIMENSION CC(MAXN), CCPRIME(MAXN), RWORK(LRW), IWORK(LIW),
     1     INFO(30), RPAR(NRPD), IPAR(31), senpar(2)
C
C The COMMON blocks /PPAR1/ and /PPAR2/ contain problem parameters.
C
      COMMON /PPAR1/ AA, EE, GG, BB, DPREY, DPRED
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY,  FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
      real t1, t2, t3, tt(2)
C

C Open output files.
c      OPEN(unit=lout,file='wdout',status='unknown')
c      OPEN(unit=lcout,file='wccout',status='unknown')
C
C Call SETPAR to set basic problem parameters.
      CALL SETPAR(senpar)
C
      ALPH = senpar(1)
      BETA = senpar(2)
C Set remaining problem parameters.
      NPARM = 2
      NEQ = NS*MX*MY
      MXNS = MX*NS
      DX = AX/REAL(MX-1)
      DY = AY/REAL(MY-1)
      DO 10 I = 1,NS
        COX(I) = DIFF(I)/DX**2
 10     COY(I) = DIFF(I)/DY**2
C
      WRITE(LOUT,20)NS
 20   FORMAT(//' Example program for DDASPK package'//
     1   ' Food web problem with NS species, NS =',I4/
     2   ' Predator-prey interaction and diffusion on a 2-D square'/)
      WRITE(LOUT,30) AA,EE,GG,BB,DPREY,DPRED, ALPH,BETA
 30   FORMAT(' Matrix parameters..  a =',E12.4,'   e =',E12.4,
     1   '   g =',E12.4/21x,' b parameter =',E12.4//
     2   ' Diffusion coefficients.. dprey =',E12.4,'   dpred =',E12.4//
     3   ' Rate parameters alpha =',E12.4,' and beta =',E12.4/)
      WRITE(LOUT,40) MX,MY,NEQ
 40   FORMAT(' Mesh dimensions (MX,MY) =',2I4,
     1   5x,' Total system size is NEQ =',I7/)
C
C Here set the flat initial guess for the predators.
	PREDIC = 1.0D5
C 
C Load values into IPAR and RPAR for sparse preconditioner.
      ML = NS*MX + 1
      MU = ML
      IPAR(1) = ML
      IPAR(2) = MU
      IPAR(3) = LENPFAC
      IPAR(4) = LENPLUFAC
      IPAR(5) = IPREMETH
      IPAR(6) = LFILILUT
      IPAR(7) = IREORDER
      IPAR(8) = ISRNORM
      IPAR(9) = NORMTYPE
      IPAR(10) = JACOUT
      IPAR(11) = JSCALCOL
      IPAR(12) = 1          ! Jacobian evaluation method
      print *, ' Jacobian evaluation for ILU preconditioner:(0,1)'
      print *, '    0 --- finite differencing '
      print *, '    1 --- ADIFOR '
      read(*,*) i
      if (i .eq. 0) ipar(12) = 0
      IPAR(13) = NRPD          ! data dependence for rpar(*)
      IPAR(30) = 0
      RPAR(1) = TOLILUT
      RPAR(2) = PERMTOL
C Check IPAR, RPAR, LENWP and LENIWP for illegal entries and long
C enough work array lengths.
      CALL DSPSETUP (NEQ, LENWP, LENIWP, RPAR, IPAR, IERR,
     .               LWPMIN, LIWPMIN)
      IF (IERR .NE. 0) THEN
         WRITE(LOUT,45) IERR
 45      FORMAT(' Error return from DSPSETUP: IERR = ',i5)
         IF (LWPMIN .GT. LENWP) THEN
            WRITE(LOUT,*) ' More WP work array length needed'
         ENDIF
         IF (LIWPMIN .GT. LENIWP) THEN
            WRITE(LOUT,*) ' More IWP work array length needed'
         ENDIF
         STOP
      ENDIF
C
C Set remaining method parameters for DDASPK.
C These include the INFO array and tolerances.
C 
      DO 50 I = 1,30
 50   INFO(I) = 0
C
C Here set INFO(11) = 1, indicating I.C. calculation requested.
      INFO(11) = 1
C
C Here set INFO(18) = 2, indicating full print for initial condition
C      calculation
      INFO(18) = 0
C
C Here set INFO(14) = 1 to get the computed initial values.
      INFO(14) = 1
C
C Here set INFO(15) = 1 to signal that a preconditioner setup routine
C is to be called in the Krylov case.
      INFO(15) = 1
C
C Here set INFO(16) = 1 to get alternative error test (on the
C differential variables only).
      INFO(16) = 1
C
C
C Here set INFO(19) = 2 to indicate that a sensitivity analysis
C is to be done with two parameters
      inso = 18
      info(inso+1) = 2
      print *, ' number of sensitivities? (0,1,2) '
      read(*,*) info(inso+1)
      Ny = neq
      neq = neq*(info(inso+1) + 1)
      if (info(inso+1) .eq. 0) go to 60
      print *, ' method for computing sensitivity residues:(0,1,3)'
      print *, '    0: central differencing '
      print *, '    1: one-side forward differencing '
      print *, '    3: Adifor '
      read (*,*) info(inso+2)
      if (info(inso+2) .eq. 2) info(inso+2) = 3     
      iwork(38) = nrpd
      info(inso+3) = 0          ! default perturbed factor
      info(inso+4) = 2          ! number of sensitivity parameters in RES
      info(inso+5) = 0          ! 0, full control 1. partial control
      print *, ' Error control: (0,1)'
      print *, '    0: full error control with sensitivity '
      print *, '    1: partial error control on state variables only '
      read(*, *) info(inso+5)
      info(inso+6) = 0          ! no derived information
      info(inso+7) = 1          ! simultaneous corrector(0), staggered(1)
      print *, ' Corrector method: (0, 1)'
      print *, '    0:  simultaneous corrector method '
      print *, '    1:  staggered corrector method'
      read(*,*) info(inso+7)
      RWORK(16) = 1.d-5
 60   continue
C
C
C Here set the tolerances.      
      RTOL = 1.0D-5
      ATOL = RTOL
      WRITE(LOUT,70)RTOL,ATOL,INFO(11),PREDIC,INFO(16)
 70   FORMAT(' Tolerance parameters.. RTOL =',E10.2,'   ATOL =',E10.2//
     1   ' Internal I.C. calculation flag INFO(11) =',I2,
     2   '   (0 = off, 1 = on)'/
     3   ' Predator I.C. guess =',E10.2//
     4   ' Alternate error test flag INFO(16) =',I2,
     5        '  (0 = off, 1 = on)')
C
C Set NOUT = number of output times.
      NOUT = 18
C
C Loop over method options: 
C METH = 0 means use INFO(12) = 0 (direct)
C METH = 1 means use INFO(12) = 1 (Krylov) without block-grouping in
C          the reaction-based factor in the preconditioner.
C METH = 2 means use INFO(12) = 1 (Krylov) with block-grouping in
C          the reaction-based factor in the preconditioner.
C A block-grouping flag JBG, communicated through IPAR, is set to
C 0 (no block-grouping) or 1 (use block-grouping) with METH = 1 or 2.
C Reset INFO(1) = 0 and INFO(11) = 1.
C
C
C Set and print various preconditioner parameters.
      INFO(12) = 1
      IWORK(27) = LENWP 
      IWORK(28) = LENIWP
      WRITE(LOUT,100) IPREMETH
 100  FORMAT(' Preconditioner flag is IPREMETH =',I3,
     1       '  (1 = ILUT, 2 = ILUTP)')
C Here call SETID to set the IWORK segment ID  indicating the
C differential and algebraic components.
      CALL SETID (MX, MY, NS, NP, 40, IWORK)
C
C Set the initial T and TOUT, and call CINIT to set initial values.
      T = 0.0D0
      TOUT = 1.0D-8
      CALL CINIT (CC, CCPRIME, PREDIC, rpar, senpar)
C
      NLI = 0
      NNI = 0
C
      WRITE(LOUT,140)
 140  FORMAT(//'   t',7X,'NSTEP  NRE  NNI  NLI  NPE  NQ',4X,
     1   'H',10X,'AVLIN',4X, 'TIME')
C
C Loop over output times, call DDASPK, and print performance data.
C The first call, with IOUT = 0, is to calculate initial values only.
C After the first call, reset INFO(11) = 0 and the initial TOUT.
C
      t1 = secnds(0.0e0)
      t2 = dtime(tt)
      DO 200 IOUT = 0,NOUT
C
         CALL DDASPK (
     1        RESWEB, NEQ, T, CC, CCPRIME, TOUT, INFO, RTOL, ATOL, 
     1        IDID, RWORK,LRW, IWORK,LIW, RPAR, IPAR, 
     1        DJACILU, DPSOLILU, senpar, h_res)
C
c        write(12,*) cc(ny/2+NY), cc(ny/2+1+NY)
        NST = IWORK(11)
        NRE = IWORK(12)
        NPE = IWORK(13)
        NNIDIF = IWORK(19) - NNI
        NNI = IWORK(19)
        NLIDIF = IWORK(20) - NLI
        NLI = IWORK(20)
        NQU = IWORK(8)
        HU = RWORK(7)
        AVLIN = 0.0D0
        IF (NNIDIF .GT. 0) AVLIN = REAL(NLIDIF)/REAL(NNIDIF)
C
        IF (METH .EQ. 0) THEN
          IMOD3 = IOUT - 3*(IOUT/3)
c          IF (IMOD3 .EQ. 0) CALL OUTWEB (T, CC, NS, MX, MY, LCOUT)
          ENDIF
C
        WRITE(LOUT,150)T,NST,NRE,NNI,NLI,NPE,NQU,HU,
     *         AVLIN,secnds(t1)
 150    FORMAT(E10.2,I5,I6,3I5,I4,E11.2,F9.4, f8.2)
C
        IF (IDID .LT. 0) THEN
          WRITE(LOUT,160)T
 160      FORMAT(//' Final time reached =',E12.4//)
          GO TO 210
          ENDIF
C
        IF (TOUT .GT. 0.9D0) TOUT = TOUT + 1.0D0
        IF (TOUT .LT. 0.9D0) TOUT = TOUT*10.0D0
        IF (IOUT .EQ. 0) THEN
c           write(10, *) (cc(j), j=ny+1, 2*ny)
          INFO(11) = 0
          TOUT = 1.0D-8
          NLI = 0
          NNI = 0
          ENDIF
 200    CONTINUE
C
 210  CONTINUE
      LENRW = IWORK(18)
      LENIW = IWORK(17)
      NST = IWORK(11)
      NRE = IWORK(12)
      NPE = IWORK(13)
      NNI = IWORK(19)
      NLI = IWORK(20)
      NPS = IWORK(21)
      NSE = IWORK(22)
      IF (NNI .GT. 0) AVLIN = REAL(NLI)/REAL(NNI)
      NETF = IWORK(14)
      NCFN = IWORK(15)
      NCFL = IWORK(16)
      WRITE(LOUT,220) LENRW,LENIW,NST,NRE+IPAR(30),NSE,NPE,NPS,NNI,
     1               NLI,AVLIN,NETF,NCFN,NCFL
 220   FORMAT(//' Final statistics for this run..'/
     1   ' RWORK size =',I8,'   IWORK size   =',I6/
     2   ' Number of time steps              =',I5/
     3   ' Number of residual evaluations    =',I5/
     3   ' Number of sensitivity evaluations =',I5/
     4   ' Number of Jac. or prec. evals.    =',I5/
     5   ' Number of preconditioner solves   =',I5/
     6   ' Number of nonlinear iterations    =',I5/
     7   ' Number of linear iterations       =',I5/
     8   ' Average Krylov subspace dimension =',F8.4/
     1   ' error test failure .............. =',I5/
     9   I3,' nonlinear conv. failures,',i5,' linear conv. failures')
C
       t3 = dtime(tt)
       print *, ' Time used: ', secnds(t1), tt(1)+tt(2)
 300   CONTINUE

      STOP
C------  End of main program for DWEB example program ------------------
      END

      SUBROUTINE SETPAR(senpar)
C-----------------------------------------------------------------------
C This routine sets the basic problem parameters, namely
C AX, AY, NS, MX, MY,  problem coefficients ACOEF, BCOEF, DIFF,
C ALPH, BETA, using parameters NP, AA, EE, GG, BB, DPREY, DPRED.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision senpar(2)
      PARAMETER (MAXS = 2)
      COMMON /PPAR1/ AA, EE, GG, BB, DPREY, DPRED
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY,  FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      AX = 1.0D0
      AY = 1.0D0
      NP = 1
      MX = 20
      MY = 20
      AA = 1.0D0
      EE = 1.0D4
      GG = 0.5D-6
      BB = 1.0D0
      DPREY = 1.0D0
      DPRED = 0.05D0
      ALPH = 50.0D0
      BETA = 100.0D0
      senpar(1) = alph
      senpar(2) = beta
      NS = 2*NP
      DO 20 J = 1,NP
        DO 10 I = 1,NP
          ACOEF(NP+I,J) = EE
          ACOEF(I,NP+J) = -GG
 10       CONTINUE
        ACOEF(J,J) = -AA
        ACOEF(NP+J,NP+J) = -AA
        BCOEF(J) = BB
        BCOEF(NP+J) = -BB
        DIFF(J) = DPREY
        DIFF(NP+J) = DPRED
 20     CONTINUE
      PI = 3.141592653589793D0
      FPI = 4.0D0*PI
C
      RETURN
C------------  End of Subroutine SETPAR  -------------------------------
      END

      SUBROUTINE SETID (MX, MY, NS, NSD, LID, IWORK)
C-----------------------------------------------------------------------
C This routine sets the ID array in IWORK, indicating which components
C are differential and which are algebraic.
C-----------------------------------------------------------------------
      DIMENSION IWORK(*)
C
      NSDP1 = NSD + 1
      DO 40 JY = 1,MY
        I00 = MX*NS*(JY-1) + LID
        DO 30 JX = 1,MX
          I0 = I00 + NS*(JX-1)  
          DO 10 I = 1,NSD
 10         IWORK(I0+I) = 1
          DO 20 I = NSDP1,NS
 20         IWORK(I0+I) = -1
 30       CONTINUE
 40     CONTINUE
C
      RETURN
C------------  End of Subroutine SETID  --------------------------------
      END

      SUBROUTINE CINIT (CC, CCPRIME, PREDIC, RPAR, senpar)
C-----------------------------------------------------------------------
C This routine computes and loads the vectors of initial values.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(*), CCPRIME(*), RPAR(*), senpar(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY,  FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
C Load CC.
      alph = senpar(1)
      beta = senpar(2)
      NPP1 = NP + 1
      DO 30 JY = 1,MY
        Y = REAL(JY-1)*DY
        ARGY = 16.0D0*Y*Y*(AY-Y)*(AY-Y)
        IYOFF = MXNS*(JY-1)
        DO 20 JX = 1,MX
          X = REAL(JX-1)*DX
          ARGX = 16.0D0*X*X*(AX-X)*(AX-X)
          IOFF = IYOFF + NS*(JX-1)
          FAC = 1.0D0 + ALPH*X*Y + BETA*SIN(FPI*X)*SIN(FPI*Y)
          DO 10 I = 1,NP
 10         CC(IOFF + I) = 10.0D0 + REAL(I)*ARGX*ARGY
          DO 15 I = NPP1,NS
 15         CC(IOFF + I) = PREDIC
 20       CONTINUE
 30     CONTINUE
C
C Load CCPRIME.
      T = 0.0D0
      CALL FWEB (T, CC, CCPRIME, RPAR, senpar)
      DO 60 JY = 1,MY
        IYOFF = MXNS*(JY-1)
        DO 50 JX = 1,MX
          IOFF = IYOFF + NS*(JX-1)
          DO 40 I = NPP1,NS
 40         CCPRIME(IOFF+I) = 0.0D0
 50     CONTINUE
 60   CONTINUE
C
      RETURN
C------------  End of Subroutine CINIT  --------------------------------
      END

      SUBROUTINE OUTWEB (T, C, NS, MX, MY, LUN)
C-----------------------------------------------------------------------
C This routine prints the values of the individual species densities
C at the current time T, to logical unit LUN.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(NS,MX,MY)
C
      WRITE(LUN,10) T
 10   FORMAT(/1X,79('-')/30X,'At time t = ',E16.8/1X,79('-') )
C
      DO 40 I = 1,NS
        WRITE(LUN,20) I
 20     FORMAT(' the species c(',i2,') values are ')
        DO 30 JY = MY,1,-1
          WRITE(LUN,25) (C(I,JX,JY),JX=1,MX)
 25       FORMAT(6(1X,G12.6))
 30       CONTINUE
        WRITE(LUN,35)
 35     FORMAT(1X,79('-'),/)
 40     CONTINUE
C
      RETURN
C------------  End of Subroutine OUTWEB  -------------------------------
      END

      SUBROUTINE RESWEB (
     *     T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR, senpar)
C-----------------------------------------------------------------------
C This routine computes the residual vector, using Subroutine FWEB
C for the right-hand sides.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*),UPRIME(*),DELTA(*),RPAR(*),IPAR(*), senpar(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      ALPH = senpar(1)
      BETA = senpar(2) 
      CALL FWEB (T, U, DELTA, RPAR, senpar)
C
      DO 30 JY = 1,MY
         IYOFF = MXNS*(JY-1)
         DO 20 JX = 1,MX
            IC0 = IYOFF + NS*(JX-1)
            DO 10 I = 1,NS
               ICI = IC0 + I
               IF (I .GT. NP) THEN
                  DELTA(ICI) = -DELTA(ICI)
               ELSE
                  DELTA(ICI) = UPRIME(ICI) - DELTA(ICI)
               ENDIF
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
C
      RETURN
C------------  End of Subroutine RESWEB  -------------------------------
      END

      SUBROUTINE FWEB (T, CC, CRATE, RPAR, senpar)
C-----------------------------------------------------------------------
C This routine computes the right-hand sides of all the equations
C and returns them in the array CRATE.
C The interaction rates are computed by calls to WEBR, and these are
C saved in RPAR(1),...,RPAR(NEQ) for use later in preconditioning.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(*), CRATE(*), RPAR(*), senpar(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY,  FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      DO 60 JY = 1,MY
        IYOFF = MXNS*(JY-1)
        IDYU = MXNS
        IF (JY .EQ. MY) IDYU = -MXNS
        IDYL = MXNS
        IF (JY .EQ. 1) IDYL = -MXNS
        DO 40 JX = 1,MX
          IC = IYOFF + NS*(JX-1) + 1
C Get interaction rates at one point (X,Y).
          CALL WEBR (T, JX, JY, CC(IC), RPAR(2+IC),senpar)
          IDXU = NS
          IF (JX .EQ. MX) IDXU = -NS
          IDXL = NS
          IF (JX .EQ. 1) IDXL = -NS
          DO 20 I = 1,NS
            ICI = IC + I - 1
C Do differencing in Y.
            DCYLI = CC(ICI) - CC(ICI-IDYL)
            DCYUI = CC(ICI+IDYU) - CC(ICI)
C Do differencing in X.
            DCXLI = CC(ICI) - CC(ICI-IDXL)
            DCXUI = CC(ICI+IDXU) - CC(ICI)
C Collect terms and load CRATE elements.
            CRATE(ICI) = COY(I)*(DCYUI - DCYLI) + COX(I)*(DCXUI - DCXLI)
     1                  + RPAR(2+ICI)
 20         CONTINUE
 40       CONTINUE
 60    CONTINUE
      RETURN
C------------  End of Subroutine FWEB  ---------------------------------
      END

      SUBROUTINE WEBR (T, JX, JY, C, CRATE, senpar)
C-----------------------------------------------------------------------
C This routine computes one block of the interaction term R of the 
C system, namely block (JX,JY), for use in preconditioning.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(*), CRATE(*), senpar(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY,  FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      alph = senpar(1)
      beta = senpar(2)
      Y = REAL(JY-1)*DY
      X = REAL(JX-1)*DX
      DO 10 I = 1,NS
 10     CRATE(I) = 0.0D0
      DO 15 J = 1,NS
        CALL DAXPY (NS, C(J), ACOEF(1,J), 1, CRATE, 1)
 15     CONTINUE
      FAC = 1.0D0 + ALPH*X*Y + BETA*SIN(FPI*X)*SIN(FPI*Y)
      DO 20 I = 1,NS
 20     CRATE(I) = C(I)*(BCOEF(I)*FAC + CRATE(I))
      RETURN
C------------  End of Subroutine WEBR  ---------------------------------
      END
