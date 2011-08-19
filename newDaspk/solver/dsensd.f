      SUBROUTINE DSENSD(
     *     QRES, NEQ, T, Y, YPRIME, QSEN, INFO, 
     *     RWORK, IWORK, IRES, RPAR, IPAR, SENPAR)
C
C***BEGIN PROLOGUE DSENSD
C***DATE WRITTEN   980813   (YYMMDD)
C***REVISION DATE  981124   (Compatible with DASPK3.0)
C
C***PURPOSE  This subroutine is auxiliary to DASPK and can
C            be called directly by the user between calls to
C            DASPK.  Its purpose is to compute the sensitivity
C            of a derived quantity specified by the user
C            in the external subroutine QRES with respect
C            to the problem parameters, given the solution and
C            sensitivities output by DASPK at time T. 
C            This version is in double precision.
C            
C-------------------------------------------------------------------
C *Arguments:
C
C  QRES:EXT       This is the name of a subroutine which you
C                 provide to define the  derived quantity Q(t,y,y',p).
C                 The subroutine you provide must have the form
C                 SUBROUTINE QRES(T, Y, YPRIME, Q, IRES, RPAR, IPAR, SENPAR)
C                 to define the derived quantity
C                 Q = Q(T, Y, YPRIME, P)
C                 If INFO(20)=3 or 4 (ADIFOR is selected), QRES is
C                 the name of ADIFOR-generated routine.
C
C  T:IN           This is the current value of the independent variable
C                 on output from DASPK.
C
C  Y(*):IN        This array contains the solution components and
C                 sensitivities at T, on output from DASPK.
C
C  YPRIME(*):IN   This array contains contains the solution and
C                 sensitivity derivatives at T, on output from DASPK.
C            
C  QSEN(*):OUT    This array of dimension NQ(Np+1), where Np is the
C                 number of parameters (defined in INFO(19) in DASPK),
C                 will contain the values of Q and sensitivities of the derived
C                 quantities with respect to the problem parameters.
C                 
C
C  INFO(*):IN     This is the vector of code options that you provided
C                 to DASPK.  Note that DSENSD will follow the options
C                 that you specified for DASPK.  In particular,
C                 if INFO(20) = 0, sensitivities of the derived quantities
C                 will be computed by second order central differences,
C                 if INFO(20) = 1, they will be computed by first order
C                 differences.  The
C                 constant cnst specified in INFO(21) will be used
C                 to define the finite difference increment.
C                 If INFO(20) = 2, you should not be using this routine
C                 and should compute the sensitivities
C                 of the derived quantities yourself, given T, Y, YPRIME
C                 from DASPK.  To use this routine successfully,
C                 you must set INFO(24)=NQ as input to DASPK.
C                 If INFO(20) = 3, sensitivities of the derived quantities
C                 will be computed by ADIFOR with seed matrix option.
C                 If INFO(20) = 4, sensitivities of the derived quantities
C                 will be computed by ADIFOR with matrix-vector product only.
C
C  RWORK(*):IN    This is the real work array on output from DASPK.
C                 Note that if you are using this routine (DSENSD)
C                 to compute sensitivities of derived quantities, 
C                 you will need to augment the length of the RWORK
C                 array which is specified in the documentation to
C                 DASPK by an additional 2*NQ locations if
C                 INFO(20) = 0 or 1.
C                                   
C  IWORK(*):IN    This is the integer work array on output from DASPK.
C
C  RPAR,IPAR:IN   These are real and integer parameter arrays which
C                 you can use for communication between your calling
C                 program and the RES and QRES
C
C  SENPAR(*):IN   This is the real work array which is used in DASPK to
C                 store problem parameters of the sensitivities.
C-------------------------------------------------------------------------
C
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       EXTERNAL QRES
       DIMENSION Y(*),YPRIME(*),QSEN(*)
       DIMENSION INFO(*),RWORK(*),IWORK(*),RPAR(*),IPAR(*),SENPAR(*)
C
C      Set pointers into IWORK as in DASPK.
C
      PARAMETER (LML=1, LMU=2, LMTYPE=4, 
     *   LIWM=1, LMXORD=3, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,
     *   LNS=9, LNSTL=10, LNST=11, LNRE=12, LNJE=13, LETF=14, LNCFN=15,
     *   LNCFL=16, LNIW=17, LNRW=18, LNNI=19, LNLI=20, LNPS=21,
     *   LNSE=22, LMITER=23, LMAXL=24, LKMP=25, LNRMAX=26, LLNWP=27,
     *   LLNIWP=28, LLOCWP=29, LLCIWP=30, LKPRIN=31,
     *   LMXNIT=32, LMXNJ=33, LMXNH=34, LLSOFF=35, LNPD = 36, 
     *   LNY=37, LNRPD=38, LICNS=41)
C
C     Set pointers into RWORK as in DASPK.
C
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4, LCJ=5, LCJOLD=6,
     *   LHOLD=7, LS=8, LROUND=9, LEPLI=10, LSQRN=11, LRSQRN=12,
     *   LEPCON=13, LSTOL=14, LEPIN=15, LPRT=16,
     *   LALPHA=21, LBETA=27, LGAMMA=33, LPSI=39, LSIGMA=45, LDELTA=51)
C
C     Set up for multiprocessing
C
      MYID = INFO(26)
      IF (INFO(27) .EQ. 0) INFO(27) = 1 ! for the multiprocessing
      NUMPROCS = INFO(27)
      NAVG = INFO(19)/NUMPROCS
      IF (MYID .LT. MOD(INFO(19), NUMPROCS)) THEN
         NP = NAVG + 1
      ELSE
         NP = NAVG
      END IF
      NY = IWORK(LNY)
      MYNEQ = IWORK(LNY)*(NP+1)
C
C     Set pointers to RWORK and IWORK segments as in DASPK
C
C
      LSAVR = LDELTA
      IF (INFO(12) .NE. 0) LSAVR = LDELTA + MYNEQ
      LE  = LSAVR + MYNEQ
      LWT = LE + MYNEQ
      LVT = LWT
      IF (INFO(16) .NE. 0) LVT = LWT + MYNEQ
      LPHI = LVT + MYNEQ
      IF (INFO(11) .EQ. 0) THEN
         LWM  = LPHI + (IWORK(LMXORD)+1)*MYNEQ
      ELSE
         IF (INFO(12) .EQ. 0) THEN
            LWM = LPHI + MAX0(IWORK(LMXORD)+1,4)*MYNEQ
         ELSE
            LWM = LPHI + MAX0(IWORK(LMXORD)+1,5)*MYNEQ
         END IF
      END IF
C      
      LSE = LWM + IWORK(LNPD)
      LQ = LSE + 4*NY
c     If INFO(20) = 2, print out error message
C
C     Find the dimension of the derived quantities
      NQ = INFO(24)
C
C     sensitivity evaluated by ADIFOR
C
      IF (INFO(20).EQ.3) THEN    ! seed matrix option
         CALL drAdfSM(
     *        T,Y,YPRIME,QSEN,IRES,RPAR,IPAR,SENPAR,QRES,NY,INFO(19),
     *        IWORK(38),RWORK(LSE))
         RETURN
      ENDIF
      IF(INFO(20).EQ.4) THEN    ! matrix-vector product option
         CALL drAdfMV(
     *        T,Y,YPRIME,QSEN,IRES,RPAR,IPAR,SENPAR,QRES,NY,INFO(19),
     *        IWORK(38),RWORK(LSE))
         RETURN
      ENDIF
      IF (INFO(20).EQ.5.OR.INFO(20).EQ.2) THEN ! matrix times vector methods
         write(*,*) 'In DSENSD: Warning!! '
         write(*,*) ' INFO(20) = 2 OR 5 is not defined for SENSD. '
         IRES = -1
         return
      END IF
C
C
      CNST = RWORK(LPRT)
      IF(INFO(20) .EQ. 0) GO TO 100
C
C     First order forward finite difference scheme.
C
C     Initial call to user defined QRES routine.
      CALL QRES(T,Y,YPRIME,QSEN,IRES,RPAR,IPAR,SENPAR)
C
C
C     Iterate on the number of parameters, NP.
C
      DO IP = 1, NP
         MYIP = MYID+1 + (IP-1)*NUMPROCS
         ISROW = IP*NY
         IDROW = IP*NQ
C
C       Determine the perturbation.
C
         VNORM = 0.0D0
         DO J = 1,NY
            VNORM = VNORM + (RWORK(LWT+J-1)/RWORK(LWT+ISROW+J-1))**2
         ENDDO
         VNORM = DSQRT(VNORM)
         IF (MYIP .GT. INFO(22)) THEN
            DEL = CNST/VNORM
         ELSE
            DEL = CNST*MAX(ABS(SENPAR(MYIP)),1.0/VNORM )
         END IF
C
C       Save and perturb Y, YP and the parameter.
C
         IF (MYIP .LE. INFO(22)) THEN
            PSAVE = SENPAR(MYIP)
            SENPAR(MYIP) = PSAVE + DEL
         END IF
C
         DO I = 1,NY
            MRKR = ISROW + I
            RWORK(LSE+I-1) = Y(I)
            Y(I) =  Y(I) + DEL*Y(MRKR)
            RWORK(LSE+NY+I-1) = YPRIME(I)
            YPRIME(I) = YPRIME(I) + DEL*YPRIME(MRKR)
         ENDDO
C
C     Call RES with perturbed values.
C
         CALL QRES(T,Y,YPRIME,RWORK(LQ),IRES,RPAR,IPAR,SENPAR)
C
C     Approx. first order sensitivity, and restore Y, YP, parameter values.
C
         DELINV = 1.0D0/DEL
         DO I = 1,NQ
            QSEN(IDROW+I)=DELINV*(RWORK(LQ+I-1) - QSEN(I))
         ENDDO 
C
C       Restore y, yprime and parameter values.  Proceed to
C       next Ny sensitivity segment.
         DO I = 1,NY
            Y(I) = RWORK(LSE+I-1)
            YPRIME(I) =  RWORK(LSE+NY+I-1)
         ENDDO
         IF (MYIP .LE. INFO(22)) SENPAR(MYIP) = PSAVE
      ENDDO
C
C     First order forward finite difference approximation complete.
C
      RETURN
C
C     Second order centered finite difference scheme.
C
 100  CONTINUE
C
C     Initial call to user defined QRES routine.
      CALL QRES(T,Y,YPRIME,QSEN,IRES,RPAR,IPAR,SENPAR)
C
C     Iterate on the number of parameters, NP.
C
      DO IP = 1, NP
         MYIP = MYID+1 + (IP-1)*NUMPROCS
         ISROW = IP*NY
         IDROW = IP*NQ
C
C       Determine the perturbation.
C
         VNORM = 0.0D0
         DO J = 1,NY
            VNORM = VNORM + (RWORK(LWT+J-1)/RWORK(LWT+ISROW+J-1))**2
         ENDDO
         VNORM = DSQRT(VNORM)
         IF (MYIP .GT. INFO(22)) THEN
            DEL = CNST/VNORM
         ELSE
            DEL = CNST*MAX(ABS(SENPAR(MYIP)),1.0/VNORM )
         END IF
C
C       Save and perturb Y, YP and the parameter.
C
         IF (MYIP .LE. INFO(22)) THEN
            PSAVE = SENPAR(MYIP)
            SENPAR(MYIP) = PSAVE + DEL
         END IF
C
         DO I = 1,NY
            MRKR = ISROW + I
            RWORK(LSE+I-1) = Y(I)
            Y(I) =  Y(I) + DEL*Y(MRKR)
            RWORK(LSE+NY+I-1) = YPRIME(I)
            YPRIME(I) = YPRIME(I) + DEL*YPRIME(MRKR)
         ENDDO
C
C     Call RES with perturbed values.
C
         CALL QRES(T,Y,YPRIME,RWORK(LQ),IRES,RPAR,IPAR,SENPAR)
C
C       Save and perturb y, yprime and the parameter in the backward
C       (-) direction.
C      
         IF (MYIP .LE. INFO(22))  SENPAR(MYIP) = PSAVE - DEL
         DO I = 1,NY
            MRKR = ISROW + I
            Y(I) =  RWORK(LSE+I-1) - DEL*Y(MRKR)
            RWORK(LSE+NY+I-1) = YPRIME(I)
            YPRIME(I) = RWORK(LSE+NY+I-1) - DEL*YPRIME(MRKR)
         ENDDO
C
C     Call RES with perturbed values.
C
         CALL QRES(T,Y,YPRIME,RWORK(LQ+NQ),IRES,RPAR,IPAR,SENPAR)
C
C       Compute centered difference approximation to the sensitivity
C
         DELINV = 1.0D0/(2.0D0*DEL)
         DO I = 1,NQ
            QSEN(IDROW+I)=DELINV*(RWORK(LQ+I-1)-RWORK(LQ+NQ+I-1))
         ENDDO 
C
C       Restore y, yprime and parameter values.  Proceed to
C       next Ny sensitivity segment.
         DO I = 1,NY
            Y(I) = RWORK(LSE+I-1)
            YPRIME(I) =  RWORK(LSE+NY+I-1)
         ENDDO
         IF (MYIP .LE. INFO(22)) SENPAR(MYIP) = PSAVE
      ENDDO
C
C     Second order centered finite difference approximation complete.
C
      RETURN
      END
C-------------------------END OF DSEND ROUTINE----------------------
      subroutine drAdfMV(
     *     T,Y,YP,DELTA,IRES,RPAR,IPAR,senpar,G_RES,NY,isenfo,
     *     nrpd, wrk)
c===================================================================
c  This routine computes the state variable and sensitivity residuals 
C  through the Adifor with Matrix-vector product form.
C     AD_SCALAR_GRADIENTs = true
c
c   How to use this routine:
C    1. put all the staffs related with 
C           SUBROUTINE RES (T,Y,YPRIME,DELTA,IRES,RPAR,IPAR,SENPAR)
C       in a file called "res.f"
C    2. create file "res.cmp" with one line as
C        res.f
C    3. create file "res4senmv.adf" as
C     AD_PROG = res.cmp
C     AD_TOP = res
C     AD_IVARS = y, yprime, senpar
C     AD_OVARS = delta
C     AD_PMAX = 1
C     AD_PREFIX = g
C     AD_SCALAR_GRADIENTS = true
C     AD_OUTPUT_DIR = resmv
C
C    4. run Adifor generate g_res(...) with
C       % Adifor AD_SCRIPT=resmv.adf
C    
C    5. Call this routine with proper parameters
C
C**** Copyright (C) 1998, Shengtai Li
C         
C
C===================================================================
      implicit none
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  cj                ! const coming from the DASPK
      real*8  delta(*)          ! residual 
      integer ires              ! error indicator 
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)        ! sensitivity parameter array
      integer ipar(*)           ! integer parameter of the problems
                                ! pass the index of the rpar for sensitivity
      external g_res            ! residual generated by the Adifor
      integer ny                ! number of equations for state variables
      integer isenfo(*)         ! information for sensitivity size>9
                                ! see routine ddsen(*) in DASPK for details 
      integer nrpd              ! size of array rpar(*)
      real*8  wrk(*)            ! work array, size >= nrpd + isenfo(1)
c === local variable
      integer i, npp1, myid, numprocs, myi
      do i = 1, nrpd+isenfo(4)
         wrk(i) = 0.0d0
      end do
      myid = isenfo(8)
      numprocs = isenfo(9)
      if (nrpd .eq. 0) then
         do i = 1, isenfo(1)
            myi = myid+1 + (i-1)*numprocs
            if (myi .le. isenfo(4)) wrk(myi) = 1.0d0
            call g_res(t,y,y(1+i*Ny),yp,yp(1+i*Ny),
     1           delta,delta(1+i*Ny),ires,rpar,ipar,senpar,wrk)
            if (myi .le. isenfo(4)) wrk(myi) = 0.0d0
         end do
      else 
         npp1 = isenfo(4) + 1
         do i = 1, isenfo(1)
            myi = myid+1 + (i-1)*numprocs
            if (myi .le. isenfo(4)) wrk(myi) = 1.0d0
            call g_res(t,y,y(1+i*Ny),yp,yp(1+i*Ny),
     1           delta,delta(1+i*Ny),ires,rpar,wrk(npp1),ipar,
     2           senpar,wrk)
            if (myi .le. isenfo(4)) wrk(myi) = 0.0d0
         end do
      end if
      return
      end
c
      subroutine drAdfSM(
     *     T,Y,YP,DELTA,IRES,RPAR,IPAR,senpar,G_RES,NY,isenfo,
     *     nrpd, wrk)
c===================================================================
c  This routine computes the state variable and sensitivity residuals 
C  through the Adifor with seed matrix input.
C     AD_SCALAR_GRADIENT = false
c
c   How to use this routine:
C    1. put all the staffs related with 
C           SUBROUTINE RES (T,Y,YPRIME,DELTA,RPAR,IPAR,SENPAR)
C       in a file called "res.f"
C    2. create file "res.cmp" with one line as
C        res.f
C    3. create file "res4sensm.adf" as
C     AD_PROG = res.cmp
C     AD_TOP = res
C     AD_IVARS = y, yprime, senpar
C     AD_OVARS = delta
C     AD_PREFIX = h
C     AD_PMAX   = # of sensitivities
C     AD_OUTPUT_DIR = res4sen
C
C    4. run Adifor generate h_res(...) with
C       % Adifor AD_SCRIPT=res4sensm.adf
C    
C    5. Call this routine with proper parameters
C
C**** Copyright (C) 1998, Shengtai Li
C         
C
C===================================================================
      implicit none
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  delta(*)          ! residual 
      integer ires              ! error inidicator
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)         ! sensitivity parameter array
      integer ipar(*)           ! integer parameter of the problems
                                ! pass the index of the rpar for sensitivity
      external g_res            ! residual generated by the Adifor
      integer ny                ! number of equations totally
      integer isenfo(*)         ! information for sensitivity size>9
                                ! see routine ddsen(*) in DASPK for details 
      integer nrpd              ! size of array rpar(*)
      real*8  wrk(*)            ! work array, 
                                ! size = isenfo(1)*(2*ny+nrpd+isenfo(4)+nq)
c === local variable
      integer i, ig_rpar, np, nq, ig_y, ig_yp,ig_delta,
     *     ipos1, ip, ipos2, length, ig_senpar
      np = isenfo(1)
      nq = isenfo(6)
      length = np*ny
C === set Pointer for the work array
      ig_senpar = 0
      ig_y  = ig_senpar + isenfo(4)*np
      ig_yp = ig_y + length
      ig_rpar = ig_yp + length
      ig_delta = ig_rpar + nrpd
C
C === Set up the SEED Matrices for the call to RES
C
      do ip = 1, np
         ipos1 = ip*Ny
         do i=1,Ny
            ipos2 = ip + (i-1)*np
            wrk(ig_y + ipos2 ) = y(i+ ipos1 )
            wrk(ig_yp + ipos2 ) = yp(i+ ipos1 )
         end do
         do i = 1, isenfo(4)
            wrk(ip + (i-1)*np) = 0.0d0
         end do
C
C === Initialize the the vector that picks out the correct dF/dp(i)
C
         if (ip .le. isenfo(4)) wrk(ip+(ip-1)*np)=1.0d0
         do i = 1, nrpd 
            wrk(ig_rpar + ip + (i-1)*np) = 0.0d0
         end do
      end do
C
C Now call the ADIFOR generated  code 
c
      if (nrpd .eq. 0) then
         call g_res(np, t, y, wrk(ig_y+1), np, yp, wrk(ig_yp+1),
     *        np, delta, wrk(ig_delta+1), np, ires,
     *        rpar, ipar, senpar, wrk, np)
      else
         call g_res(np, t, y, wrk(ig_y+1), np, yp, wrk(ig_yp+1),
     *        np, delta, wrk(ig_delta+1), np, ires,
     *        rpar, wrk(ig_rpar+1), np, ipar, senpar, wrk, np)
      end if
c      
      do ip=1,np
         do i=1,nq
            delta(i+ip*Nq)=wrk(ig_delta + ip+(i-1)*Np)
         end do
      end do
c
      return
      end
c

