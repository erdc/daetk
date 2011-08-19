      SUBROUTINE RES (
     *     T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR, senpara)
C-----------------------------------------------------------------------
C This routine computes the residual vector, using Subroutine FWEB
C for the right-hand sides.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*),UPRIME(*),DELTA(*),RPAR(*),IPAR(*),senpara(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      CALL FWEB (T, U, DELTA, RPAR,senpara)
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
      SUBROUTINE FWEB (T, CC, CRATE, RPAR, senpara)
C-----------------------------------------------------------------------
C This routine computes the right-hand sides of all the equations
C and returns them in the array CRATE.
C The interaction rates are computed by calls to WEBR, and these are
C saved in RPAR(1),...,RPAR(NEQ) for use later in preconditioning.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(*), CRATE(*), RPAR(*), senpara(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
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
          CALL WEBR (T, JX, JY, CC(IC), RPAR(2+IC), senpara)
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

      SUBROUTINE WEBR (T, JX, JY, C, CRATE, senpara)
C-----------------------------------------------------------------------
C This routine computes one block of the interaction term R of the 
C system, namely block (JX,JY), for use in preconditioning.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION senpara(*)
      DIMENSION C(*), CRATE(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      ALPH = senpara(1)
      BETA = senpara(2)
c      ALPH = 50.0D0
c      BETA = 100.0D0
      Y = REAL(JY-1)*DY
      X = REAL(JX-1)*DX
      DO 10 I = 1,NS
 10     CRATE(I) = 0.0D0
      DO 15 J = 1,NS
        CALL DAXPY (NS, C(J), ACOEF(1,J), 1, CRATE, 1)
 15     CONTINUE
      FAC = 1.0D0 + ALPH*X*Y + BETA*SIN(FPI*X)*SIN(FPI*Y)
      DO 20 I = 1,NS
         CRATE(I) = C(I)*(BCOEF(I)*FAC + CRATE(I))
c         CRATE(I) = 100.d0*(BCOEF(I)*FAC + CRATE(I))
 20      continue
      RETURN
C------------  End of Subroutine WEBR  ---------------------------------
      END
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
