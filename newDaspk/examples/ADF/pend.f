       program pend 
*  
*      *** Demo program for DASPK ***
*
* Single Pendulum problem with sensitivity equation coupling 
* run  for  0.0 <= t <= 1.0. 
*
* x1' = x3                      x1(0) = 0.5
* x2' = x4                      x2(0) = -sqrt(p1^2 - x1^2)
* x3' = -x1*x5                  x3(0) = 0.0
* x4' = -x2*x5 + g              x4(0) = 0.0
* 0   = x1*x3 + x2*x4
*
*
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
       EXTERNAL RES, i_res, h_res, j_res
       DIMENSION Y(10),RWORK(300),IWORK(100),YPRIME(10) 
       DIMENSION INFO(30),RPAR(1),IPAR(20), senpar(1)
* 
* Specify Error tolerances, length of iwork and rwork, NEQ and the
* interval of integeration.
*
       RTOL = 1.0d-6
       ATOL = 1.0d-6
       LRW = 300
       LIW = 100
       T = 0.0
       TOUT = 1.0
       neq = 5
*
* Initialize y, yprime, senpar
*
       call init(y,yprime,senpar)
*
* Initialize the INFO vector to 0.
*
       DO 10 I = 1,30
          INFO(I) = 0
  10   CONTINUE
*
* Get solution at intermediate steps.
*
       INFO(3)  = 1
*
* Finite difference approx. (FDA) of Jacobian. (1 = no,0 = yes)
*
       info(5) = 2
       print *, '   0 ------ finite difference '
       print *, '   2 ------ ADIFOR with SparcLinC option '
       print *, '   3 ------ ADIFOR with seed matrix option '
       print *, ' please input the Jacobian method: (0,2,3)'
       read(*,*)info(5)
C
C Here set INFO(14) = 1 to get the computed initial values.
       INFO(14) = 1
*
* Set x5 as a algebraic variable, does not do the error test on it.
*
       info(16) = 1
       NEQ = 5
       do i = 1, 4
          iwork(40+i) = 1
          iwork(40+5+i) = 0
       end do
       iwork(41) = 2
       iwork(42) = 2

       iwork(45) = -1
       iwork(50) = 1
       iwork(38) = 0
*
*
* Make the IC's consistent and relay this info. to DASSLSO.
*
c      call consistic(y,yprime)
      INFO(11) = 1
      PRINT *, 'initial condition options(3,4,5)'
      read(*,*) info(11)
      inso = 18
      info(inso+1) = 1
      print *, ' number of sensitivities? (0,1) '
      read(*,*) info(inso+1)
      Ny = neq
      neq = neq*(info(inso+1) + 1)
      print *, '    0: central differencing '
      print *, '    1: one-side forward differencing '
      print *, '    3: Adifor input'
      print *, ' method for computing sensitivity residuals:(0,1,3)'
      read (*,*) info(inso+2)
      info(inso+3) = 1          ! default perturbed factor
      info(inso+4) = 0          ! no sensitivity parameter in the RES
      info(inso+5) = 0          ! 0, full control 1. partial control
      info(inso+6) = 0          ! no derived information
      info(inso+7) = 1          ! simultaneous corrector(0), staggered(1)
      print *, '    0:  simultaneous corrector method '
      print *, '    1:  staggered corrector method'
      print *, ' Corrector method: (0, 1)'
      read(*,*) info(inso+7)

      print *, ' Initial values are:'
      do i = 1,neq
         print*, 'y(', i, ')=', y(i), ', yprime(i)= ',yprime(i)
      enddo
      print *, 'y(1)**2 + y(2)**2 - 1.0         = ', 
     *     y(1)*y(1) + y(2)*y(2) - 1.0
      print *, 'y(1)*y(3)+y(2)*y(4)             = ', 
     *     y(1)*y(3) + y(2)*y(4)
      print *, 'y(3)**2 + y(4)**2 - y(5) - y(2) = ',
     *     y(3)**2 + y(4)**2 - (y(1)*y(1)+y(2)*y(2))*y(5)-y(2)  
      print *, 'input the end time TOUT:'
      read *, TOUT
*     
 20   IF(T .LT. TOUT) THEN 
*     
         if (info(5) .eq. 2) then 
            CALL DDASPK(
     *           RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *           LRW,IWORK,LIW,RPAR,IPAR,i_res, psol,senpar, H_res)
         else 
            CALL DDASPK(
     *           RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *           LRW,IWORK,LIW,RPAR,IPAR,j_res, psol,senpar, H_res)
         end if

         if (info(11).ne.0)  then
            info(11) = 0
            write(*,99)t, y(1)*y(1) + y(2)*y(2) - 1.0, 
     *        y(1)*y(3) + y(2)*y(4),
     *        y(3)**2 + y(4)**2 - (y(1)*y(1)+y(2)*y(2))*y(5)-y(2)
            pause
         end if

         write(20,99)t, y(1), y(2), y(6), y(7)
         write(21,99)t, y(1)*y(1) + y(2)*y(2) - 1.0, 
     *        y(1)*y(3) + y(2)*y(4),
     *        y(3)**2 + y(4)**2 - (y(1)*y(1)+y(2)*y(2))*y(5)-y(2)
         if (idid.gt.0) GOTO 20
      ENDIF 
 99   format(1x, 5e15.6)
      PRINT *
*
* Output IWORK and RWORK vector values
*
      PRINT*
      WRITE(*,*)'TOTAL # OF STEPS:',IWORK(11)
      PRINT*
      WRITE(*,*)'TOTAL # OF ERROR TEST FAILURES:',IWORK(14)
      PRINT*
      WRITE(*,*)'TOTAL # OF CONV. TEST FAILURES:',IWORK(15)
      PRINT*
      WRITE(*,*)'TOTAL # OF RES function evals:',IWORK(12)
*     
      END
*------------------------------------------------------------------------*
      SUBROUTINE RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR, senpar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION Y(*), YPRIME(*), DELTA(*), RPAR(*), IPAR(*),senpar(*)
*     
      DELTA(1) = YPRIME(1) - Y(3)
      DELTA(2) = YPRIME(2) - Y(4)
      DELTA(3) = YPRIME(3) + Y(1)*Y(5)
      DELTA(4) = YPRIME(4) + Y(2)*Y(5) + 1.0D0
      DELTA(5) = cj*(Y(1)*Y(3) + Y(2)*Y(4)) 
*     print *,' delta(5) =', delta(5)
*     delta(5) = y(3)*y(3) + y(4)*y(4) - y(5) - y(2)
*     
      RETURN
      END 
*------------------------------------------------------------------------*
*------------------------------------------------------------------------*
      SUBROUTINE CONSISTIC(Y,YPRIME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*), YPRIME(*)
*     
*     Set yprime to the value of the RHS, for consistency.
*     
      YPRIME(1) = Y(3)
      YPRIME(2) = Y(4)
      YPRIME(3) = -Y(1)*Y(5)
      YPRIME(4) = -Y(2)*Y(5) - 1.0D0
*     
      return
      end
*------------------------------------------------------------------------* 
      SUBROUTINE INIT(Y,YPRIME,senpar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION Y(*), YPRIME(*), senPAR(*)
*     
*     Set y, yprime, IPAR, RPAR to 0, initially.
*     
      do i = 1,10
         y(i) = 0.0
         yprime(i) = 0.0
      enddo
*     
*     
*     Parameter values - placed somewhere in senPAR.
*     
      senpar(1) = 1.0d0
*     
*     Specify the initial values for the SP problem.
*     
      Y(1) = 5.0D-1
      temp =  sqrt(senpar(1)**2 - y(1)**2) 
      Y(2) = -temp
      y(3) = 10.0
      y(4) = 10.0
C     
      y(5) = 10.0
      y(7) = -senpar(1)/temp
c     ... the state variables that satisfy the constraints.
c     y(3) = 11.82969400622030
c     y(4) = 6.829877018922193
*     y(5) = y(3)*y(3) + y(4)*y(4) - y(2) 
*     
*     Y(8) = 2.113537329186674
*     y(9) = -7.886251345948634
*     
      return
      end
*------------------------------------------------------------------------* 


