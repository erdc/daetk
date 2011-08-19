       program crane
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
       EXTERNAL RES, i_res, h_res
       DIMENSION Y(30),RWORK(1000),IWORK(200),YPRIME(30) 
       DIMENSION INFO(30),RPAR(1),IPAR(1), senpar(2), g(3)
* 
* Specify Error tolerances, length of iwork and rwork, NEQ and the
* interval of integeration.
*
       RTOL = 1.0d-4
       ATOL = 1.0d-4
       LRW = 1000
       LIW = 200
       T = 0.0
       TOUT = 1.0
       neq = 10
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
* Jacobian evaluation (0,1,2,3)
*
       info(5) = 2
C
C Here set INFO(14) = 1 to get the computed initial values.
       INFO(14) = 1
*
* Set x5 as a algebraic variable, does not do the error test on it.
*
       info(16) = 1
       NEQ = 10
       do i = 1, 8
          iwork(40+i) = 1
          iwork(40+neq+i) = 0
       end do
c       iwork(41) = 2
c       iwork(42) = 2
c       iwork(45) = 2
c       iwork(46) = 2
c       iwork(48) = 2

       iwork(49) = -1           !index-1
       iwork(50) = -1           !index-2
       iwork(60) = 1
       iwork(38) = 0
*
* Make the IC's consistent and relay this info. to DASSLSO.
*
*      call consistic(y,yprime)
      INFO(11) = 1
      PRINT *, 'initial condition options(0,1,3,4,5)'
      read(*,*) info(11)
      inso = 18
      info(inso+1) = 1
      print *, ' number of sensitivities? (0,1,2) '
      read(*,*) info(inso+1)
      Ny = neq
      neq = neq*(info(inso+1) + 1)
      print *, ' method for computing sensitivity residues:(0,1,3)'
      print *, '    0: central differencing '
      print *, '    1: one-side forward differencing '
      print *, '    3: Adifor input'
      read (*,*) info(inso+2)
      info(inso+3) = 1          ! default perturbed factor
      info(inso+4) = 2          ! number of sensitivity parameter in the RES
      info(inso+5) = 0          ! 0, full control 1. partial control
      info(inso+6) = 0          ! no derived information
      info(inso+7) = 1          ! simultaneous corrector(0), staggered(1)

      do i = 1,neq
         print*, 'y(', i, ')=', y(i), ', yprime(i)= ',yprime(i)
      enddo
      g(1) = y(3)*y(3)-(y(1)-y(4))**2-y(2)**2
      g(2) = y(3)*y(7)-(y(1)-y(4))*(y(5)-y(8))-y(2)*y(6)
      g(3) = y(3)*yprime(7)+yprime(3)*y(7)-(y(1)-y(4))*
     *     (yprime(5)-yprime(8)) - (yprime(1)-yprime(4))*(y(5)-y(8))
     *     -y(2)*yprime(6)-yprime(2)*y(6)
      print *, ' g(1) =', g(1)
      print *, ' g(2) =', g(2)
      print *, ' g(3) =', g(3)
      print*, 'input the TOUT:'
      read*, TOUT
*
 20   IF(T .LT. TOUT.and.idid.ge.0) THEN 
*
         CALL DDASPK(
     *        RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *        LRW,IWORK,LIW,RPAR,IPAR,i_res, psol,senpar, H_res)
         print*, 'info(11) = ', info(11), ', idid = ', idid
         if (info(11).ne.0)  then
            info(11) = 0
            tout = 1.0
            write(*,*)'time = ',T  
            do i = 1,neq
               print*, 'y(', i, ')=', y(i), ',yprime(i)=',yprime(i)
            enddo
         g(1) = y(3)*y(3)-(y(1)-y(4))**2-y(2)**2
         g(2) = y(3)*y(7)-(y(1)-y(4))*(y(5)-y(8))-y(2)*y(6)
         g(3) = y(3)*yprime(7)+yprime(3)*y(7)-(y(1)-y(4))*
     *     (yprime(5)-yprime(8)) - (yprime(1)-yprime(4))*(y(5)-y(8))
     *     -y(2)*yprime(6)-yprime(2)*y(6)
         print *, ' g(1) =', g(1)
         print *, ' g(2) =', g(2)
         print *, ' g(3) =', g(3)
            pause 
         end if
*     
         write(20,*) y(1), y(2), y(3)
         GOTO 20
      ENDIF 
      PRINT*
*
* Output IWORK and RWORK vector values
*
        WRITE(*,*)'DASSLSO SOLUTION FOR THE SINGLE '
        write(*,*)'PENDULUM PROBLEM REQUIRED:'
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
       cg = 9.81
       cm1 = 3000.0
       cm2 = 100.0
       ccj = 6.24
       cm = 100.0
       cc1 = 0.2
       cc2 = 0.2*0.15
       cc3 = 0.15
       cu1 = senpar(1)
       cu2 = senpar(2)

       delta(1) = yprime(1) - y(5) 
       delta(2) = yprime(2) - y(6)
       delta(3) = yprime(3) - y(7)
       delta(4) = yprime(4) - y(8)
       delta(5) = yprime(5) + 1.0/cm2*y(10)*dsin(y(9))
       delta(6) = yprime(6) + 1.0/cm2*(y(10)*dcos(y(9)) - cm*cg)
       delta(7) = yprime(7) + 1.0/ccj*(cc2*y(7)+cc3*cu2-cc3*cc3*y(10))
       delta(8) = yprime(8) + 1.0/cm1*(cc1*y(8) - cu1 -y(10)*sin(y(9)))
       delta(9) = y(9) - atan((y(1)-y(4))/y(2))
       delta(10)= cj*(y(3)*y(7)-(y(1)-y(4))*(y(5)-y(8))-y(2)*y(6))

       RETURN
       END 
*------------------------------------------------------------------------*
*------------------------------------------------------------------------*
       SUBROUTINE CONSISTIC(Y,YPRIME)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION Y(*), YPRIME(*)
*
* Set yprime to the value of the RHS, for consistency.
*
*
       cg = 9.81
       cm1 = 3000.0
       cm2 = 100.0
       ccj = 6.24
       cm = 100.0
       cc1 = 0.2
       cc2 = 0.2*0.15
       cc3 = 0.15
       cu1 = 421.8
       cu2 = -147.8

       yprime(1) = y(5) 
       yprime(2) = y(6)
       yprime(3) = y(7)
       yprime(4) = y(8)
       yprime(5) = - 1.0/cm2*y(10)*dsin(y(9))
       yprime(6) = - 1.0/cm2*(y(10)*dcos(y(9)) - cm*cg)
       yprime(7) = - 1.0/ccj*(cc2*y(7)+cc3*cu2-cc3*cc3*y(10))
       yprime(8) = - 1.0/cm1*(cc1*y(8) - cu1 -y(10)*sin(y(9)))
       return
       end
*------------------------------------------------------------------------* 
       SUBROUTINE INIT(Y,YPRIME,senpar)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
       DIMENSION Y(*), YPRIME(*), senPAR(*)
*
* Set y, yprime, IPAR, RPAR to 0, initially.
*
       do i = 1,30
         y(i) = 0.0
         yprime(i) = 0.0
       enddo
*
*
* Parameter values - placed somewhere in senPAR.
*
       senpar(1) = 421.8
       senpar(2) = -147.8
*
* Specify the initial values for the SP problem.
*
       y(1) = 3.0
       y(2) = 5.0
       y(3) = 5.0
       y(4) = 3.0
c       y(1) = 2.0
c       y(2) = 6.0
c       y(3) = 9.0
c       y(4) = 1.0
       y(5) = 1.0
       y(6) = 1.0
       y(7) = 2.0
       y(8) = 2.0


       y(9) = 1.0
       y(10) = 0.0
       return
       end
*------------------------------------------------------------------------* 


