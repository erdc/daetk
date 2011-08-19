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
c       cu1 = 421.8
c       cu2 = -147.8

       delta(1) = yprime(1) - y(5) 
       delta(2) = yprime(2) - y(6)
       delta(3) = yprime(3) - y(7)
       delta(4) = yprime(4) - y(8)
       delta(5) = yprime(5) + 1.0/cm2*y(10)*dsin(y(9))
       delta(6) = yprime(6) + 1.0/cm2*(y(10)*dcos(y(9)) - cm*cg)
       delta(7) = yprime(7) + 1.0/ccj*(cc2*y(7)+cc3*cu2-cc3*cc3*y(10))
       delta(8) = yprime(8) + 1.0/cm1*(cc1*y(8) - cu1 -y(10)*sin(y(9)))
       delta(9) = y(9) - atan((y(1)-y(4))/y(2))
       a1 = (y(1)-y(4))*(y(5)-y(8))
       delta(10)= cj*(y(3)*y(7)-a1-y(2)*y(6))
       RETURN
       END 
