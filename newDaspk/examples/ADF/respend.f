       SUBROUTINE RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,senpar)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
       DIMENSION Y(*), YPRIME(*), DELTA(*),RPAR(*),IPAR(*),senpar(*)
*
          DELTA(1) = YPRIME(1) - Y(3)
          DELTA(2) = YPRIME(2) - Y(4)
          DELTA(3) = YPRIME(3) + Y(1)*Y(5)
          DELTA(4) = YPRIME(4) + Y(2)*Y(5) + 1.0D0
          DELTA(5) = cj*(Y(1)*Y(3) + Y(2)*Y(4)) 
*          delta(5) = y(3)*y(3) + y(4)*y(4) - y(5) - y(2)
*
       RETURN
       END 
