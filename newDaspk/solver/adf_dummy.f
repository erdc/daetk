c
C === dummy routines with ADIFOR with SparsLinC option
C     used when you cannot access the SparsLinC library from ADIFOR
C
      subroutine xspini
      return
      end
      subroutine dspxsq(indvec,ajac,num,iwk,joutlen,ires)
      integer indvec, num, iwk, joutlen, ires
      real*8  ajac
      return
      end
      subroutine  dspsd(indvec, i, a, n)
      integer indvec, i, n
      real*8 a
      return
      end
      subroutine ehrpt
      return
      end
 
