cccccccccc FORTRAN subroutine linbin2009.f cccccccccc

c Obtains bin counts for univariate data via the linear 
c binning strategy. If "trun=0" then weight from end 
c observations is given to corresponding end grid points. 
c If "trun=1" then end observations are truncated.

c Last changed: 20 MAR 2009

      subroutine linbin(X,n,a,b,M,trun,w,gcounts)     
      double precision X(*),a,b,w(*),gcounts(*),lxi,delta,rem
      integer n,M,i,li,trun

c     Initialize grid counts to zero

      do 10 i=1,M
         gcounts(i) = dble(0)
10    continue

      delta = (b-a)/(M-1)
      do 20 i=1,n
         lxi = ((X(i)-a)/delta) + 1

c        Find integer part of "lxi"

         li = lxi 

         rem = lxi - li
         if (li.ge.1.and.li.lt.M) then
            gcounts(li) = gcounts(li) + w(i)*(1-rem)
            gcounts(li+1) = gcounts(li+1) + w(i)*rem
         endif

         if (li.lt.1.and.trun.eq.0) then
            gcounts(1) = gcounts(1) + w(i)
         endif

         if (li.ge.M.and.trun.eq.0) then
            gcounts(M) = gcounts(M) + w(i)
         endif
 20   continue

      return
      end

cccccccccc End of linbin2009.f cccccccccc
