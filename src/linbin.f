cccccccccc FORTRAN subroutine linbin.f cccccccccc

c Obtains bin counts for univariate data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 27/01/95

      subroutine linbin(X,n,a,b,M,trun,gcounts)     
      double precision X(*),a,b,gcounts(*),lxi,delta,rem
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
            gcounts(li) = gcounts(li) + (1-rem)
            gcounts(li+1) = gcounts(li+1) + rem
         elseif (li.lt.1.and.trun.eq.0) then
            gcounts(1) = gcounts(1) + 1
         elseif (li.ge.M) then
            if (li.eq.M.or.trun.eq.0) then 
               gcounts(M) = gcounts(M) + 1
            endif
         endif
20    continue

      return
      end

cccccccccc End of linbin.f cccccccccc
