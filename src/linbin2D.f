cccccccccc FORTRAN subroutine linbin2D.f cccccccccc

c Obtains bin counts for bivariate data
c via the linear binning strategy. In this version
c observations outside the mesh are ignored. 

c Last changed: 25 AUG 1995

      subroutine lbtwod(X,n,a1,a2,b1,b2,M1,M2,gcounts)  
      integer n,M1,M2,i,li1,li2,ind1,ind2,ind3,ind4
      double precision X(*),a1,a2,b1,b2,gcounts(*)
      double precision lxi1,lxi2,delta1,delta2,rem1,rem2

c     Initialize grid counts to zero

      do 10 i = 1,(M1*M2)
         gcounts(i) = dble(0)
10    continue

      delta1 = (b1 - a1)/(M1 - 1)
      delta2 = (b2 - a2)/(M2 - 1)
      do 20 i = 1,n
         lxi1 = ((X(i) - a1)/delta1) + 1
         lxi2 = ((X(n+i) - a2)/delta2) + 1

c        Find the integer part of "lxi1" and "lxi2"

         li1 = lxi1
         li2 = lxi2
         rem1 = lxi1 - li1
         rem2 = lxi2 - li2 
         if (li1.ge.1) then              
            if (li2.ge.1) then
               if (li1.lt.M1) then
                  if (li2.lt.M2) then
                     ind1 = M1*(li2-1) + li1
                     ind2 = M1*li2 + li1
                     ind3 = M1*(li2-1) + li1 + 1
                     ind4 = M1*li2 + li1 + 1
                     gcounts(ind1) = gcounts(ind1)+(1-rem1)*(1-rem2)
                     gcounts(ind2) = gcounts(ind2)+rem1*(1-rem2)
                     gcounts(ind3) = gcounts(ind3)+(1-rem1)*rem2
                     gcounts(ind4) = gcounts(ind4)+rem1*rem2
                  endif
               endif  
            endif 
         endif
20    continue

      return
      end

cccccccccc End of linbin2D.f cccccccccc
