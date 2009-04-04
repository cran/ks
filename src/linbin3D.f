cccccccccc FORTRAN subroutine linbin3D.f cccccccccc

c Obtains bin counts for trivariate data
c via the linear binning strategy. In this version
c observations outside the mesh are ignored. 

c Last changed: 28 JUL 2005

      subroutine lbthrd(X,n,a1,a2,a3,b1,b2,b3,M1,M2,M3,gcounts)  
      integer n,M1,M2,M3,i,li1,li2,li3,ind1,ind2,ind3,ind4
      integer ind5,ind6,ind7,ind8
      double precision X(*),a1,a2,a3,b1,b2,b3,gcounts(*)
      double precision lxi1,lxi2,lxi3,delta1,delta2,delta3
      double precision rem1,rem2,rem3

c     Initialize grid counts to zero

      do 10 i = 1,(M1*M2*M3)
         gcounts(i) = dble(0)
10    continue

      delta1 = (b1 - a1)/(M1 - 1)
      delta2 = (b2 - a2)/(M2 - 1)
      delta3 = (b3 - a3)/(M3 - 1)
      do 20 i = 1,n
         lxi1 = ((X(i)-a1)/delta1) + 1
         lxi2 = ((X(n+i)-a2)/delta2) + 1
         lxi3 = ((X(2*n+i)-a3)/delta3) + 1

c        Find the integer part of "lxi1","lxi2" and "lxi3"

         li1 = lxi1
         li2 = lxi2
         li3 = lxi3
         rem1 = lxi1 - li1
         rem2 = lxi2 - li2 
         rem3 = lxi3 - li3
         if (li1.ge.1) then
            if (li1.lt.M1) then
               if (li2.ge.1) then
                  if (li2.lt.M2) then
                     if (li3.ge.1) then
                        if (li3.lt.M3) then
                           ind1 = li1+M1*(li2-1)+M1*M2*(li3-1)
                           ind2 = li1+1+M1*(li2-1)+M1*M2*(li3-1)
                           ind3 = li1+M1*li2+M1*M2*(li3-1)
                           ind4 = li1+1+M1*li2+M1*M2*(li3-1)
                           ind5 = li1+M1*(li2-1)+M1*M2*li3
                           ind6 = li1+1+M1*(li2-1)+M1*M2*li3
                           ind7 = li1+M1*li2+M1*M2*li3
                           ind8 = li1+1+M1*li2+M1*M2*li3
             gcounts(ind1) = gcounts(ind1)+(1-rem1)*(1-rem2)*(1-rem3)
             gcounts(ind2) = gcounts(ind2)+rem1*(1-rem2)*(1-rem3)
             gcounts(ind3) = gcounts(ind3)+(1-rem1)*rem2*(1-rem3)
             gcounts(ind4) = gcounts(ind4)+rem1*rem2*(1-rem3)
             gcounts(ind5) = gcounts(ind5)+(1-rem1)*(1-rem2)*rem3
             gcounts(ind6) = gcounts(ind6)+rem1*(1-rem2)*rem3
             gcounts(ind7) = gcounts(ind7)+(1-rem1)*rem2*rem3
             gcounts(ind8) = gcounts(ind8)+rem1*rem2*rem3                
                        endif
                     endif
                  endif
               endif  
            endif 
         endif
20    continue

      return
      end

cccccccccc End of linbin3D.f cccccccccc


