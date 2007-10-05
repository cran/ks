cccccccccc FORTRAN subroutine linbin4D.f cccccccccc

c Obtains bin counts for quadrivariate data
c via the linear binning strategy. In this version
c observations outside the mesh are ignored. 

c Last changed: 31 AUG 2005

      subroutine lbfoud(X,n,a1,a2,a3,a4,b1,b2,b3,b4,M1,M2,M3,M4,gcounts)  
      integer n,M1,M2,M3,M4,i,li1,li2,li3,li4,ind1,ind2,ind3,ind4
      integer ind5,ind6,ind7,ind8,ind9,ind10,ind11,ind12,ind13,ind14
      integer ind15,ind16
      double precision X(*),a1,a2,a3,a4,b1,b2,b3,b4,gcounts(*)
      double precision lxi1,lxi2,lxi3,lxi4,delta1,delta2,delta3,delta4
      double precision rem1,rem2,rem3,rem4

c     Initialize grid counts to zero

      do 10 i = 1,(M1*M2*M3*M4)
         gcounts(i) = dble(0)
10    continue

      delta1 = (b1 - a1)/(M1 - 1)
      delta2 = (b2 - a2)/(M2 - 1)
      delta3 = (b3 - a3)/(M3 - 1)
      delta4 = (b4 - a4)/(M4 - 1)
      do 20 i = 1,n
         lxi1 = ((X(i)-a1)/delta1) + 1
         lxi2 = ((X(n+i)-a2)/delta2) + 1
         lxi3 = ((X(2*n+i)-a3)/delta3) + 1
         lxi4 = ((X(3*n+i)-a4)/delta4) + 1

c        Find the integer part of "lxi1","lxi2", "lxi3" and "lxi4"

         li1 = lxi1
         li2 = lxi2
         li3 = lxi3
         li4 = lxi4
         rem1 = lxi1 - li1
         rem2 = lxi2 - li2 
         rem3 = lxi3 - li3
         rem4 = lxi4 - li4

         if (li1.ge.1) then
            if (li1.lt.M1) then
               if (li2.ge.1) then
                  if (li2.lt.M2) then
                     if (li3.ge.1) then
                        if (li3.lt.M3) then
                           if (li4.gt.1) then
                              if(li4.lt.M4) then
       ind1 = li1+M1*(li2-1)+M1*M2*(li3-1)+M1*M2*M3*(li4-1)
       ind2 = li1+1+M1*(li2-1)+M1*M2*(li3-1)+M1*M2*M3*(li4-1)
       ind3 = li1+M1*li2+M1*M2*(li3-1)+M1*M2*M3*(li4-1)
       ind4 = li1+1+M1*li2+M1*M2*(li3-1)+M1*M2*M3*(li4-1)
       ind5 = li1+M1*(li2-1)+M1*M2*li3+M1*M2*M3*(li4-1)
       ind6 = li1+1+M1*(li2-1)+M1*M2*li3+M1*M2*M3*(li4-1)
       ind7 = li1+M1*li2+M1*M2*li3+M1*M2*M3*(li4-1)
       ind8 = li1+1+M1*li2+M1*M2*li3+M1*M2*M3*(li4-1)
       
       ind9 = li1+M1*(li2-1)+M1*M2*(li3-1)+M1*M2*M3*li4
       ind10 = li1+1+M1*(li2-1)+M1*M2*(li3-1)+M1*M2*M3*li4
       ind11 = li1+M1*li2+M1*M2*(li3-1)+M1*M2*M3*li4
       ind12 = li1+1+M1*li2+M1*M2*(li3-1)+M1*M2*M3*li4
       ind13 = li1+M1*(li2-1)+M1*M2*li3+M1*M2*M3*li4
       ind14 = li1+1+M1*(li2-1)+M1*M2*li3+M1*M2*M3*li4
       ind15 = li1+M1*li2+M1*M2*li3+M1*M2*M3*li4
       ind16 = li1+1+M1*li2+M1*M2*li3+M1*M2*M3*li4
       
       gcounts(ind1) = gcounts(ind1)+(1-rem1)*(1-rem2)*(1-rem3)*(1-rem4)
       gcounts(ind2) = gcounts(ind2)+rem1*(1-rem2)*(1-rem3)*(1-rem4)
       gcounts(ind3) = gcounts(ind3)+(1-rem1)*rem2*(1-rem3)*(1-rem4)
       gcounts(ind4) = gcounts(ind4)+rem1*rem2*(1-rem3)*(1-rem4)
       gcounts(ind5) = gcounts(ind5)+(1-rem1)*(1-rem2)*rem3*(1-rem4)
       gcounts(ind6) = gcounts(ind6)+rem1*(1-rem2)*rem3*(1-rem4)
       gcounts(ind7) = gcounts(ind7)+(1-rem1)*rem2*rem3*(1-rem4)
       gcounts(ind8) = gcounts(ind8)+rem1*rem2*rem3*(1-rem4)

       gcounts(ind9)  = gcounts(ind9)+(1-rem1)*(1-rem2)*(1-rem3)*rem4
       gcounts(ind10) = gcounts(ind10)+rem1*(1-rem2)*(1-rem3)*rem4
       gcounts(ind11) = gcounts(ind11)+(1-rem1)*rem2*(1-rem3)*rem4
       gcounts(ind12) = gcounts(ind12)+rem1*rem2*(1-rem3)*rem4
       gcounts(ind13) = gcounts(ind13)+(1-rem1)*(1-rem2)*rem3*rem4
       gcounts(ind14) = gcounts(ind14)+rem1*(1-rem2)*rem3*rem4
       gcounts(ind15) = gcounts(ind15)+(1-rem1)*rem2*rem3*rem4
       gcounts(ind16) = gcounts(ind16)+rem1*rem2*rem3*rem4
                              endif
                           endif 
                        endif
                     endif
                  endif
               endif  
            endif 
         endif
20    continue

      return
      end

cccccccccc End of linbin4D.f cccccccccc


