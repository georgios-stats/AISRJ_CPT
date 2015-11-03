module fixedupdates_mod

      implicit none

      private
      public      :: updateH,updateS
      public      :: updateBeta,updateAlpha

contains

      subroutine updateH(k,h,alpha,beta,s,y,en,jup)

            implicit none

            integer, intent(in)                 :: k
            double precision, intent(inout)     :: h(0:k)
            double precision, intent(in)        :: alpha
            double precision, intent(in)        :: beta
            double precision, intent(in)        :: s(0:k+1)
            integer, intent(in)                 :: en
            integer, intent(in)                 :: y(en)
            integer, intent(in)                 :: jup

            integer                 :: count_y
            integer                 :: i
            double precision        :: hNew
            double precision        :: u
            double precision        :: logAccRat

            ! PROPOSE THE CANDIDATE h
            call rnguniform(u)
            u = -0.5d0+u
            hNew = exp(u)*h(jup)

            ! COMPUTE THE LOG ACCEPTANCE RATIO
            count_y = 0
            do i = 1,en
                  if (y(i).ge.s(jup).and.y(i).lt.s(jup+1)) then
                        count_y = count_y+1
                  endif
            enddo

            logAccRat = alpha*(log(hNew)-log(h(jup))) &
                              -beta*(hNew-h(jup))      &
                              -(hNew-h(jup))*(s(jup+1)-s(jup)) &
                              +(log(hNew)-log(h(jup)))*count_y

            ! ACCEPT/REJECT
            call rnguniform(u)
            if (logAccRat.gt.log(u)) then
                  h(jup) = hNew
            endif

      end subroutine updateH

      subroutine updateS(k,h,alpha,beta,s,y,en,jup)

            implicit none

            integer, intent(in)                 :: k
            double precision, intent(inout)     :: s(0:k+1)
            double precision, intent(in)        :: alpha
            double precision, intent(in)        :: beta
            double precision, intent(in)        :: h(0:k)
            integer, intent(in)                 :: en
            integer, intent(in)                 :: y(en)
            integer, intent(in)                 :: jup

            double precision                    :: sNew
            double precision                    :: logAccRat
            double precision                    :: u
            integer                             :: i
            integer                             :: count_y1
            integer                             :: count_y2
            integer                             :: countNew_y1
            integer                             :: countNew_y2

            ! PROPOSE THE CANDIDATE s
            call rnguniform(u)
            sNew = s(jup-1) +(s(jup+1)-s(jup-1))*u

            ! COMPUTE THE LOG ACCEPTANCE RATIO
            count_y1 = 0
            count_y2 = 0
            countNew_y1 = 0
            countNew_y2 = 0
            do i = 1,en
                  if (y(i).ge.s(jup-1).and.y(i).lt.s(jup)) then
                        count_y1 = count_y1+1
                  endif
                  if (y(i).ge.s(jup-1).and.y(i).lt.sNew) then
                        countNew_y1 = countNew_y1+1
                  endif
                  if (y(i).ge.s(jup).and.y(i).lt.s(jup+1)) then
                        count_y2 = count_y2+1
                  endif
                  if (y(i).ge.sNew.and.y(i).lt.s(jup+1)) then
                        countNew_y2 = countNew_y2+1
                  endif
            enddo

            logAccRat = log(h(jup-1))*(countNew_y1-count_y1) &
                              +log(h(jup))*(countNew_y2-count_y2) &
                              +(h(jup)-h(jup-1))*(sNew-s(jup)) &
                              -log(s(jup)-s(jup-1))+log(sNew-s(jup-1)) &
                              +log(s(jup+1)-sNew)-log(s(jup+1)-s(jup))

            ! ACCEPT/REJECT
            call rnguniform(u)
            if (logAccRat.gt.log(u)) then
                  s(jup) = sNew
            endif

      end subroutine updateS

      subroutine updateBeta(k,h,alpha,beta,ee,ff)

            implicit none

            integer, intent(in)                 :: k
            double precision, intent(in)        :: h(0:k)
            double precision, intent(in)        :: alpha
            double precision, intent(in)        :: ee
            double precision, intent(in)        :: ff
            double precision, intent(inout)     :: beta

            double precision                    :: par1
            double precision                    :: par2

            par1 = (k+1.d0)*alpha+ee
            par2 = sum(h(0:k))+ff

            call rnggamma(beta,par1)
            beta = beta/par2

      end subroutine updateBeta

      subroutine updateAlpha(k,h,alpha,cc,dd,aspr,beta)

            implicit none

            integer, intent(in)                 :: k
            double precision, intent(in)        :: h(0:k)
            double precision, intent(in)        :: cc
            double precision, intent(in)        :: dd
            double precision, intent(in)        :: beta
            double precision, intent(in)        :: aspr
            double precision, intent(inout)     :: alpha

            double precision                    :: alphaNew
            double precision                    :: u
            double precision                    :: logAccRat
            double precision                    :: dlgama

            call rnguniform(u)

            alphaNew = alpha*aspr**(u-0.5d0)

            logAccRat = (alphaNew-alpha)* &
                                    ((k+1.d0)*log(beta)+sum(log(h(0:k)))-dd) &
                              -(k+1.d0)*(dlgama(alphaNew)-dlgama(alpha)) &
                              +cc*(log(alphaNew)-log(alpha))

            call rnguniform(u)

            if (logAccRat.gt.log(u)) then
                  alpha = alphaNew
            endif

      end subroutine updateAlpha

end module fixedupdates_mod


