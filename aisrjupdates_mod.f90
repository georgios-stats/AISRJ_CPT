

! Copyrigtht 2012 Georgios Karagiannis
!
! This file is part of AISRJ_CPT.
!
! AISRJ_CPT is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation version 2 of the License.
!
! AISRJ_CPT is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with AISRJ_CPT.  If not, see <http://www.gnu.org/licenses/>.

! ----------------------------------------------------------------------

! REFERENCES:
! -----------
!
! Karagiannis, G., & Andrieu, C. (2013).
! Annealed importance sampling reversible jump MCMC algorithms.
! Journal of Computational and Graphical Statistics, 22(3), 623-648.
!
! Georgios Karagiannis
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email : Georgios.Karagiannis@pnnl.gov
! Email (current): georgios-stats@gmail.com
!
! Christophe Andrieu
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email: C.Andrieu@bristol.ac.uk


! ----------------------------------------------------------------------

module aisrjupdates_mod

      implicit none

      private
      public      :: AnnImpWeight
      public      :: AisRjSplitSweep,gimel

contains

      ! ================================================================
      ! ANNEALING IMPORTANCE WEIGHTS
      ! ================================================================

      subroutine AnnImpWeight(logIW,k,lambda,h,alpha,beta,s,Leng,jstar,y,en)

            ! The probability 1/(k+1) of killing a change point is considered 
            ! here.

            implicit none

            integer, intent(in)                 :: k
            integer, intent(in)                 :: jstar
            double precision, intent(inout)     :: logIW
            double precision, intent(in)        :: s(0:k+2)
            double precision, intent(in)        :: h(0:k+1)
            double precision, intent(in)        :: Leng
            double precision, intent(in)        :: alpha
            double precision, intent(in)        :: beta
            double precision, intent(in)        :: lambda
            integer, intent(in)                 :: en
            integer, intent(in)                 :: y(en)

            integer                             :: I01
            integer                             :: I12
            integer                             :: I02
            integer                             :: i
            double precision                    :: dlgama
            double precision                    :: w1
            double precision                    :: w2

            I01 = 0
            I12 = 0
            I02 = 0
            do i = 1,en
                  if (s(jstar).le.y(i).and.y(i).lt.s(jstar+1)) then
                        I01 = I01+1
                  else if (s(jstar+1).le.y(i).and.y(i).lt.s(jstar+2)) then
                        I12 = I12+1
                  end if
            end do
            I02 = I01+I12

            w1 = (s(jstar+1)-s(jstar))/(s(jstar+2)-s(jstar))
            w2 = 1.d0-W1

            logIW = log(h(jstar))*I01 +log(h(jstar+1))*I12 &
                              -(w1*log(h(jstar))+w2*log(h(jstar+1)))*I02 &
                        -h(jstar)*(s(jstar+1)-s(jstar)) -h(jstar+1)*(s(jstar+2)-s(jstar+1)) &
                              +h(jstar)**w1*h(jstar+1)**w2*(s(jstar+2)-s(jstar)) &
                  ! model prior
                        +log(lambda) -log(k+1.d0) &
                  ! positions prior
                        +log(2.d0*(k+1.d0)*(2.d0*k+3.d0)) -2.d0*log(Leng) &
                        +log(s(jstar+1)-s(jstar)) &
                        +log(s(jstar+2)-s(jstar+1)) &
                        -log(s(jstar+2)-s(jstar)) &
                  ! heights prior
                        +alpha*log(beta) -dlgama(alpha) &
                        +(alpha-1)*( w2*log(h(jstar)) +w1*log(h(jstar+1)) ) &
                        -beta*( h(jstar) +h(jstar+1) -h(jstar)**w1*h(jstar+1)**w2 ) &
                  ! auxiliary values (position+kill+[height])
                        +log(Leng) -log(k+1.d0) &
                  ! jacobian
                        +2.d0*log(h(jstar)+h(jstar+1)) &
                        -w1*log(h(jstar)) &
                        -w2*log(h(jstar+1))

      end subroutine AnnImpWeight

      ! ================================================================
      ! TEMPERATURE FUNCTION
      ! ================================================================

      subroutine gimel(gt,t,Tau)

            implicit none

            integer, intent(in)           :: t
            integer, intent(in)           :: Tau
            double precision, intent(out) :: gt

            if (t.lt.0.or.t.gt.Tau) stop 'gimel>>t.lt.0.or.t.gt.Tau'

            gt = 1.d0*t/Tau

      end subroutine gimel

      ! ================================================================
      ! ANNEALING IMPORTANCE SWEEP
      ! ================================================================

      subroutine AisRjSplitSweep(k,lambda,h,alpha,cc,dd,aspr,beta,ee,ff,s,jstar,y,en,gt)

            implicit none

            integer, intent(in)                 :: k
            integer, intent(inout)              :: jstar
            double precision, intent(inout)     :: s(0:k+2)
            double precision, intent(inout)     :: h(0:k+1)
            double precision, intent(inout)     :: alpha
            double precision, intent(inout)     :: beta
            double precision, intent(in)        :: lambda
            double precision, intent(in)        :: gt
            double precision, intent(in)        :: aspr
            double precision, intent(in)        :: cc
            double precision, intent(in)        :: dd
            double precision, intent(in)        :: ee
            double precision, intent(in)        :: ff
            integer, intent(in)                 :: en
            integer, intent(in)                 :: y(en)

            integer           :: d
            integer           :: i_bl
            integer           :: n_bl
            integer           :: ran_bl(1:2*k+6)
            integer           :: bl
            integer           :: jj
            integer           :: i
            integer           :: jup
            integer           :: I10
            integer           :: I01
            integer           :: I12
            integer           :: I02
            integer           :: I23
            integer           :: I10New
            integer           :: I01New
            integer           :: I12New
            integer           :: I02New
            integer           :: I23New
            double precision  :: u
            double precision  :: logAccRat
            double precision  :: w1
            double precision  :: w2
            double precision  :: dlgama
            double precision  :: w1New
            double precision  :: w2New
            double precision  :: alphaNew
            double precision  :: hNew
            double precision  :: sNew
            double precision  :: par1
            double precision  :: par2
            double precision  :: PrJ(0:k)

            ! SET THE PARAMETERS
            n_bl = 2*k+6

            ! GENERATE A RANDOM BLOCK
            do jj = 1,n_bl
                  ran_bl(jj) = jj
            end do
            call genprm(ran_bl(1:n_bl),n_bl)

            ! DO THE RANDOM PERMUTATION SWEEP 
            do i_bl = 1, n_bl

            ! CHOOSE A BLOCK TO UPDATE
            bl = ran_bl(i_bl)

            ! UPDATE ALPHA
            ! ====================================================

            if (bl.eq.1) then

            call rnguniform(u)
            alphaNew = alpha*aspr**(u-0.5d0)

            w1 = (s(jstar+1)-s(jstar))/(s(jstar+2)-s(jstar))
            w2 = 1.d0-w1

            logAccRat = (alphaNew-alpha)*( &
                              sum(log(h(0:jstar-1))) &
                              +((1-gt)*w1+gt)*log(h(jstar)) &
                              +((1-gt)*w2+gt)*log(h(jstar+1)) &
                              +sum(log(h(jstar+2:k+1))) &
                              +log(beta)*(k+1+gt) -dd &
                                          ) &
                        -(dlgama(alphaNew)-dlgama(alpha))*(k+1+gt) &
                        +cc*(log(alphaNew)-log(alpha))

            call rnguniform(u)
            if (logAccRat.ge.log(u)) alpha = alphaNew

            ! UPDATE BETA
            ! ====================================================

            else if (bl.eq.2) then

            w1 = (s(jstar+1)-s(jstar))/(s(jstar+2)-s(jstar))
            w2 = 1.d0-w1

            par1 = alpha*(k+1+gt)+ee
            par2 = sum(h(0:jstar-1)) &
                  +(1-gt)*h(jstar)**w1*h(jstar+1)**w2 +gt*h(jstar) +gt*h(jstar+1) &
                  +sum(h(jstar+2:k+1)) &
                  +ff

            call rnggamma(beta,par1)
            beta = beta/par2


            ! UPDATE HEIGHTS
            ! ====================================================

            else if (3.le.bl.and.bl.le.k+4) then

            jup = bl-3 ! UPDATE THE jup-th H

            ! PROPOSE THE CANDIDATE h
            call rnguniform(u)
            u = -0.5d0+u
            hNew = exp(u)*h(jup)

            ! COMPUTE THE ACCEPTANCE RATIO

            if (jup.le.jstar-1.or.jstar+2.le.jup) then

                  I01 = 0
                  do i = 1,en
                        if (y(i).ge.s(jup).and.y(i).lt.s(jup+1)) then
                              I01 = I01+1
                        end if
                  end do

                  logAccRat = (log(hNew)-log(h(jup)))*I01 &
                              -(hNew-h(jup))*(s(jup+1)-s(jup)) &
                              -beta*(hNew-h(jup))      &
                              +alpha*(log(hNew)-log(h(jup)))

            else if (jup.eq.jstar) then

                  I01 = 0
                  I12 = 0
                  I02 = 0
                  do i = 1,en
                        if (y(i).ge.s(jstar).and.y(i).lt.s(jstar+1)) then
                              I01 = I01+1
                        else if (y(i).ge.s(jstar+1).and.y(i).lt.s(jstar+2)) then
                              I12 = I12+1
                        end if
                  end do
                  I02 = I01+I12

                  w1 = (s(jstar+1)-s(jstar))/(s(jstar+2)-s(jstar))
                  W2 = 1.d0-w1

                  logAccRat = (log(hNew)-log(h(jstar)))*(&
                                    (1-gt)*w1*I02 +gt*I01 +alpha*(1-gt)*w1+(alpha-1)*gt+1 &
                                          ) &
                              -(hNew**w1-h(jstar)**w1)*h(jstar+1)**w2*(1-gt)*(s(jstar+2)-s(jstar)+beta) &
                              -(hNew-h(jstar))*gt*(s(jstar+1)-s(jstar)+beta) &
                              -( log(hNew+h(jstar+1))-log(h(jstar)+h(jstar+1)) ) *2 *(1-gt)

            else if (jup.eq.jstar+1) then

                  I01 = 0
                  I12 = 0
                  I02 = 0
                  do i = 1,en
                        if (y(i).ge.s(jstar).and.y(i).lt.s(jstar+1)) then
                              I01 = I01+1
                        else if (y(i).ge.s(jstar+1).and.y(i).lt.s(jstar+2)) then
                              I12 = I12+1
                        end if
                  end do
                  I02 = I01+I12

                  w1 = (s(jstar+1)-s(jstar))/(s(jstar+2)-s(jstar))
                  W2 = 1.d0-w1

                  logAccRat = (log(hNew)-log(h(jstar+1)))*(&
                                          (1-gt)*w2*I02 +gt*I12 &
                                          +alpha*(1-gt)*w2 +(alpha-1)*gt+1 &
                                          ) &
                              -(hNew**w2-h(jstar+1)**w2)*h(jstar)**w1*(1-gt)*(s(jstar+2)-s(jstar)+beta) &
                              -(hNew-h(jstar+1))*gt*(s(jstar+2)-s(jstar+1)+beta) &
                              -( log(hNew+h(jstar))-log(h(jstar+1)+h(jstar)) ) *2 *(1-gt)

            end if

            ! ACCEPT/REJECT
            call rnguniform(u)
            if (logAccRat.gt.log(u)) then
                  h(jup) = hNew
            end if

            ! UPDATE POSITIONS
            ! ====================================================

            else if (k+5.le.bl.and.bl.le.2*k+5) then

                  jup = bl-k-4

                  ! PROPOSE THE CANDIDATE s
                  call rnguniform(u)
                  sNew = s(jup-1) +(s(jup+1)-s(jup-1))*u

                  ! COMPUTE THE LOG ACCEPTANCE RATIO

                  if (jup.eq.jstar) then

                        w1 = (s(jstar+1)-s(jstar)) / (s(jstar+2)-s(jstar))
                        w2 = 1.d0-w1
                        w1New = (s(jstar+1)-sNew) / (s(jstar+2)-sNew)
                        w2New = 1.d0-w1New

                        I10 = 0
                        I01 = 0
                        I12 = 0
                        I02 = 0
                        I10New = 0
                        I01New = 0
                        I12New = 0
                        I02New = 0
                        do i = 1,en
                              if (s(jstar-1).le.y(i).and.y(i).lt.s(jstar)) then
                                    I10 = I10+1
                              else if (s(jstar).le.y(i).and.y(i).lt.s(jstar+1)) then
                                    I01 = I01+1
                              else if (s(jstar+1).le.y(i).and.y(i).lt.s(jstar+2)) then
                                    I12 = I12+1
                              end if
                              if (s(jstar-1).le.y(i).and.y(i).lt.sNew) then
                                    I10New = I10New+1
                              else if (sNew.le.y(i).and.y(i).lt.s(jstar+1)) then
                                    I01New = I01New+1
                              else if (s(jstar+1).le.y(i).and.y(i).lt.s(jstar+2)) then
                                    I12New = I12New+1
                              end if
                        end do
                        I02 = I01+I12
                        I02New = I01New+I12New

                        par2 = (alpha-1)*(1-gt)*w1New*log(h(jstar)) &
                                    +(alpha-1)*(1-gt)*w2New*log(h(jstar+1)) &
                              -beta*(1-gt)*h(jstar)**w1New*h(jstar+1)**w2New &
                              +log(sNew-s(jstar-1)) &
                                    +(1-gt)*log(s(jstar+2)-sNew) &
                                    +gt*log(s(jstar+1)-sNew) &
                              +log(h(jstar-1))*I10New &
                                    +log(h(jstar))*((1-gt)*w1New*I02New +gt*I01New) &
                                    +log(h(jstar+1))*((1-gt)*w2New*I02New) &
                              -h(jstar-1)*sNew &
                                    -(1-gt) *h(jstar)**w1New *h(jstar+1)**w2New *(s(jstar+2)-sNew) &
                                    -gt *h(jstar)*(-sNew) &
                              +(1-gt)*w1New*log(h(jstar)) +(1-gt)*w2New*log(h(jstar+1))

                        par1 = (alpha-1)*(1-gt)*w1*log(h(jstar)) &
                                    +(alpha-1)*(1-gt)*w2*log(h(jstar+1)) &
                              -beta*(1-gt)*h(jstar)**w1*h(jstar+1)**w2 &
                              +log(s(jstar)-s(jstar-1)) &
                                    +(1-gt)*log(s(jstar+2)-s(jstar)) &
                                    +gt*log(s(jstar+1)-s(jstar)) &
                              +log(h(jstar-1))*I10 &
                                    +log(h(jstar))*((1-gt)*w1*I02 +gt*I01) &
                                    +log(h(jstar+1))*((1-gt)*w2*I02) &
                              -h(jstar-1)*s(jstar) &
                                    -(1-gt)*h(jstar)**w1*h(jstar+1)**w2*(s(jstar+2)-s(jstar)) &
                                    -gt*h(jstar)*(-s(jstar)) &
                              +(1-gt)*w1*log(h(jstar)) +(1-gt)*w2*log(h(jstar+1))

                        logAccRat = par2-par1

                  else if (jup.eq.jstar+1) then

                        w1 = (s(jstar+1)-s(jstar))/(s(jstar+2)-s(jstar))
                        w2 = 1.d0-w1
                        w1New = (sNew-s(jstar))/(s(jstar+2)-s(jstar))
                        w2New = 1.d0-w1New

                        I01 = 0
                        I12 = 0
                        I02 = 0
                        I01New = 0
                        I12New = 0
                        I02New = 0
                        do i = 1,en
                              if (s(jstar).le.y(i).and.y(i).lt.s(jstar+1)) then
                                    I01 = I01+1
                              else if (s(jstar+1).le.y(i).and.y(i).lt.s(jstar+2)) then
                                    I12 = I12+1
                              end if
                              if (s(jstar).le.y(i).and.y(i).lt.sNew) then
                                    I01New = I01New+1
                              else if (sNew.le.y(i).and.y(i).lt.s(jstar+2)) then
                                    I12New = I12New+1
                              end if
                        end do
                        I02 = I01+I12
                        I02New = I01New+I12New

                        par2 = (alpha-1)*(1-gt)*w1New*log(h(jstar)) &
                                    +(alpha-1)*(1-gt)*w2New*log(h(jstar+1)) &
                              -beta*(1-gt)*h(jstar)**w1New*h(jstar+1)**w2New &
                              +gt*log(sNew-s(jstar)) +gt*log(s(jstar+2)-sNew) &
                              +log(h(jstar))*((1-gt)*w1New*I02New+gt*I01New) &
                                    +log(h(jstar+1))*((1-gt)*w2New*I02New+gt*I12New) &
                              -(1-gt)*h(jstar)**w1New*h(jstar+1)**w2New*(s(jstar+2)-s(jstar)) &
                                    -gt*h(jstar)*(sNew-s(jstar)) &
                                    -gt*h(jstar+1)*(s(jstar+2)-sNew) &
                              +(1-gt)*w1New*log(h(jstar)) +(1-gt)*w2New*log(h(jstar+1))

                        par1 = (alpha-1)*(1-gt)*w1*log(h(jstar)) &
                                    +(alpha-1)*(1-gt)*w2*log(h(jstar+1)) &
                              -beta*(1-gt)*h(jstar)**w1*h(jstar+1)**w2 &
                              +gt*log(s(jstar+1)-s(jstar)) +gt*log(s(jstar+2)-s(jstar+1)) &
                              +log(h(jstar))*((1-gt)*w1*I02+gt*I01) &
                                    +log(h(jstar+1))*((1-gt)*w2*I02+gt*I12) &
                              -(1-gt)*h(jstar)**w1*h(jstar+1)**w2*(s(jstar+2)-s(jstar)) &
                                    -gt*h(jstar)*(s(jstar+1)-s(jstar)) &
                                    -gt*h(jstar+1)*(s(jstar+2)-s(jstar+1)) &
                              +(1-gt)*w1*log(h(jstar)) +(1-gt)*w2*log(h(jstar+1))

                        logAccRat = par2-par1

                  else if (jup.eq.jstar+2) then

                        w1 = (s(jstar+1)-s(jstar)) / (s(jstar+2)-s(jstar))
                        w2 = 1.d0-w1
                        w1New = (s(jstar+1)-s(jstar)) / (sNew-s(jstar))
                        w2New = 1.d0-w1New

                        I01 = 0
                        I12 = 0
                        I02 = 0
                        I23 = 0
                        I01New = 0
                        I12New = 0
                        I02New= 0
                        I23New = 0
                        do i = 1,en
                              if (s(jstar).le.y(i).and.y(i).lt.s(jstar+1)) then
                                    I01 = I01+1
                              else if (s(jstar+1).le.y(i).and.y(i).lt.s(jstar+2)) then
                                    I12 = I12+1
                              else if (s(jstar+2).le.y(i).and.y(i).lt.s(jstar+3)) then
                                    I23 = I23+1
                              end if
                              if (s(jstar).le.y(i).and.y(i).lt.s(jstar+1)) then
                                    I01New = I01New+1
                              else if (s(jstar+1).le.y(i).and.y(i).lt.sNew) then
                                    I12New = I12New+1
                              else if (sNew.le.y(i).and.y(i).lt.s(jstar+3)) then
                                    I23New = I23New+1
                              end if
                        end do
                        I02 = I01+I12
                        I02New = I01New+I12New

                        par2 = (alpha-1)*(1-gt)*w1New*log(h(jstar)) &
                                    +(alpha-1)*(1-gt)*w2New*log(h(jstar+1)) &
                              -beta*(1-gt)*h(jstar)**w1New*h(jstar+1)**w2New &
                              +(1-gt)*log(sNew-s(jstar)) &
                                    +gt*log(sNew-s(jstar+1)) &
                                    +log(s(jstar+3)-sNew) &
                              +log(h(jstar))*((1-gt)*w1New*I02New) &
                                    +log(h(jstar+1))*((1-gt)*w2New*I02New+gt*I12New) &
                                    +log(h(jstar+2))*I23New &
                              -(1-gt)*h(jstar)**w1New*h(jstar+1)**w2New*(sNew-s(jstar)) &
                                    -gt*h(jstar+1)*(sNew-s(jstar+1)) &
                                    -h(jstar+2)*(s(jstar+3)-sNew) &
                              +(1-gt)*w1New*log(h(jstar)) +(1-gt)*w2New*log(h(jstar+1))

                        par1 = (alpha-1)*(1-gt)*w1*log(h(jstar)) &
                                    +(alpha-1)*(1-gt)*w2*log(h(jstar+1)) &
                              -beta*(1-gt)*h(jstar)**w1*h(jstar+1)**w2 &
                              +(1-gt)*log(s(jstar+2)-s(jstar)) &
                                    +gt*log(s(jstar+2)-s(jstar+1)) &
                                    +log(s(jstar+3)-s(jstar+2)) &
                              +log(h(jstar))*((1-gt)*w1*I02) &
                                    +log(h(jstar+1))*((1-gt)*w2*I02+gt*I12) &
                                    +log(h(jstar+2))*I23 &
                              -(1-gt)*h(jstar)**w1*h(jstar+1)**w2*(s(jstar+2)-s(jstar)) &
                                    -gt*h(jstar+1)*(s(jstar+2)-s(jstar+1)) &
                                    -h(jstar+2)*(s(jstar+3)-s(jstar+2)) &
                              +(1-gt)*w1*log(h(jstar)) +(1-gt)*w2*log(h(jstar+1))

                        logAccRat = par2-par1

                  else

                        I10 = 0
                        I01 = 0
                        I10New = 0
                        I01New = 0
                        do i = 1,en
                              if (s(jup-1).le.y(i).and.y(i).lt.s(jup)) then
                                    I10 = I10+1
                              else if (s(jup).le.y(i).and.y(i).lt.s(jup+1)) then
                                    I01 = I01+1
                              end if
                              if (s(jup-1).le.y(i).and.y(i).lt.sNew) then
                                    I10New = I10New+1
                              else if (sNew.le.y(i).and.y(i).lt.s(jup+1)) then
                                    I01New = I01New+1
                              end if
                        end do

                        logAccRat = log(h(jup-1))*(I10New-I10) &
                                    +log(h(jup))*(I01New-I01) &
                                    +(h(jup)-h(jup-1))*(sNew-s(jup)) &
                                    +log(sNew-s(jup-1))-log(s(jup)-s(jup-1)) &
                                    +log(s(jup+1)-sNew)-log(s(jup+1)-s(jup))

                  end if

                  ! ACCEPT/REJECT
                  call rnguniform(u)
                  if (logAccRat.gt.log(u)) s(jup)=sNew

            else if (bl.eq.2*k+6) then

            ! COMPUTE THE PROBABILITY UP TO A NORMCONST
                  do jup = 0,k
                        call logRhoDen(PrJ(jup), &
                                          k,lambda, &
                                          h(0:k+1),alpha,cc,dd,beta,ee,ff, &
                                          s(0:k+2),jup, &
                                          y(1:en),en,gt)
                  end do
                  ! FIND THE MAXIMUM AND CORRECT (NATURAL CORRECTION)
                  par1 = PrJ(0)
                  do jup = 1,k
                        if (par1.lt.PrJ(jup)) par1 = PrJ(jup)
                  end do
                  do jup = 0,k
                        PrJ(jup) = PrJ(jup)-par1
                  end do
                  ! COMUTE THE CUMULATIVE UP TO A NORMCONST
                  PrJ(0) = exp(PrJ(0))
                  do jup = 1,k
                        PrJ(jup) = PrJ(jup-1)+exp(PrJ(jup))
                  end do
                  ! SET THE NORMCONST OF THE CUMULATIVE
                  par1 = PrJ(k)

                  ! ARTIFICIAL CORRECTION
                  if (par1.eq.0.d0) then
                        do jup = 0,k
                              PrJ(jup) = jup+1.d0
                              par1 = k+1.d0
                        end do
                  end if

            ! SAMPLE THE JSTAR
                  call rnguniform(u)
                  jup = 0
                  do 
                        if (PrJ(jup).ge.u*par1) exit
                        if (jup.ge.k+1) then
                              stop 'AisRjSplitSweep::(jstar.ge.k+1)'
                        end if

                        jup = jup+1
                  end do

                  jstar = jup

            end if

            end do

      end subroutine AisRjSplitSweep

      ! ================================================================
      ! RHO DISTRUBUTION UP TO A NORMALIZING CONSTANT
      ! ================================================================

      subroutine logRhoDen(logRho,k,lambda,h,alpha,cc,dd,beta,ee,ff,s,jstar,y,en,gt)

            double precision, intent(out) :: logRho
            integer, intent(in)           :: k
            double precision, intent(in)  :: lambda
            double precision, intent(in)  :: h(0:k+1)
            double precision, intent(in)  :: alpha
            double precision, intent(in)  :: cc
            double precision, intent(in)  :: dd
            double precision, intent(in)  :: beta
            double precision, intent(in)  :: ee
            double precision, intent(in)  :: ff
            double precision, intent(in)  :: s(0:k+2)
            integer, intent(in)           :: jstar
            double precision, intent(in)  :: gt
            integer, intent(in)           :: en
            integer, intent(in)           :: y(en)

            integer           :: i
            integer           :: j
            integer           :: countY
            integer           :: countY1
            integer           :: countY2
            double precision  :: w1,w2
            double precision  :: dlgama

            logRho = 0.d0

            do j = 0,jstar-1
            ! counter
                  countY = 0
                  do i = 1,en
                        if (y(i).ge.s(j).and.y(i).lt.s(j+1)) then
                              countY = countY+1
                        end if
                  end do

                  logRho = logRho &
            ! likelihood
                        +log(h(j))*countY &
                        -h(j)*(s(j+1)-s(j)) &
            ! prior
                        +log(s(j+1)-s(j)) &
                        -beta*h(j) &
                        +(alpha-1)*log(h(j))
            end do

            ! counter
            countY1 = 0
            countY2 = 0
            do i = 1,en
                  if (y(i).ge.s(jstar).and.y(i).lt.s(jstar+1)) then
                        countY1 = countY1+1
                  else if (y(i).ge.s(jstar+1).and.y(i).lt.s(jstar+2)) then
                        countY2 = countY2+1
                  end if
            end do

            w1 = (s(jstar+1)-s(jstar))/(s(jstar+2)-s(jstar))
            w2 = 1-w1

            logRho = logRho &
                        +log(h(jstar)) &
                              *( &
                                    (1-gt)*w1*(countY1+countY2) &
                                          +gt*countY1 &
                              ) &
                        +log(h(jstar+1)) &
                              *( &
                                    (1-gt)*w2*(countY1+countY2) &
                                          +gt*countY2 &
                              ) &
                        -(1-gt)*h(jstar)**w1 &
                              *h(jstar+1)**w2 &
                              *(s(jstar+2)-s(jstar)) &
                        -gt*h(jstar)*(s(jstar+1)-s(jstar)) &
                        -gt*h(jstar+1)*(s(jstar+2)-s(jstar+1)) &
                        +(1-gt)*log(s(jstar+2)-s(jstar)) &
                        +gt*log(s(jstar+1)-s(jstar)) &
                        +gt*log(s(jstar+2)-s(jstar+1)) &
                        +(k+1+gt)*alpha*log(beta) &
                        -(k+1+gt)*dlgama(alpha) &
                        -beta*(1-gt)*h(jstar)**w1*h(jstar+1)**w2 &
                        -beta*gt*h(jstar) &
                        -beta*gt*h(jstar+1) &
                        +(alpha-1)*((1-gt)*w1+gt)*log(h(jstar)) &
                        +(alpha-1)*((1-gt)*w2+gt)*log(h(jstar+1)) 

            do j = jstar+2,k+1
            ! counter
                  countY = 0
                  do i = 1,en
                        if (y(i).ge.s(j).and.y(i).lt.s(j+1)) then
                              countY = countY+1
                        end if
                  end do

                  logRho = logRho &
            ! likelihood
                        +log(h(j))*countY &
                        -h(j)*(s(j+1)-s(j)) &
            ! prior
                        +log(s(j+1)-s(j)) &
                        -beta*h(j) &
                        +(alpha-1)*log(h(j))
            end do
            logRho = logRho &
                        -dd*alpha+(cc-1)*log(alpha) &
                        -ff*beta+(ee-1)*log(beta)

            logRho = logRho &
                        +(1-gt)*w1*log(h(jstar)) &
                        +(1-gt)*w2*log(h(jstar+1)) &
                        -2*(1-gt)*log(h(jstar)+h(jstar+1))

      end subroutine logRhoDen


end module aisrjupdates_mod


