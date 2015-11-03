! CONTACT DETAILS :

! Georgios Karagiannis (Copyright 2012)
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email (current): Georgios.Karagiannis@pnnl.gov

! Christophe Andrieu
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email: C.Andrieu@bristol.ac.uk

! ----------------------------------------------------------------------

subroutine rnggamma(x,alpha)

implicit none

double precision, intent(out) :: x
double precision, intent(in)  :: alpha

double precision      :: u,p,v,y,yy(1)
double precision      :: e,b,d,c

if (alpha.eq.1.d0) then

      call rnguniform(u)
      x = -log(u)
      return

else if (alpha.lt.1.d0) then

      e = exp(1.d0)
      b = (e+alpha)/e

      do
            call rnguniform(u)
            p = b*u
            if (p.gt.1.d0) then
                  x = -log((b-p)/alpha)
                  call rnguniform(u)
                  if (u.le.x**(alpha-1)) return
            else
                  x = p**(1/alpha)
                  call rnguniform(u)
                  if (u.le.exp(-x)) return
            end if
      end do

else if (alpha.gt.1.d0) then

      d = alpha-1.d0/3
      c = 1.d0/sqrt(9.d0*d)
      do
            do
                  call rngnormal(yy,1)
                  y = yy(1)
                  v = (1+c*y)**3
                  if (v.gt.0.d0) exit
            end do
            call rnguniform(u)
            if (u.lt.1-0.0331*y**4) then
                  x = d*v
                  return
            end if
            if (log(u).lt.0.5d0*y**2+d*(1-v+log(v))) then
                  x = d*v
                  return
            end if
      end do
end if

end subroutine rnggamma


