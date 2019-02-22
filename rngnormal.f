c	The Author of this code is Prof Peter Green (University of Bristol), 
c	and I borrowed this bit of code from the Nmix program.

c	FROM PJG GAUSS4.F

	subroutine rngnormal(z,n)
c..	generates an n-vector of i.i.d. N(0,1) r.v.'s z
c	by Box-Mueller method
	implicit none
	integer n,n1,i
	double precision z(n)
	double precision u,v
	n1 = n-1
	do 1 i = 1,n1,2
	call rnguniform(u)
	u = sqrt(-2.0d0*log(u))
	call rnguniform(v)
	v = 6.28318530717959d0*v
	z(i) = u*sin(v)
1	z(i+1) = u*cos(v)
	if(mod(n,2).eq.0) return
	call rnguniform(u)
	call rnguniform(v)
	z(n) = sqrt(-2.0d0*log(u))*sin(6.28318530717959d0*v)
	return
	end

	function rgauss()
	implicit none
	double precision rgauss
	integer n
	double precision z(1)
	n = 1
	call rngnormal(z,n)
	rgauss = z(1)
	return
	end
