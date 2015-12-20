
# --------------------------------------------------------------------------------
# 
# Copyrigtht 2012 Georgios Karagiannis
# 
# This file is part of AISRJ_CPT.
# 
# AISRJ_CPT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# AISRJ_CPT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with AISRJ_CPT.  If not, see <http://www.gnu.org/licenses/>.
# 
# --------------------------------------------------------------------------------

# References:
# 
# Karagiannis, G., & Andrieu, C. (2013). 
# Annealed importance sampling reversible jump MCMC algorithms. 
# Journal of Computational and Graphical Statistics, 22(3), 623-648.
# 
# Georgios Karagiannis
# School of Mathematics, University of Bristol
# University Walk, Bristol, BS8 1TW, UK
# Email : Georgios.Karagiannis@pnnl.gov
# Email (current): georgios-stats@gmail.com
# 
# Christophe Andrieu
# School of Mathematics, University of Bristol
# University Walk, Bristol, BS8 1TW, UK
# Email: C.Andrieu@bristol.ac.uk


FC = gfortran  #  choose a fortran compiler
FFLAGS = -O3   #  choose a flag

# Generate the executable

aisrj.exe: algama.o genprm.o uniformrng.o rngnormal.o rnggamma.o fixedupdates_mod.o aisrjupdates_mod.o aisrjcpt_pro.o
	$(FC) $(FFLAGS) -o aisrj.exe algama.o uniformrng.o genprm.o rngnormal.o rnggamma.o fixedupdates_mod.o aisrjupdates_mod.o aisrjcpt_pro.o

algama.o: algama.f
	$(FC) $(FFLAGS) -c algama.f

uniformrng.o: uniformrng.f
	$(FC) $(FFLAGS) -c uniformrng.f
   
genprm.o: uniformrng.o genprm.f
	$(FC) $(FFLAGS) -c genprm.f
   
rngnormal.o: uniformrng.o rngnormal.f
	$(FC) $(FFLAGS) -c rngnormal.f
   
rnggamma.o: uniformrng.o rnggamma.f90
	$(FC) $(FFLAGS) -c rnggamma.f90
   
fixedupdates_mod.mod: fixedupdates_mod.o fixedupdates_mod.f90 rnggamma.o rngnormal.o genprm.o uniformrng.o algama.o
	$(FC) $(FFLAGS) -c fixedupdates_mod.f90 
fixedupdates_mod.o: fixedupdates_mod.f90 rnggamma.o rngnormal.o genprm.o uniformrng.o algama.o
	$(FC) $(FFLAGS) -c fixedupdates_mod.f90  
   
aisrjupdates_mod.mod: aisrjupdates_mod.o aisrjupdates_mod.f90 rnggamma.o rngnormal.o genprm.o uniformrng.o algama.o
	$(FC) $(FFLAGS) -c aisrjupdates_mod.f90 
aisrjupdates_mod.o: aisrjupdates_mod.f90 rnggamma.o rngnormal.o genprm.o uniformrng.o algama.o
	$(FC) $(FFLAGS) -c aisrjupdates_mod.f90
   
aisrjcpt_pro.o: fixedupdates_mod.mod aisrjupdates_mod.mod rnggamma.o rngnormal.o genprm.o uniformrng.o algama.o aisrjcpt_pro.f90
	$(FC) $(FFLAGS) -c aisrjcpt_pro.f90

# Clean auxiliary files

clean:
	rm aisrjupdates_mod.o fixedupdates_mod.o rnggamma.o rngnormal.o genprm.o uniformrng.o algama.o aisrjcpt_pro.o fixedupdates_mod.mod aisrjupdates_mod.mod


