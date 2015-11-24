
# Contact details :

# Georgios Karagiannis © 2012
# School of Mathematics, University of Bristol
# University Walk, Bristol, BS8 1TW, UK
# Email (current): Georgios.Karagiannis@pnnl.gov

# Christophe Andrieu
# School of Mathematics, University of Bristol
# University Walk, Bristol, BS8 1TW, UK
# Email: C.Andrieu@bristol.ac.uk

FC = gfortran  #  choose a fortran compiler
FFLAGS = -O3   #  choose a flag

# Generate the executable

aisrj.exe: algama.o genprm.o mt19937ar.o rngnormal.o rnggamma.o fixedupdates_mod.o aisrjupdates_mod.o aisrjcpt_pro.o
	$(FC) $(FFLAGS) -o aisrj.exe algama.o mt19937ar.o genprm.o rngnormal.o rnggamma.o fixedupdates_mod.o aisrjupdates_mod.o aisrjcpt_pro.o

algama.o: algama.f
	$(FC) $(FFLAGS) -c algama.f

mt19937ar.o: mt19937ar.f
	$(FC) $(FFLAGS) -c mt19937ar.f
   
genprm.o: mt19937ar.o genprm.f
	$(FC) $(FFLAGS) -c genprm.f
   
rngnormal.o: mt19937ar.o rngnormal.f
	$(FC) $(FFLAGS) -c rngnormal.f
   
rnggamma.o: mt19937ar.o rnggamma.f90
	$(FC) $(FFLAGS) -c rnggamma.f90
   
fixedupdates_mod.mod: fixedupdates_mod.o fixedupdates_mod.f90 rnggamma.o rngnormal.o genprm.o mt19937ar.o algama.o
	$(FC) $(FFLAGS) -c fixedupdates_mod.f90 
fixedupdates_mod.o: fixedupdates_mod.f90 rnggamma.o rngnormal.o genprm.o mt19937ar.o algama.o
	$(FC) $(FFLAGS) -c fixedupdates_mod.f90  
   
aisrjupdates_mod.mod: aisrjupdates_mod.o aisrjupdates_mod.f90 rnggamma.o rngnormal.o genprm.o mt19937ar.o algama.o
	$(FC) $(FFLAGS) -c aisrjupdates_mod.f90 
aisrjupdates_mod.o: aisrjupdates_mod.f90 rnggamma.o rngnormal.o genprm.o mt19937ar.o algama.o
	$(FC) $(FFLAGS) -c aisrjupdates_mod.f90
   
aisrjcpt_pro.o: fixedupdates_mod.mod aisrjupdates_mod.mod rnggamma.o rngnormal.o genprm.o mt19937ar.o algama.o aisrjcpt_pro.f90
	$(FC) $(FFLAGS) -c aisrjcpt_pro.f90

# Clean auxiliary files

clean:
	rm aisrjupdates_mod.o fixedupdates_mod.o rnggamma.o rngnormal.o genprm.o mt19937ar.o algama.o aisrjcpt_pro.o fixedupdates_mod.mod aisrjupdates_mod.mod


