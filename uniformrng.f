c Copyrigtht 2012 Georgios Karagiannis
c
c This file is part of AISRJ_CPT.
c
c AISRJ_CPT is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation version 2 of the License.
c
c AISRJ_CPT is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with AISRJ_CPT.  If not, see <http://www.gnu.org/licenses/>.

c ----------------------------------------------------------------------

c REFERENCES:
c -----------
c
c Karagiannis, G., & Andrieu, C. (2013).
c Annealed importance sampling reversible jump MCMC algorithms.
c Journal of Computational and Graphical Statistics, 22(3), 623-648.
c
c Georgios Karagiannis
c School of Mathematics, University of Bristol
c University Walk, Bristol, BS8 1TW, UK
c Email : Georgios.Karagiannis@pnnl.gov
c Email (current): georgios-stats@gmail.com
c
c Christophe Andrieu
c School of Mathematics, University of Bristol
c University Walk, Bristol, BS8 1TW, UK
c Email: C.Andrieu@bristol.ac.uk
c
      subroutine init_seedrng()
      call random_seed()
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on (0,1)-real-interval
c-----------------------------------------------------------------------
      subroutine rnguniform(rnd)
      double precision rnd
      call random_number(rnd)
      return
      end


