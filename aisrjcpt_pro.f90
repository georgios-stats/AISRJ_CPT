! CONTACT DETAILS :

! Georgios Karagiannis
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email (current): Georgios.Karagiannis@pnnl.gov

! Christophe Andrieu
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email: C.Andrieu@bristol.ac.uk

! ----------------------------------------------------------------------

program aisrjcpt_pro

      use fixedupdates_mod, only      : updateH,updateS,updateBeta,updateAlpha
      use aisrjupdates_mod, only      : AnnImpWeight
      use aisrjupdates_mod, only      : AisRjSplitSweep,gimel

      implicit none

!      PARAMETERS OF THE ALGORITHM

      integer                 :: nsweep
      integer                 :: nburnin
      integer                 :: Tau

!      DATASET

      integer                 :: enmax
      parameter               ( enmax = 200 )
      integer                 :: en
      integer                 :: y(enmax)

!      FIXED PARAMETERS
      integer                 :: kmaxmax
      parameter               ( kmaxmax = 40 )
      integer                 :: kmax
      double precision        :: Leng
      double precision        :: lambda
      double precision        :: cc
      double precision        :: dd
      double precision        :: ee
      double precision        :: ff
      double precision        :: aspr

!      RANDOM PARAMETERS

      integer                 :: k
      integer                 :: k_ais
      integer                 :: k1
      integer                 :: k2
      double precision        :: h(0:kmaxmax)
      double precision        :: s(0:kmaxmax+1)
      double precision        :: alpha
      double precision        :: beta
      double precision        :: h_ais(0:kmaxmax)
      double precision        :: s_ais(0:kmaxmax+1)
      double precision        :: alpha_ais
      double precision        :: beta_ais
      double precision        :: h1,h2,s0

      double precision        :: ExpAccPrRJ
      double precision        :: ExpAccPrSplit
      double precision        :: ExpAccPrMerge
      integer                 :: nSplit
      integer                 :: nMerge

!      ALGORITHM VARIABLES

      double precision        :: pb(0:kmaxmax)
      double precision        :: pd(0:kmaxmax)
      integer                 :: Qtype

!      OTHER

      character(len=1)        :: move
      integer                 :: i
      integer                 :: j
      integer                 :: iter
      integer                 :: jup
      integer                 :: jstar
      integer                 :: seed
      double precision        :: u
      double precision        :: logIW
      double precision        :: AccPr
      integer                 :: t
      double precision        :: logAIW
      double precision        :: gt
      double precision        :: Dgt
      double precision        :: pr_k(0:kmaxmax)
      double precision        :: cost

!     ESTIMATES

      integer                 :: countK

!     ARGUMENTS

      integer                 :: ia
      integer                 :: iargc
      character(len=20)       :: word
      character(len=14)       :: specT
      integer                 :: istdT
      integer                 :: kseed

!      SET THE PARAMETERS OF THE ALGORITHM
!      ===================================

      nsweep = 2*10**5
      nburnin = max(10,nsweep/1000)
      Tau = 1

      ! ARGUMENTS
      do ia = 1,iargc()
            call getarg(ia,word)
            if(word(1:5).eq.'-Tau=') then
                  read(word,'(5x,i14)') Tau
            else if(word(1:8).eq.'-NSweep=') then
                  read(word,'(8x,i14)') nsweep
                  nburnin = max(10,nsweep/10**3)
            else if(word(1:9).eq.'-NBurnin=') then
                  read(word,'(9x,i14)') nburnin
            end if
      end do

      write(0,*)  '  '
      write(0,*)  '! LOG =============================================='
      write(0,*)  '  '
      write(0,*)  '! ALGORITHM PARAMETERS ============================='
      write(0,*)  ' '
      write(0,*)  'BURN IN AREA             : ', nburnin
      write(0,*)  'NUMBER OF ITERATIONS     : ', nsweep
      write(0,*)  'Tau                      : ', Tau
      write(0,*)  '  '

!     SET THE DATASET
!     ===============

      data en / 191 /
      data (y(i), i = 1,191) / 74,   231,   354,   356,   480, 492,      &
   496,   506,   722,   802,   814,   847,                   &
   913,  1145,  1971,  2011,  2023,  2052,  2242,  2339,  2404,  2590,  2613,  2705,&
  2902,  3333,  3349,  3503,  3598,  3623,  3642,  3720,  3922,  3958,  4068,  4344,&
  4360,  4448,  4673,  4726,  4743,  5281,  5468,  5502,  5603,  5644,  5783,  5825,&
  5826,  6076,  6156,  6159,  6483,  6539,  6570,  6666,  6736,  6777,  6870,  6894,&
  6985,  7128,  7144,  7171,  7315,  7360,  7366,  7574,  7603,  7715,  7758,  7951,&
  8085,  8505,  8600,  8725,  8759,  8886,  9104,  9106,  9106,  9484,  9520,  9535,&
  9566,  9781,  9792,  9929,  9933,  9948, 10020, 10116, 10240, 10290, 10410, 10613,&
 10789, 10844, 10937, 10996, 11311, 11370, 11431, 11432, 11445, 11634, 11979, 11999,&
 12080, 12366, 12480, 12588, 12776, 13009, 13037, 13059, 13120, 13198, 13297, 13623,&
 13898, 13952, 14169, 14282, 14314, 14702, 14853, 15214, 15526, 15880, 16187, 16462,&
 16540, 16557, 17762, 18406, 18873, 19744, 19792, 19915, 20371, 20869, 20918, 21049,&
 21231, 21486, 21680, 21904, 22470, 22932, 23160, 23966, 24483, 26126, 26180, 26506,&
 27818, 28166, 28911, 29128, 29248, 29523, 29543, 29609, 29901, 29905, 30273, 30580,&
 30916, 30935, 31264, 31594, 31906, 32442, 32587, 32662, 33026, 33063, 33082, 33238,&
 33285, 33414, 35044, 35073, 35290, 35297, 35315, 36673, 39039, 39991, 40623/

!     SET THE PARAMETERS OF THE MODEL
!     ===============================

      kmax = 30
      Leng = 40907.d0
      lambda = 3.d0
      cc = 2.d0
      dd = 2.d0
      ee = 1.d0
      ff=en/Leng
      aspr=2.d0

      write(0,*)  '  '
      write(0,*)  'FIXED PARAMETERS ==================================='
      write(0,*)  '  '
      write(0,*)  'Rnage                    : ', 0*Leng, ' to ', Leng
      write(0,*)  'Prior kmax               : ', nburnin
      write(0,*)  'Prior c                  : ', cc
      write(0,*)  'Prior d                  : ', dd
      write(0,*)  'Prior e                  : ', ee
      write(0,*)  'Prior f                  : ', ff
      write(0,*)  '  '

!     OPEN THE FILES
!     =================================

      ! SPECIFICATION
      istdT = 14-int(log10(0.5+Tau))
      write(specT,'(i14)') Tau

      open(1,file='./results/k.T'//specT(istdT:14))
      open(2,file='./results/s.T'//specT(istdT:14))
      open(3,file='./results/h.T'//specT(istdT:14))
      open(4,file='./results/alpha.T'//specT(istdT:14))
      open(5,file='./results/beta.T'//specT(istdT:14))
      open(6,file='./results/AccPrRJ.T'//specT(istdT:14))
      open(7,file='./results/logAIW.T'//specT(istdT:14))
      open(8,file='./results/move.T'//specT(istdT:14))

!     START THE RANDOM NUMBER GENERATOR
!     =================================

      call system_clock(count=seed)
      call init_genrand(seed)
      do i = 1,10; call rnguniform(u); end do

      write(0,*)  '  '
      write(0,*)  'Random number generator ============================'
      write(0,*)  '  '
      write(0,*)  'seed                     : ', seed
      write(0,*)  '  '

!     PRIORS OF THE MODEL
!     ===========================

      pr_k(0) = exp(-lambda)
      do j = 1, kmax
            pr_k(j) = pr_k(j-1)*lambda/j
      end do

!     MODEL PROPOSAL DISTRIBUTION
!     ===========================

      Qtype = 2

      ! ARGUMENTS
      do ia = 1,iargc()
            call getarg(ia,word)
            if(word(1:7).eq.'-Qtype=') then
                  read(word,'(7x,i1)') Qtype
            end if
      end do

      if (Qtype.eq.0) then ! ONLY FIXED MOVES
            do j = 0,kmax
                  pd(j) = 0.d0
                  pb(j) = 0.d0
            end do
      else if (Qtype.eq.1) then ! REASONABLE MOVES
            pb(0) = 1.0d0
            pd(0) = 0.0d0
            do j = 1,kmax-1
                  pb(j) = 0.5d0
                  pd(j) = 0.5d0
            end do
            pb(kmax) = 0.0d0
            pd(kmax) = 1.0d0
      else if (Qtype.eq.2) then ! 'OPTIMISED' MOVES
            cost = min(1.d0/min(1.d0,pr_k(1)/pr_k(0)),1.d0/min(1.d0,pr_k(kmax-1)/pr_k(kmax)))
            do j = 1,kmax-1
                  cost = min(cost,1.d0/(min(1.d0,pr_k(j+1)/pr_k(j))+min(1.d0,pr_k(j-1)/pr_k(j))))
            end do
            cost = 0.9d0*cost
            pd(0) = 0.0
            do j = 0,kmax-1
                  pb(j) = cost*min(1.d0,pr_k(j+1)/pr_k(j))
                  pd(j+1) = cost*min(1.d0,pr_k(j)/pr_k(j+1))
            end do
            pb(kmax) = 0.d0
      else 
            stop
      end if

      write(0,*)  '  '
      write(0,*)  'MODER PROPOSALS ===================================='
      write(0,*)  '  '
      write(0,*)  'Type                     : ', Qtype
      write(0,*)  '  '

! ======================================================================
! SET THE SEEDS
! ======================================================================

      ! RANDOM PARAMETERS

      kseed = 3
      ! ARGUMENTS
      do ia = 1,iargc()
            call getarg(ia,word)
            if(word(1:7).eq.'-Kseed=') then
                  read(word,'(7x,i2)') kseed
            end if
      end do

      k = kseed
      s(0) = 0.d0
      do j = 1,k
            s(j) = j/(k+1.0)*Leng
      end do
      s(k+1) = 1.d0*Leng
      h(0:k) = dble(en/Leng)
      alpha = 1.d0
      beta = 1.d0*Leng/en

      write(0,*)  '  '
      write(0,*)  'SEEDS =============================================='
      write(0,*)  '  '
      write(0,*)  'k                        : ', k
      do j = 0,k
      write(0,'("s(",i3") : ",3X,f15.6,3X,", h(",i3") : ",f15.6)')  j,s(j),j,h(j)
      end do
      write(0,'("s(",i3") : ",3X,f15.6,3X)')  k+1,s(k+1)
      write(0,*)  'alpha                    : ', alpha
      write(0,*)  'beta                     : ', beta

! ======================================================================
! SET THE COUNTERS
! ======================================================================

      countK = 0
      ExpAccPrRJ = 0.d0
      ExpAccPrSplit = 0.d0
      nSplit = 0
      ExpAccPrMerge = 0.d0
      nMerge = 0

! ======================================================================
! SWEEP (BLOCKWISE MCMC, SYSTEMATIC SCAN)
! ======================================================================

      do iter = -nburnin,nsweep

      if(mod(iter,nsweep/10).eq.0) then
            write(0,'(i5,$)') (nsweep-iter)/(nsweep/10)
      end if

      ! ================================================================
      ! BLOCK 1ST      ! FIXED UPDATE (SYSTEMATIC SCAN)
      ! ================================================================

      ! POSITIONS
      if (k.ne.0) then
            call rnguniform(u)
            jup = 1+int(k*u)
            call updateS(k,h(0:k),alpha,beta,s(0:k+1),y(1:en),en,jup)
            move = 'P'
      end if

      ! HEIGHTS
      call rnguniform(u)
      jup = int((k+1)*u)
      call updateH(k,h(0:k),alpha,beta,s(0:k+1),y(1:en),en,jup)
      move = 'H'

      ! ALPHA
      call updateAlpha(k,h(0:k),alpha,cc,dd,aspr,beta)
      move = 'A'

      ! BETA
      call updateBeta(k,h(0:k),alpha,beta,ee,ff)
      move = 'B'

      ! ================================================================
      ! BLOCK 2ND      ! RJ UPDATE (RANDOM SCAN)
      ! ================================================================

      call rnguniform(u)
      if (u.le.pb(k)) then

            ! direction of the move
            k1 = k
            k2 = k+1

      ! BIRTH RJ MOVE
      ! ================================================================

            ! UPDATE TEMPERATURE (t=0)
            t = 0
            call gimel(gt,0,Tau)
            call gimel(Dgt,0+1,Tau)
            Dgt = Dgt-gt

            ! POSITIONS
            call rnguniform(u)
            s0 = u*Leng
            ! JSTAR
            jstar = 0
            do
                  if (s(jstar+1).gt.s0) exit
                  jstar = jstar+1
            end do
            ! HEIGHTS
            call rnguniform(u)

            h1 = exp( &
                        log(h(jstar)) &
                        -(s(jstar+1)-s0) &
                              /(s(jstar+1)-s(jstar)) &
                                    *(log(1.d0-u)-log(u)) &
                        )
            h2 = exp( log(h1) +log(1.d0-u) -log(u) )

            ! REARRANGE

            k_ais = k               ! k
            do j = 0,jstar-1
                  s_ais(j) = s(j)   ! s(0:jstar-1)
                  h_ais(j) = h(j)   ! h(0:jstar-1)
            end do
            s_ais(jstar) = s(jstar) ! s(jstar)
            h_ais(jstar) = h1       ! h(jstar)
            s_ais(jstar+1) = s0     ! s(jstar+1)
            h_ais(jstar+1) = h2     ! h(jstar+1)
            do j = jstar+2,k+1
                  s_ais(j) = s(j-1) ! s(k+2)
                  h_ais(j) = h(j-1) ! h(k+2)
            end do
            s_ais(k+2) = s(k+1)     ! s(k+2)
            alpha_ais = alpha       ! alpha
            beta_ais = beta         ! beta

            ! IW (t=0)
            call AnnImpWeight(logIW,k_ais,lambda, &
                              h_ais(0:k_ais+1),alpha_ais,beta_ais, &
                              s_ais(0:k_ais+2),Leng,jstar, &
                              y(1:en),en)
            ! AIW (t=0)
            logAIW = Dgt*logIW

            ! AIS STEP
            do t = 1, Tau-1

                  ! UPDATE TEMPERATURE (t)
                  call gimel(gt,t,Tau)
                  call gimel(Dgt,t+1,Tau)
                  Dgt = Dgt-gt

                  ! AIS SWEEP (t)
                  call AisRjSplitSweep(k_ais,lambda, &
                                          h_ais(0:k_ais+1), &
                                          alpha_ais,cc,dd,aspr, &
                                          beta_ais,ee,ff, &
                                          s_ais(0:k_ais+2),jstar, &
                                          y(1:en),en,gt)

                  ! IW (t)
                  call AnnImpWeight(logIW,k_ais,lambda, &
                                    h_ais(0:k_ais+1),alpha_ais,beta_ais, &
                                    s_ais(0:k_ais+2),Leng,jstar, &
                                    y(1:en),en)

                  ! AIW (t)
                  logAIW = logAIW+Dgt*logIW

            end do

            ! AISRJ ACCEPTANCE PROBABILITY
            AccPr = min(1.d0,pd(k_ais+1)/pb(k_ais)*exp(logAIW))

            ! AISRJ ACCEPT/REJECT
            call rnguniform(u)
            if (AccPr.ge.u) then
                  k = k_ais+1
                  do j = 0,k
                        s(j) = s_ais(j)
                        h(j) = h_ais(j)
                  end do
                  s(k+1) = s_ais(k+1)
                  alpha = alpha_ais
                  beta = beta_ais
            end if

            ! COUNT THE MOVE
            move = 'S'

      else if (u.le.pb(k)+pd(k)) then

      ! DEATH RJ MOVE
      ! ================================================================

            ! direction of the move
            k1 = k
            k2 = k-1

            ! UPDATE TEMPERATURE (t=0)
            t = Tau
            call gimel(gt,Tau,Tau)
            call gimel(Dgt,Tau-1,Tau)
            Dgt = gt-Dgt

            ! MATCH THE DIMs STEP (t=0)
            ! jstar
            call rnguniform(u)
            jstar = int(k*u)  ! [0:k-1]

            ! REARRANGE (t=0)
            k_ais = k-1
            do j = 0,k
                  s_ais(j) = s(j)
                  h_ais(j) = h(j)
            end do
            s_ais(k+1) = s(k+1)
            alpha_ais = alpha
            beta_ais = beta

            ! IW  (t=0)
            call AnnImpWeight(logIW,k_ais,lambda, &
                              h_ais(0:k_ais+1),alpha_ais,beta_ais, &
                              s_ais(0:k_ais+2),Leng,jstar, &
                              y(1:en),en)

            ! AIW (t=0)
            logAIW = -Dgt*logIW

            ! AIS STEP
            do t = Tau-1,1,-1

                  ! UPDATE TEMPERATURE (t)
                  call gimel(gt,t,Tau)
                  call gimel(Dgt,t-1,Tau)
                  Dgt = gt-Dgt

                  ! AIS SWEEP (t)
                  call AisRjSplitSweep(k_ais,lambda, &
                                          h_ais(0:k_ais+1), &
                                          alpha_ais,cc,dd,aspr, &
                                          beta_ais,ee,ff, &
                                          s_ais(0:k_ais+2),jstar, &
                                          y(1:en),en,gt)

                  ! IW (t)
                  call AnnImpWeight(logIW,k_ais,lambda,&
                                    h_ais(0:k_ais+1),alpha_ais,beta_ais,&
                                    s_ais(0:k_ais+2),Leng,jstar, &
                                    y(1:en),en)
                  ! AIW (t)
                  logAIW = logAIW-Dgt*logIW

            end do

            ! AISRJ ACCEPTANCE PROBABILITY
            AccPr = min( 1.d0, pb(k_ais)/pd(k_ais+1) *exp(logAIW) )

            call rnguniform(u)
            if (AccPr.ge.u) then
                  k = k_ais
                  do j = 0,jstar-1
                        s(j) = s_ais(j)
                        h(j) = h_ais(j)
                  end do
                  s(jstar) = s_ais(jstar)
                  h(jstar) = exp( &
                                    (s_ais(jstar+1)-s_ais(jstar)) &
                                          /(s_ais(jstar+2)-s_ais(jstar)) &
                                          *log(h_ais(jstar)) &
                                    +(s_ais(jstar+2)-s_ais(jstar+1)) &
                                          /(s_ais(jstar+2)-s_ais(jstar)) &
                                          *log(h_ais(jstar+1)) &
                                    )
                  do j = jstar+1,k
                        s(j) = s_ais(j+1)
                        h(j) = h_ais(j+1)
                  end do
                  s(k+1) = s_ais(k+2)
                  ! CLEAN
                  s(k+2) = 0.d0
                  h(k+1) = 0.d0
                  alpha = alpha_ais
                  beta = beta_ais
            end if

            move = 'M'

      else

            ! direction of the move
            k1 = k
            k2 = k

            move = 'F'

      end if

      ! SAVE THE SAMPLE
      ! ================================================================

      if (iter.gt.0) then

            ! K
            write(1,*) k
            ! S
            write(2,100) s(0:k+1)
            ! H
            write(3,101) h(0:k)
            ! ALPHA
            write(4,*) alpha
            ! BETA
            write(5,*) beta
            ! ACCEPTACE PROBABILITIES
            write(6,*) AccPr
            ! LOG SCALE ANNEALING IMPORTANCE WEIGHTS
            write(7,*) logAIW
            ! MOVE
            write(8,*) k1, k2
                  
            ! COUNT ACCEPTACE PROBABILITIES
            if (move.eq.'S'.or.move.eq.'M') ExpAccPrRJ = ExpAccPrRJ+AccPr
            if (move.eq.'S') then
                  ExpAccPrSplit = ExpAccPrSplit+AccPr
                  nSplit = nSplit+1
            else if (move.eq.'M') then
                  ExpAccPrMerge = ExpAccPrMerge+AccPr
                  nMerge = nMerge+1
            end if

            ! COUNT K
            countK = countK+k

      end if

      end do

      ! COMPUTE THE EXPECTATIONS
      ExpAccPrRJ = ExpAccPrRJ/(nSplit+nMerge)
      ExpAccPrSplit = ExpAccPrSplit/nSplit
      ExpAccPrMerge = ExpAccPrMerge/nMerge

      ! PRINT OUTPUT
      write(0,*)  '  '
      write(0,*)  '! ESTIMATES =============================='
      write(0,*)  ' '
      write(0,*)  'Estimated average k                       : ', countK/dble(nsweep)
      write(0,*)  ' '
      write(0,*)  'Estimated Acceptance probability          : ', ExpAccPrRJ
      write(0,*)  'Estimated Acceptance probability (Split)  : ', ExpAccPrSplit
      write(0,*)  'Estimated Acceptance probability (Merge)  : ', ExpAccPrMerge
      write(0,*)  ' '

      write(0,*)  'DONE'

100      format(30f15.5)
101      format(30f15.10)

end program aisrjcpt_pro


