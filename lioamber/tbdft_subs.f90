#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module tbdft_subs
   implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_init(M_in, Nuc, open_shell)
! This subroutine initialize the variables for TBDFT calculations. Also the
! file gamma.in is readed.
   use tbdft_data, only: MTB, MTBDFT, end_bTB, Iend_TB, rhoa_TBDFT, rhob_TBDFT,&
                         gammaW, n_biasTB, basTB, n_atTB,n_atperbias,linkTB,   &
                         VbiasTB, rhofirst_TB,tbdft_calc
   implicit none
   logical, intent(in)       :: open_shell
   integer, intent(in)       :: M_in
   integer, intent(in)       :: Nuc(M_in)
   LIODBLE, allocatable :: rhoTB_real(:,:,:)
   integer                   :: tot_at
   integer                   ::  ii, jj,kk,ll,pp,rr
   TDCOMPLEX :: liocmplx

   MTBDFT = MTB+M_in

   allocate(rhoa_TBDFT(MTBDFT,MTBDFT))
   if (open_shell) allocate (rhob_TBDFT(MTBDFT,MTBDFT))

   open(unit = 1001, file = 'gamma.in')

   n_atTB = int(MTB / n_biasTB)

   if (mod(MTB,n_biasTB) /= 0) then
      print*,"MTB most be multiple of the number of bias used"
      stop
   endif

   if (mod(MTB,2) /= 0) then
      print*,"MTB most be multiple of 2"
      stop
   endif

   allocate(VbiasTB(n_biasTB))
   read(1001,*) VbiasTB       !Reading the potential applied to each bias
   read(1001,*) n_atperbias   !Reading number of atoms coupledto one bias
   read(1001,*) end_bTB       !Number of basis coupled per atom
   tot_at=n_biasTB*n_atperbias
   allocate(linkTB(n_biasTB, n_atperbias))
   allocate(gammaW(n_atperbias*end_bTB))
   allocate(basTB(end_bTB))
   allocate(Iend_TB(n_biasTB, end_bTB*n_atperbias))

   do ii = 1, n_biasTB
      read(1001,*) linkTB(ii,:) ! Reading the atoms coupled
   enddo

   read(1001,*) basTB           ! Reading the basis coupled in the LIO order

   do ii = 1, end_bTB * n_atperbias
      read(1001,*) gammaW(ii)   ! Reading the weight of gamma per basis.
   enddo
   close(1001)

  !The index are stored into Iend_TB. First we look for the coupling electrode,
  !then we look for the index which belong to those atoms in the reading order.
   do jj = 1, n_biasTB
      rr = 0
   do kk = 1, n_atperbias
      ll = 0
      pp = 1
      do ii = 1, M_in
         if (linkTB(jj,kk) == Nuc(ii)) then
            ll = ll +1
            if (ll == basTB(pp)) then
               pp = pp +1
               rr = rr +1
               Iend_TB(jj,rr) = ii
            endif
         endif
      enddo
   enddo
   enddo

   if (tbdft_calc == 3) then
      if (open_shell) then
         allocate(rhoTB_real(MTBDFT,MTBDFT,2), rhofirst_TB(MTBDFT,MTBDFT,2))
      else
         allocate(rhoTB_real(MTBDFT,MTBDFT,1), rhofirst_TB(MTBDFT,MTBDFT,1))
      endif

      open(unit = 1001, file = 'rhofirstTB')

      if (open_shell) then
         do jj = 1, MTBDFT
         do ii = 1, MTBDFT
            read(1001,*) rhoTB_real(ii,jj,1), rhoTB_real(ii,jj,2)
         enddo
         enddo
      else
         do jj = 1, MTBDFT
         do ii = 1, MTBDFT
            read(1001,*) rhoTB_real(ii,jj,1)
         enddo
         enddo
      endif

      close(1001)

      do jj = 1,MTBDFT
      do ii = 1,MTBDFT
         if (ii == jj) then
            rhofirst_TB(ii,jj,1) = liocmplx(rhoTB_real(ii,jj,1), 0.0d0)
         else
            rhofirst_TB(ii,jj,1) = liocmplx(0.5d0 * rhoTB_real(ii,jj,1), 0.0d0)
         endif
      enddo
      enddo

      if (open_shell) then
         do jj = 1,MTBDFT
         do ii = 1,MTBDFT
            if (ii == jj) then
               rhofirst_TB(ii,jj,2) = liocmplx(rhoTB_real(ii,jj,2), 0.0d0)
            else
               rhofirst_TB(ii,jj,2) = liocmplx(0.5d0 * rhoTB_real(ii,jj,2), 0.0d0)
            endif
         enddo
         enddo
      endif
      write(*,*) "RHOFIRST_TB READ"
   endif

end subroutine tbdft_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_td_init (M_in,rho, rho_0, thrddim)
   ! This subroutine initialize matrices to store rho matrix during TD.
   use tbdft_data, only: MTB, MTBDFT, rhoa_TBDFT, rhob_TBDFT, rhold_AOTB,  &
                         rhonew_AOTB
   implicit none
   integer  , intent(in)  :: M_in
   integer  , intent(in)  :: thrddim
   TDCOMPLEX, intent(in)  :: rho_0(M_in,M_in,thrddim)
   TDCOMPLEX, intent(out) :: rho(MTBDFT,MTBDFT,thrddim)
   integer :: ii, jj
   TDCOMPLEX :: liocmplx

   allocate(rhold_AOTB(MTBDFT,MTBDFT,thrddim), &
            rhonew_AOTB(MTBDFT,MTBDFT,thrddim))

   rhold_AOTB  = liocmplx(0.0d0,0.0D0)
   rhonew_AOTB = liocmplx(0.0d0,0.0D0)

   do jj = 1, MTBDFT
   do ii = 1, MTBDFT
      if (ii == jj) then
         rho(ii,jj,1) = liocmplx(rhoa_TBDFT(ii,jj), 0.0D0)
      else
         rho(ii,jj,1) = liocmplx(0.5D0 * rhoa_TBDFT(ii,jj), 0.0D0)
      endif
   enddo
   enddo

   ! Open shell option
   if (thrddim == 2) then
      do jj = 1, MTBDFT
      do ii = 1, MTBDFT
         if (ii == jj) then
            rho(ii,jj,2) = liocmplx(rhob_TBDFT(ii,jj), 0.0D0)
         else
            rho(ii,jj,2) = liocmplx(0.5D0 * rhob_TBDFT(ii,jj), 0.0D0)
         endif
      enddo
      enddo
   endif

   do jj = 1, M_in
   do ii = 1, M_in
      rho(ii+MTB,jj+MTB,:) = rho_0(ii,jj,:)
   enddo
   enddo

end subroutine tbdft_td_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine getXY_TBDFT(M_in,x_in,y_in,xmat,ymat)
! This subroutine modified the X_in and Y_in with the TB elements, just if
! tbdft_calc /=0. If not, it doesn't change these matrices.

   use tbdft_data, only: MTB, tbdft_calc

   implicit none
   integer     , intent(in)  :: M_in
   LIODBLE, intent(in)  :: x_in(M_in,M_in)
   LIODBLE, intent(in)  :: y_in(M_in,M_in)
   LIODBLE, intent(out) :: xmat(M_in+MTB,M_in+MTB)
   LIODBLE, intent(out) :: ymat(M_in+MTB,M_in+MTB)
   integer :: ii, jj

   xmat = 0.0d0
   ymat = 0.0d0

   if (tbdft_calc == 0) then
      xmat = x_in
      ymat = y_in
   else
      do ii = 1, MTB
         xmat(ii,ii) = 1.0d0
         ymat(ii,ii) = 1.0d0
      enddo

      do jj = 1, M_in
      do ii = 1, M_in
         xmat(MTB+ii, MTB+jj) = x_in(ii,jj)
         ymat(MTB+ii, MTB+jj) = y_in(ii,jj)
      enddo
      enddo
   endif

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine construct_rhoTBDFT(M, rho, rho_0 ,rho_TBDFT, niter, open_shell)
! This subroutine initialize the TBDFT density matrix in the first step. After
! that, it just correct the TB section of the density (divide it by 2).

   use tbdft_data , only: MTB, MTBDFT

   implicit none
   logical     , intent(in)  :: open_shell
   integer     , intent(in)  :: M
   integer     , intent(in)  :: niter
   LIODBLE, intent(in)  :: rho_0(M,M)
   LIODBLE, intent(in)  :: rho_TBDFT(MTBDFT,MTBDFT)
   LIODBLE, intent(out) :: rho(MTBDFT,MTBDFT)
   integer      :: ii, jj
   LIODBLE :: ocup

   ocup = 1.0d0
   if (open_shell) ocup = 0.5d0
   if (niter /= 1) then

      do ii = 1   , MTBDFT
      do jj = ii+1, MTBDFT
         rho(ii,jj) = 0.5D0 * rho_TBDFT(ii,jj)
         rho(jj,ii) = rho(ii,jj)
      enddo
      enddo

   else if (niter == 1) then
      rho = 0.0D0
      do ii = 1, MTB
         rho(ii,ii) = ocup
      enddo
      rho(MTB+1:MTB+M,MTB+1:MTB+M) = rho_0
   endif

end subroutine construct_rhoTBDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine build_chimera_TBDFT (M_in, fock_in, fock_TBDFT, spin_id)
   ! This subroutine adds the TB elements to the Hamiltonian.
   use tbdft_data, only: MTBDFT, MTB, Iend_TB, end_bTB, alfaTB, alfaTB2,       &
                         betaTB, gammaTB, gammaW, n_biasTB, n_atperbias,       &
                         n_atTB, VbiasTB, tbdft_calc

   integer     , intent(in)  :: M_in
   integer     , intent(in)  :: spin_id
   LIODBLE, intent(in)  :: fock_in (M_in, M_in)
   LIODBLE, intent(out) :: fock_TBDFT (MTBDFT, MTBDFT)
   LIODBLE :: V_aux, fermi_level
   integer      :: ii, jj, kk, link

   fock_TBDFT(:,:) = 0.0D0

   if (spin_id == 1) fermi_level = alfaTB
   if (spin_id == 2) fermi_level = alfaTB2

   V_aux = 1.0D0
   if (tbdft_calc == 1) V_aux = 0.0d0

   do jj = 1, end_bTB * n_atperbias
   do ii = 1, n_biasTB
      link = n_atTB * ii
      fock_TBDFT(link, Iend_TB(ii,jj)+MTB) = gammaW(jj) * gammaTB
      fock_TBDFT(Iend_TB(ii,jj)+MTB,link)  = gammaW(jj) * gammaTB
   enddo
   enddo

   do ii = 1, n_biasTB
   do jj = 1, n_atTB
      kk =jj + ((ii-1) * (n_atTB))
      fock_TBDFT(kk,kk) = fermi_level + V_aux * VbiasTB(ii)
      if (jj < n_atTB) then
         fock_TBDFT(kk,kk+1) = betaTB
         fock_TBDFT(kk+1,kk) = betaTB
      endif
   enddo
   enddo

   fock_TBDFT(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)

end subroutine build_chimera_TBDFT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine extract_rhoDFT (M_in, rho_in, rho_out)
   ! This subroutine separate the DFT part of the density from the TBDFT density
   ! matrix.

   use tbdft_data, only: MTBDFT, MTB

   implicit none
   integer     , intent(in)  :: M_in
   LIODBLE, intent(in)  :: rho_in(MTBDFT,MTBDFT)
   LIODBLE, intent(out) :: rho_out(M_in,M_in)
   integer :: ii, jj

   rho_out=0.0D0

   do jj = 1, M_in
   do ii = 1, M_in
      rho_out(ii,jj) = rho_in(MTB+ii,MTB+jj)
   enddo
   enddo

end subroutine extract_rhoDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine chimeraTBDFT_evol(M_in, fock_in, fock_TBDFT, istep, spin_id)
   ! This subroutine modify and add the TB section of the Hamiltonian during TD.

   use tbdft_data, only: MTBDFT, MTB, Iend_TB, end_bTB, alfaTB, alfaTB2,       &
                         betaTB, gammaTB, start_tdtb, end_tdtb, gammaW,        &
                         n_atTB,n_biasTB, n_atperbias, VbiasTB, tbdft_calc

   integer     , intent(in)  :: M_in
   integer     , intent(in)  :: istep, spin_id
   LIODBLE, intent(in)  :: fock_in(M_in, M_in)
   LIODBLE, intent(out) :: fock_TBDFT(MTBDFT, MTBDFT) ! Temporary dimensions
   LIODBLE :: pi = 4.0D0 * atan(1.0D0)
   LIODBLE :: lambda, t_step, f_t, fermi_level
   integer      :: ii,jj,kk, link

   if (spin_id == 1) fermi_level = alfaTB
   if (spin_id == 2) fermi_level = alfaTB2

   lambda = 1.0d0 / real(end_tdtb - start_tdtb)
   f_t = 1.0D0
   if (tbdft_calc == 1) then
      if (istep < start_tdtb) then
         f_t = 0.0D0
      else if ((istep >= start_tdtb) .and. (istep < end_tdtb)) then
         t_step = real(istep - start_tdtb)
         f_t    = (-cos(pi * lambda * t_step) + 1.0D0) / 2.0D0
      else if (istep >= end_tdtb) then
         f_t = 1.0D0
      endif
   endif
   fock_TBDFT(:,:) = 0.0D0

   do jj = 1, end_bTB * n_atperbias
   do ii = 1, n_biasTB
      link = n_atTB * ii
      fock_TBDFT(link, Iend_TB(ii,jj)+MTB) = gammaW(jj) * gammaTB
      fock_TBDFT(Iend_TB(ii,jj)+MTB,link)  = gammaW(jj) * gammaTB
   enddo
   enddo

   do ii = 1, n_biasTB
   do jj = 1, n_atTB
      kk = jj + ((ii -1) * (n_atTB))
      fock_TBDFT(kk,kk) = fermi_level + f_t * VbiasTB(ii)
      if (jj < n_atTB) then
         fock_TBDFT(kk,kk+1) = betaTB
         fock_TBDFT(kk+1,kk) = betaTB
      endif
   enddo
   enddo

   fock_TBDFT(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)

end subroutine chimeraTBDFT_evol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine TB_current (M_in,delta_rho, overlap, TB_electrode, TB_M)
   ! This subroutine calculates the charge difference in each TB electrode and the
   ! DFT part.
   use tbdft_data, only:MTBDFT, MTB,n_atTB,n_biasTB

   implicit none
   integer     , intent(in)    :: M_in
   LIODBLE, intent(in)    :: overlap(M_in, M_in)
   LIODBLE, intent(in)    :: delta_rho(MTBDFT,MTBDFT)
   LIODBLE, intent(inout) :: TB_M
   LIODBLE, intent(inout) :: TB_electrode(n_biasTB)
   integer      :: ii, jj, kk
   LIODBLE :: qe

   do ii = 1, n_biasTB
   do jj = 1, n_atTB
      kk = jj + ((ii-1) * (n_atTB))
      TB_electrode(ii) = TB_electrode(ii) + delta_rho(kk,kk)
   enddo
   enddo

   do ii = 1, M_in
   do jj = 1, M_in
      qe = delta_rho(MTB+ii,MTB+jj) * overlap(ii, jj)
      TB_M = qe + TB_M
   enddo
   enddo

end subroutine TB_current

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_scf_output(open_shell)
  ! This subroutine calculate the charge of the TB electrodes after a SCF
  ! calculation.

   use tbdft_data, only: rhoa_TBDFT, rhob_TBDFT, MTBDFT, n_biasTB, n_atTB,&
                         tbdft_calc

   implicit none
   logical, intent(in) :: open_shell
   LIODBLE        :: rho_aux(MTBDFT, MTBDFT)
   LIODBLE        :: chargeTB(n_biasTB), chargeTB_spin(n_biasTB)
   integer             :: ii, jj,kk

   if (tbdft_calc == 0) return

   if (open_shell) then
      rho_aux = rhoa_TBDFT + rhob_TBDFT
      chargeTB_spin = 0.0d0
   else
      rho_aux = rhoa_TBDFT
   endif

   chargeTB = n_atTB

   do ii = 1, n_biasTB
   do jj = 1, n_atTB
      kk = jj + ((ii-1) * (n_atTB))
      chargeTB(ii) = chargeTB(ii) - rho_aux(kk,kk)
      if (open_shell) chargeTB_spin(ii) = chargeTB_spin(ii) +                  &
                                          rhob_TBDFT(kk,kk) - rhoa_TBDFT(kk,kk)
   enddo
   enddo

   open(unit = 20202,  file = 'mullikenTB')
   if (open_shell) then
      do ii = 1, n_biasTB
         write(20202,*) "Mulliken Spin TB  electrode", ii, chargeTB_spin(ii)
      enddo
   endif
   do ii = 1, n_biasTB
      write(20202,*) "Mulliken TB  electrode", ii, chargeTB(ii)
   enddo
   close(20202)
end subroutine tbdft_scf_output
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_td_output(M_in, thrddim, rho_aux, overlap, istep, Iz, natom, &
                           Nuc, open_shell)
! This subroutine store the charge and charge differnce of each TD step into
! currentTB, and mullikenTB files.

   use tbdft_data, only: rhold_AOTB, rhonew_AOTB, MTB, MTBDFT,n_atTB, n_biasTB,&
                         tbdft_calc
   use properties, only: mulliken
   implicit none

   logical     , intent(in) :: open_shell
   integer     , intent(in) :: M_in, istep, thrddim
   integer     , intent(in) :: natom
   integer     , intent(in) :: Nuc(M_in)
   integer     , intent(in) :: Iz(natom)
   LIODBLE, intent(in) :: overlap(M_in,M_in)
   TDCOMPLEX   , intent(in) :: rho_aux(MTBDFT, MTBDFT,thrddim)
   integer      :: ii, jj, kk
   LIODBLE :: I_TB_elec(n_biasTB), I_TB_elec_a(n_biasTB)
   LIODBLE :: I_TB_M
   LIODBLE :: chargeTB(n_biasTB), chargeTB_a(n_biasTB)
   LIODBLE :: chargeM_TB, chargeM_TB_a
   LIODBLE :: qe(natom), qs(natom)
   LIODBLE :: rhoscratch(MTBDFT,MTBDFT,thrddim)

   if (tbdft_calc == 0) return
   I_TB_M    = 0.0d0
   I_TB_elec = 0.0d0

   if (istep == 1) then
      open(unit = 10101, file = 'currentTB')
      open(unit = 20202, file = 'mullikenTB')

   else
      if (tbdft_calc == 1) then
         rhoscratch = real(rhonew_AOTB - rhold_AOTB)
      else if (tbdft_calc > 1) then
         rhoscratch = real(rhonew_AOTB)
      endif

      call TB_current(M_in,rhoscratch(:,:,1), overlap, I_TB_elec, I_TB_M)
      if (open_shell) then
         I_TB_elec_a = I_TB_elec
         call TB_current(M_in,rhoscratch(:,:,2), overlap, I_TB_elec, I_TB_M)

         do ii = 1, n_biasTB
            write(10101,*) "Current Alpha TB electrode", ii, I_TB_elec_a(ii)
            write(10101,*) "Current Beta  TB electrode", ii, I_TB_elec(ii) -   &
                                                             I_TB_elec_a(ii)
         enddo
      endif

      do ii = 1, n_biasTB
         write(10101,*) "Current TB electrode", ii, I_TB_elec(ii)
      enddo
      write(10101,*) "Current DFT part M", I_TB_M

      chargeTB = n_atTB
      do ii = 1, n_biasTB
      do jj = 1, n_atTB
         kk = jj + ((ii-1) * (n_atTB))
         chargeTB(ii) = chargeTB(ii) - real(rho_aux(kk,kk,1),COMPLEX_SIZE)
      enddo
      enddo

      if (open_shell) then
         chargeTB_a = 0.0d0
         do ii = 1, n_biasTB
         do jj = 1, n_atTB
            kk = jj + ((ii-1) * (n_atTB))
            chargeTB_a(ii) = chargeTB_a(ii)-real(rho_aux(kk,kk,2),COMPLEX_SIZE)
         enddo
         enddo

         do ii = 1, n_biasTB
              write(20202,*) "Mulliken Spin TB electrode", ii,-chargeTB_a(ii)+ &
                                                          chargeTB(ii) - n_atTB
          enddo
         chargeTB = chargeTB + chargeTB_a

      endif

      chargeM_TB = 0.0D0
      rhoscratch = real(rho_aux)

      if (open_shell) then
         call mulliken(rhoscratch(MTB+1:MTB+M_in,MTB+1:MTB+M_in,1), &
                       rhoscratch(MTB+1:MTB+M_in,MTB+1:MTB+M_in,2), &
                       overlap, Nuc, Iz, qe, qs)

         chargeM_TB_a = 0.0d0
         do ii = 1, natom
            chargeM_TB_a = chargeM_TB_a + qs(ii)
         enddo

         write(20202,*) "Mulliken Spin DFT  part M", chargeM_TB_a
      else
         call mulliken(rhoscratch(MTB+1:MTB+M_in,MTB+1:MTB+M_in,1), &
                       overlap, Nuc, Iz, qe)
      endif

      do ii = 1, natom
            chargeM_TB = chargeM_TB + qe(ii)
      enddo

      do ii = 1, n_biasTB
         write(20202,*) "Mulliken TB electrode", ii, chargeTB(ii)
      enddo
      write(20202,*) "Mulliken DFT  part M", chargeM_TB

   endif
end subroutine tbdft_td_output
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine write_rhofirstTB(M_in, OPEN)
   ! This subroutine write the rho matrix after a SCF in rhofirstTB file,
   ! for being
   ! used wtith TB-DLVN.
   use tbdft_data , only: tbdft_calc, rhoa_TBDFT, rhob_TBDFT

   implicit none
   integer, intent(in) :: M_in
   logical, intent(in) :: OPEN
   integer :: ii, jj

   if (tbdft_calc /= 2) return

   open(unit = 1001, file = 'rhofirstTB')
   if (OPEN) then
      do jj = 1, M_in
      do ii = 1, M_in
         write(1001,*) rhoa_TBDFT(ii,jj), rhob_TBDFT(ii,jj)
      enddo
      enddo
   else
      do jj = 1, M_in
      do ii = 1, M_in
         write(1001,*) rhoa_TBDFT(ii,jj)
      enddo
      enddo
   endif

   close(1001)

   write(*,*) 'RHOFIRST_TB writed'

end subroutine write_rhofirstTB
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine transport_TB(M, dim3, rho_aux ,Ymat, istep, OPEN, rho_aop, rho_bop)
   ! This subroutine add the driving term in TD during a DLVN calculation.
   use tbdft_data      , only: MTB, MTBDFT, rhofirst_TB, driving_rateTB,       &
                               tbdft_calc, rhonew_AOTB
   use typedef_cumat   , only: cumat_x
   use typedef_operator, only: operator

   implicit none
   logical        , intent(in)  :: OPEN
   integer        , intent(in)  :: M, dim3, istep
   type(cumat_x)  , intent(in)  :: Ymat
   TDCOMPLEX      , intent(out) :: rho_aux(MTBDFT,MTBDFT,dim3)
   type (operator), intent(in)  :: rho_aop
   type (operator), intent(in), optional :: rho_bop

   LIODBLE    :: rho_real(MTBDFT,MTBDFT,dim3)
   LIODBLE    :: scratchgamma
   TDCOMPLEX       :: rhoscratch(MTBDFT,MTBDFT,dim3)
   TDCOMPLEX       :: liocmplx
   integer         :: ii, jj

   if (tbdft_calc /= 3) return

   rhoscratch   = liocmplx(0.0D0,0.0D0)
   rho_real     = 0.0d0
   scratchgamma = 0.0d0

   call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
   if (OPEN) call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))

   if (istep>200.and.istep<=1200) then
      scratchgamma = driving_rateTB * &
                     dexp(-0.0001D0 * (dble(istep - 1200)**2))
   else if (istep <= 200) then
      scratchgamma = 0.0d0
   else if (istep > 1200) then
      scratchgamma = driving_rateTB
   endif

   do jj = 1, MTB
   do ii = 1, MTB
      rhoscratch(ii,jj,1) = rho_aux(ii,jj,1) - rhofirst_TB(ii,jj,1)
      if (OPEN) then
         rhoscratch(ii,jj,2) = rho_aux(ii,jj,2) - rhofirst_TB(ii,jj,2)
      endif
   enddo
   enddo

   do jj = MTB+1, MTB+M
   do ii = 1, MTB
      rhoscratch(ii,jj,1) = real(0.5d0,COMPLEX_SIZE/2) * &
                            (rho_aux(ii,jj,1) - rhofirst_TB(ii,jj,1))
      rhoscratch(jj,ii,1) = real(0.5d0,COMPLEX_SIZE/2) * &
                            (rho_aux(jj,ii,1) - rhofirst_TB(jj,ii,1))
      if (OPEN) then
         rhoscratch(ii,jj,2) =  real(0.5d0,COMPLEX_SIZE/2) * &
                                (rho_aux(ii,jj,2) - rhofirst_TB(ii,jj,2))
         rhoscratch(jj,ii,2) =  real(0.5d0,COMPLEX_SIZE/2) * &
                                (rho_aux(jj,ii,2) - rhofirst_TB(jj,ii,2))
      endif
   enddo
   enddo

   rhoscratch  = real(scratchgamma,COMPLEX_SIZE/2) * rhoscratch
   rhonew_AOTB = rhoscratch
   rho_aux     = rhoscratch

   call Ymat%change_base(rho_aux(:,:,1), 'dir')
   if (OPEN) call Ymat%change_base(rho_aux(:,:,2), 'dir')

end subroutine transport_TB
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_calibration(E, fock_aop, rho_aop, fock_bop, rho_bop)
   ! This subroutine calculate iteratively the fermi energy in TB electrodes,
   ! necessary to concentrate a charge equal to TB_charge_ref in the DFT part.

   use garcha_mod      , only: Smat, Iz, natom, Pmat_vec
   use basis_data      , only: M, Nuc
   use SCF_aux         , only: fix_densmat
   use tbdft_data      , only: alfaTB, TB_q_tot, TB_charge_ref, TB_q_told
   use typedef_operator, only: operator
   use properties      , only: mulliken

   implicit none
   type(operator), intent(inout)           :: rho_aop, fock_aop
   type(operator), intent(inout), optional :: rho_bop, fock_bop
   LIODBLE, intent(inout) :: E
   LIODBLE :: Q_old, Q_new, delta_Q
   LIODBLE :: Ef_old, Ef_new
   LIODBLE :: q(natom)
   LIODBLE :: escale_f
   LIODBLE :: rho_aux(M,M)
   integer :: niter, ii
   logical :: converged=.false.

   escale_f=1.0d0

   write(*,'(A)') "INITIATING TB CHARGE CONVERGENCY"

   niter = 0

   ! Initializations.
   Q_old = 0.0d0
   Ef_old = alfaTB

   do while (.not.converged.and.niter<=TB_q_tot)
      niter = niter + 1

      Ef_new = alfaTB

      call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
      call spunpack('L', M, Pmat_vec, rho_aux)
      call fix_densmat(rho_aux)

      call mulliken(rho_aux, Smat, Nuc, Iz, q)

      Q_new = 0.0d0
      do ii=1,natom
         Q_new = Q_new + q(ii)
      enddo

      if (niter == 1) Q_old = Q_new

      if ((Q_new-TB_charge_ref>0.0d0.and.Q_old-TB_charge_ref<0.0d0) .or. &
          (Q_new-TB_charge_ref<0.0d0.and.Q_old-TB_charge_ref>0.0d0)) then
         escale_f = escale_f * 0.1d0
         Q_new  = Q_old
         Ef_new = Ef_old
      end if

      delta_Q = Q_new - TB_charge_ref

      if (delta_Q > 0.0d0) then
         alfaTB = Ef_new + 0.1d0 * escale_f
      else if (delta_Q < 0.0d0) then
         alfaTB = Ef_new - 0.1d0 * escale_f
      end if

      if (abs(delta_Q)<TB_q_told) converged = .true.

      write(*,'(A)') "---------------------------------------------------------"
      write(*,'(A,I4)') "TB Charge Convergence Step =", niter
      write(*,'(A,F16.10,A)') "Fermi level = ", Ef_new, " A.U."
      write(*,'(A,F16.10,A)') "New charge = ", Q_new, " A.U."
      Write(*,'(A,ES9.2,A,ES9.2)') "Charge difference = ",delta_Q,             &
                                   " Tolerance = ", TB_q_told
      write(*,'(A)') "---------------------------------------------------------"

      Q_old  = Q_new
      Ef_old = Ef_new

   end do

      if (converged) then
         write(*,'(A)')"THE CONVERGENCE HAS FINISHED SUCCESFULLY"
         write(*,'(A)')"----------------------------------------"
         write(*,'(A,F16.10,A)')"Best Fermi energy = ", Ef_new, " A.U."
         write(*,'(A,F16.10,A)')"Best charge = ", Q_new, " A.U."
      else
         write(*,*)"NO CONVERGNCE ACHIEVED, MORE STEPS ARE NEEDED"
         write(*,*)"------------------------------------"
         write(*,'(A,F16.10,A)')"Best Fermi energy =", Ef_new, " A.U."
         write(*,'(A,F16.10,A)')"Best charge =", Q_new, " A.U."
      end if

end subroutine tbdft_calibration

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module tbdft_subs
