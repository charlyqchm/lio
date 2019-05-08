!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module tbdft_subs
   implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_init(M_in, Nuc, natom, open_shell)

   use tbdft_data, only: MTB, MTBDFT, end_bTB, Iend_TB, rhoa_TBDFT, rhob_TBDFT,&
                         gammaW, n_biasTB, basTB, n_atTB,n_atperbias,linkTB,   &
                         VbiasTB

   implicit none

   logical, intent(in) :: open_shell
   integer, intent(in) :: M_in, natom
   integer, intent(in) :: Nuc(M_in)
   integer :: tot_at
   integer ::  ii, jj,kk,ll,pp,rr

   MTBDFT = MTB+M_in

   allocate(rhoa_TBDFT(MTBDFT,MTBDFT))
   if (open_shell) allocate (rhob_TBDFT(MTBDFT,MTBDFT))


   open(unit=1001, file='gamma.in')

   read(1001,*) n_biasTB

   n_atTB = int(MTB/n_biasTB)

   if(mod(MTB,n_biasTB)/=0) then
      print*,"MTB most be multiple of the number of bias used"
      stop
   end if

   allocate(VbiasTB(n_biasTB))
   read(1001,*) VbiasTB
   read(1001,*) n_atperbias
   read(1001,*) end_bTB
   tot_at=n_biasTB*n_atperbias
   allocate(linkTB(n_biasTB, n_atperbias))
   allocate(gammaW(n_atperbias*end_bTB))
   allocate(basTB(end_bTB))
   allocate(Iend_TB(n_biasTB, end_bTB*n_atperbias))

   do ii=1, n_biasTB
      read(1001,*) linkTB(ii,:)
   end do

   read(1001,*) basTB

   do ii = 1, end_bTB*n_atperbias
      read(1001,*) gammaW(ii)
   enddo
   close(1001)

!TB: Los indices que se guardan en Iend_TB se buscan primero por el electrodo
!    al que se acoplan, luego buscan los indices que corresponden dichos atomos
!    en el orden que fueron leidos.

   do jj = 1, n_biasTB
      rr=0
   do kk = 1, n_atperbias
      ll=0
      pp=1
      do ii = 1, M_in
         if (linkTB(jj,kk)==Nuc(ii)) then
            ll=ll+1
            if(ll==basTB(pp)) then
               pp=pp+1
               rr=rr+1
               Iend_TB(jj,rr) = ii
            end if
         end if
      end do
   end do
   end do

end subroutine tbdft_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_td_init (M_in,rho, rho_0, thrddim)

   use tbdft_data, only: MTB, MTBDFT, rhoa_TBDFT, rhob_TBDFT, rhold_AOTB,  &
                         rhonew_AOTB
   implicit none
   integer        , intent(in) :: M_in
   integer        , intent(in) :: thrddim
#ifdef TD_SIMPLE
   complex(kind=4), intent(in)  :: rho_0(M_in,M_in,thrddim)
   complex(kind=4), intent(out) :: rho(MTBDFT,MTBDFT,thrddim)
#else
   complex(kind=8), intent(in)  :: rho_0(M_in,M_in,thrddim)
   complex(kind=8), intent(out) :: rho(MTBDFT,MTBDFT,thrddim)
#endif
   integer :: ii, jj

   allocate(rhold_AOTB(MTBDFT,MTBDFT,thrddim), &
            rhonew_AOTB(MTBDFT,MTBDFT,thrddim))

   rhold_AOTB  = 0.0d0
   rhonew_AOTB = 0.0d0

   do jj = 1, MTBDFT
   do ii = 1, MTBDFT
      if (ii == jj) then
         rho(ii,jj,1) = cmplx(rhoa_TBDFT(ii,jj), 0.0D0)
      else
         rho(ii,jj,1) = cmplx(rhoa_TBDFT(ii,jj)/2.0d0, 0.0D0)
      end if
   end do
   end do

   ! Open shell option
   if (thrddim == 2) then
      do jj = 1, MTBDFT
      do ii = 1, MTBDFT
         if (ii == jj) then
            rho(ii,jj,2) = cmplx(rhob_TBDFT(ii,jj), 0.0D0)
         else
            rho(ii,jj,2) = cmplx(rhob_TBDFT(ii,jj)/2.0d0, 0.0D0)
         end if
      end do
      end do
   end if

   do jj = 1, M_in
   do ii = 1, M_in
      rho(ii+MTB,jj+MTB,:) = rho_0(ii,jj,:)
   end do
   end do

end subroutine tbdft_td_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine getXY_TBDFT(M_in,x_in,y_in,xmat,ymat)

   use tbdft_data, only: MTB, MTBDFT, n_biasTB

   implicit none
   integer     , intent(in)  :: M_in
   real(kind=8), intent(in)  :: x_in(M_in,M_in)
   real(kind=8), intent(in)  :: y_in(M_in,M_in)
   real(kind=8), intent(out) :: xmat(MTBDFT,MTBDFT)
   real(kind=8), intent(out) :: ymat(MTBDFT,MTBDFT)
   integer :: ii, jj

   xmat = 0.0d0
   ymat = 0.0d0

   do ii = 1, MTB
      xmat(ii,ii) = 1.0d0
      ymat(ii,ii) = 1.0d0
   end do

   do jj = 1, M_in
   do ii = 1, M_in
      xmat(MTB+ii, MTB+jj) = x_in(ii,jj)
      ymat(MTB+ii, MTB+jj) = y_in(ii,jj)
   end do
   end do
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine construct_rhoTBDFT(M, rho, rho_0 ,rho_TBDFT, niter, open_shell)

   use tbdft_data , only: MTB, MTBDFT

   implicit none
   logical     , intent(in)  :: open_shell
   integer     , intent(in)  :: M
   integer     , intent(in)  :: niter
   real(kind=8), intent(in)  :: rho_0(M,M)
   real(kind=8), intent(in)  :: rho_TBDFT(MTBDFT,MTBDFT)
   real(kind=8), intent(out) :: rho(MTBDFT,MTBDFT)
   integer      :: ii, jj
   real(kind=8) :: ocup

   ocup = 1.0d0
   if (open_shell) ocup = 0.5d0

   if (niter/=1) then

      do ii = 1   , MTBDFT
      do jj = ii+1, MTBDFT
         rho(ii,jj) = rho_TBDFT(ii,jj) / 2
         rho(jj,ii) = rho(ii,jj)
      end do
      end do

   else if (niter == 1) then

      rho = 0.0D0
      do ii = 1, MTB
         rho(ii,ii) = ocup
      end do

      rho(MTB+1:MTB+M,MTB+1:MTB+M) = rho_0

   end if

end subroutine construct_rhoTBDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine build_chimera_TBDFT (M_in,fock_in, fock_TBDFT, natom)

   use tbdft_data, only: MTBDFT, MTB, Iend_TB, end_bTB, alfaTB, betaTB, &
                         gammaTB, Vbias_TB, gammaW, n_biasTB, n_atperbias,n_atTB

   integer     , intent(in)  :: M_in
   integer     , intent(in)  :: natom
   real(kind=8), intent(in)  :: fock_in (M_in, M_in)
   real(kind=8), intent(out) :: fock_TBDFT (MTBDFT, MTBDFT)
   integer :: ii, jj, kk, link

   fock_TBDFT(:,:) = 0.0D0

   do jj = 1, end_bTB*n_atperbias
   do ii = 1, n_biasTB
      link= MTB-n_biasTB+ii
      fock_TBDFT(link, Iend_TB(ii,jj)+MTB) = gammaW(jj) * gammaTB
      fock_TBDFT(Iend_TB(ii,jj)+MTB,link) = gammaW(jj) * gammaTB
   end do
   end do

   do ii = 1,n_biasTB
   do jj = 1,n_atTB-1
      kk=jj+((ii-1)*(n_atTB-1))
      fock_TBDFT(kk,kk) = alfaTB
      if (jj<n_atTB-1) then
         fock_TBDFT(kk,kk+1) = betaTB
         fock_TBDFT(kk+1,kk) = betaTB
      end if
   end do
      link = MTB-n_biasTB+ii
      kk   = ii*(n_atTB-1)
      fock_TBDFT(link,link) = alfaTB
      fock_TBDFT(link,kk) = betaTB
      fock_TBDFT(kk,link) = betaTB
   end do

   fock_TBDFT(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)

end subroutine build_chimera_TBDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine extract_rhoDFT (M_in, rho_in, rho_out)

   use tbdft_data, only: MTBDFT, MTB

   implicit none
   integer     , intent(in)  :: M_in
   real(kind=8), intent(in)  :: rho_in(MTBDFT,MTBDFT)
   real(kind=8), intent(out) :: rho_out(M_in,M_in)
   integer              :: ii, jj

   rho_out=0.0D0

   do jj = 1, M_in
   do ii = 1, M_in
      rho_out(ii,jj) = rho_in(MTB+ii,MTB+jj)
   end do
   end do

end subroutine extract_rhoDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine chimeraTBDFT_evol(M_in,fock_in, fock_TBDFT, natom, istep)

   use tbdft_data, only: MTBDFT, MTB, Iend_TB, end_bTB, alfaTB, betaTB,        &
                         gammaTB, Vbias_TB, start_tdtb, end_tdtb, gammaW,      &
                         n_atTB,n_biasTB, n_atperbias, VbiasTB

   integer     , intent(in)  :: M_in
   integer     , intent(in)  :: natom
   integer     , intent(in)  :: istep
   real(kind=8), intent(in)  :: fock_in(M_in, M_in)
   real(kind=8), intent(out) :: fock_TBDFT(MTBDFT, MTBDFT) !temporal dimensions
   real(kind=8) :: pi = 4.0D0 * atan(1.0D0)
   real(kind=8) :: lambda, t_step, f_t
   integer      :: ii,jj,kk, link

   lambda = 1.0d0 / real(end_tdtb - start_tdtb)

   if (istep < start_tdtb) then
      f_t = 0.0D0
   else if ((istep >= start_tdtb) .and. (istep < end_tdtb)) then
      t_step = real(istep - start_tdtb)
      f_t    = (-cos(pi * lambda * t_step) + 1.0D0) / 2.0D0
   else if (istep >= end_tdtb) then
      f_t = 1.0D0
   end if

   fock_TBDFT(:,:) = 0.0D0

   do jj = 1, end_bTB*n_atperbias
   do ii = 1, n_biasTB
      link= MTB-n_biasTB+ii
      fock_TBDFT(link, Iend_TB(ii,jj)+MTB) = gammaW(jj) * gammaTB
      fock_TBDFT(Iend_TB(ii,jj)+MTB,link) = gammaW(jj) * gammaTB
   end do
   end do

   do ii = 1,n_biasTB
   do jj = 1,n_atTB-1
      kk=jj+((ii-1)*(n_atTB-1))
      fock_TBDFT(kk,kk) = alfaTB+f_t*VbiasTB(ii)
      if (jj<n_atTB-1) then
         fock_TBDFT(kk,kk+1) = betaTB
         fock_TBDFT(kk+1,kk) = betaTB
      end if
   end do
      link = MTB-n_biasTB+ii
      kk   = ii*(n_atTB-1)
      fock_TBDFT(link,link) = alfaTB+f_t*VbiasTB(ii)
      fock_TBDFT(link,kk) = betaTB
      fock_TBDFT(kk,link) = betaTB
   end do

   fock_TBDFT(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)

end subroutine chimeraTBDFT_evol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine TB_current (M_in,rhold,rhonew, overlap, TB_electrode, TB_M)
   use tbdft_data, only:MTBDFT, MTB,n_atTB,n_biasTB

   implicit none
   integer        , intent(in)    :: M_in
   real(kind=8)   , intent(in)    :: overlap(M_in, M_in)
#ifdef TD_SIMPLE
   complex(kind=4), intent(in)    :: rhold(MTBDFT,MTBDFT)
   complex(kind=4), intent(in)    :: rhonew(MTBDFT,MTBDFT)
#else
   complex(kind=8), intent(in)    :: rhold(MTBDFT,MTBDFT)
   complex(kind=8), intent(in)    :: rhonew(MTBDFT,MTBDFT)
#endif
   real(kind=8)   , intent(inout) :: TB_M
   real(kind=8)   , intent(inout) :: TB_electrode(n_biasTB)
   integer      :: ii, jj, kk
   real(kind=8) :: qe
   real(kind=8), allocatable :: delta_rho(:,:)

   allocate(delta_rho(MTBDFT,MTBDFT))

   delta_rho = real(rhonew) - real(rhold)

   do ii=1,n_biasTB
   do jj = 1, n_atTB-1
      kk=jj+((ii-1)*(n_atTB-1))
      TB_electrode(ii) = TB_electrode(ii) + delta_rho(kk,kk)
   end do
      TB_electrode(ii) = TB_electrode(ii) + delta_rho(MTB-n_biasTB+ii,       &
                                                      MTB-n_biasTB+ii)
   end do

   do ii = 1,M_in
   do jj = 1,M_in
      qe = delta_rho(MTB+ii,MTB+jj) * overlap(ii, jj)
      TB_M = qe + TB_M
   enddo
   enddo

end subroutine TB_current

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_scf_output(M_in, open_shell)
   use tbdft_data, only: rhoa_TBDFT, rhob_TBDFT, MTBDFT, MTB, n_biasTB, n_atTB

   implicit none
   logical, intent(in) :: open_shell
   integer, intent(in) :: M_in
   real(kind=8) :: rho_aux(MTBDFT, MTBDFT)
   real(kind=8) :: chargeTB(n_biasTB)
   integer      :: ii, jj,kk


   if (open_shell) then
      rho_aux = rhoa_TBDFT + rhob_TBDFT
   else
      rho_aux = rhoa_TBDFT
   end if

   chargeTB = n_atTB

   do ii=1,n_biasTB
   do jj = 1, n_atTB-1
      kk=jj+((ii-1)*(n_atTB-1))
      chargeTB(ii) = chargeTB(ii) - rho_aux(kk,kk)
   end do
      chargeTB(ii) = chargeTB(ii) - rho_aux(MTB-n_biasTB+ii,MTB-n_biasTB+ii)
   end do

   open(unit=20202, file='mullikenTB')
   do ii=1,n_biasTB
      write(20202,*) "Mulliken TB  electro", ii, chargeTB(ii)
   end do
   close(20202)
end subroutine tbdft_scf_output

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine tbdft_td_output(M_in, thrddim, rho_aux, overlap, istep, Iz, natom, &
                           Nuc, open_shell)
   use tbdft_data, only: rhold_AOTB, rhonew_AOTB, MTB, MTBDFT,n_atTB, n_biasTB
   implicit none

   logical        , intent(in) :: open_shell
   integer        , intent(in) :: M_in, istep, thrddim
   integer        , intent(in) :: natom
   integer        , intent(in) :: Nuc(M_in)
   integer        , intent(in) :: Iz(natom)
   real(kind=8)   , intent(in)  :: overlap(M_in,M_in)
#ifdef TD_SIMPLE
   complex(kind=4), intent(in)  :: rho_aux(MTBDFT, MTBDFT,thrddim)
#else
   complex(kind=8), intent(in) :: rho_aux(MTBDFT, MTBDFT,thrddim)
#endif
   integer      :: ii,jj,kk
   real(kind=8) :: I_TB_elec(n_biasTB)
   real(kind=8) :: I_TB_M
   real(kind=8) :: chargeTB(n_biasTB)
   real(kind=8) :: chargeM_TB
   real(kind=8) :: orb_charge, tot_orb_charge
   real(kind=8) :: qe(natom)
   real(kind=8) :: rhoscratch(M_in,M_in)

   I_TB_M    = 0.0d0
   I_TB_elec = 0.0d0

   if (istep == 1) then
      open(unit=10101,file='currentTB')
      open(unit=20202,file='mullikenTB')

   else
      call TB_current(M_in,rhold_AOTB(:,:,1),rhonew_AOTB(:,:,1), overlap, &
                         I_TB_elec, I_TB_M)
      if (open_shell) then
         call TB_current(M_in,rhold_AOTB(:,:,2),rhonew_AOTB(:,:,2), overlap, &
                         I_TB_elec, I_TB_M)
      end if

      do ii=1, n_biasTB
         write(10101,*) "Current TB electrode", ii, I_TB_elec(ii)
      end do
      write(10101,*) "Current DFT part M", I_TB_M

      chargeTB = n_atTB
      do ii=1,n_biasTB
      do jj = 1, n_atTB-1
         kk=jj+((ii-1)*(n_atTB-1))
         chargeTB(ii) = chargeTB(ii) - rho_aux(kk,kk,1)
      end do
         chargeTB(ii) = chargeTB(ii) - rho_aux(MTB-n_biasTB+ii,              &
                                               MTB-n_biasTB+ii,1)
      end do

      if (open_shell) then
         do ii=1,n_biasTB
         do jj = 1, n_atTB-1
            kk=jj+((ii-1)*(n_atTB-1))
            chargeTB(ii) = chargeTB(ii) - rho_aux(kk,kk,2)
         end do
            chargeTB(ii) = chargeTB(ii) - rho_aux(MTB-n_biasTB+ii,           &
                                                  MTB-n_biasTB+ii,2)
         end do
      end if

      chargeM_TB = 0.0D0
      do ii = 1,natom
         qe(ii) = Iz(ii)
      enddo

      rhoscratch = real(rho_aux(MTB+1:MTB+M_in,MTB+1:MTB+M_in,1))

      call mulliken_calc(natom, M_in, rhoscratch, overlap, Nuc, qe)

      do ii = 1, natom
         chargeM_TB = chargeM_TB + qe(ii)
      enddo

      if (open_shell) then
         rhoscratch = real(rho_aux(MTB+1:MTB+M_in,MTB+1:MTB+M_in,2))

         call mulliken_calc(natom, M_in, rhoscratch, overlap, Nuc, qe)

         do ii = 1,natom
            chargeM_TB = chargeM_TB + qe(ii)
         enddo
      end if

      do ii=1,n_biasTB
         write(20202,*) "Mulliken TB electrode",ii,chargeTB(ii)
      end do
      write(20202,*) "Mulliken DFT  part M", chargeM_TB

   endif
end subroutine tbdft_td_output

end module tbdft_subs