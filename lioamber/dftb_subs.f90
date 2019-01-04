!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_subs
   implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_init(M, OPEN)

   use dftb_data, only: MTB, MDFTB, end_bTB, IendA_TB,IendB_TB , rho_aDFTB, rho_bDFTB,    &
                        gammaW, dftb_transport, rhofirst_TB, coup_atoms_A,     &
                        coup_atoms_B, at_coupled

   implicit none

   integer, intent(in) :: M
   logical, intent(in) :: OPEN
   real*8, allocatable :: rhoTB_real(:,:,:)
   integer   ::  ii,jj

   MDFTB=2*MTB+M
   allocate(rho_aDFTB(MDFTB,MDFTB))
   if (OPEN) allocate (rho_bDFTB(MDFTB,MDFTB))

   if (dftb_transport/=0) then
      if (OPEN) then
         allocate(rhoTB_real(MDFTB,MDFTB,2), rhofirst_TB(MDFTB,MDFTB,2))
      else
         allocate(rhoTB_real(MDFTB,MDFTB,1),rhofirst_TB(MDFTB,MDFTB,1))
      end if
   end if

   open(unit=1001,file='gamma.in')
   read(1001,*) at_coupled
   allocate(coup_atoms_A(at_coupled),coup_atoms_B(at_coupled))
   allocate(IendA_TB(at_coupled,end_bTB),IendB_TB(at_coupled,end_bTB))
   allocate(gammaW(at_coupled,end_bTB))

   read(1001,*) coup_atoms_A
   read(1001,*) coup_atoms_B

   do jj=1, at_coupled
   do ii=1, end_bTB
      read(1001,*) gammaW(jj, ii)
   enddo
   enddo

   close(1001)

   if (dftb_transport==2) then
      open(unit=1001,file='rhofirstTB')


      if(OPEN) then
         do jj=1, MDFTB
         do ii=1, MDFTB
            read(1001,*) rhoTB_real(ii,jj,1),rhoTB_real(ii,jj,2)
         enddo
         enddo
      else
         do jj=1, MDFTB
         do ii=1, MDFTB
            read(1001,*) rhoTB_real(ii,jj,1)
         enddo
         enddo
      end if

      close(1001)

      do jj=1,MDFTB
      do ii=1,MDFTB
         if (ii==jj) then
            rhofirst_TB(ii,jj,1) = dcmplx(rhoTB_real(ii,jj,1),0.0d0)
         else
            rhofirst_TB(ii,jj,1) = dcmplx(0.5d0*rhoTB_real(ii,jj,1),0.0d0)
         end if
      end do
      end do

      if(OPEN) then
         do jj=1,MDFTB
         do ii=1,MDFTB
            if (ii==jj) then
               rhofirst_TB(ii,jj,2) = dcmplx(rhoTB_real(ii,jj,2),0.0d0)
            else
               rhofirst_TB(ii,jj,2) = dcmplx(0.5d0*rhoTB_real(ii,jj,2),0.0d0)
            end if
         end do
         end do
      end if
      write(*,*) "RHOFIRST_TB READ"
   end if


end subroutine dftb_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dftb_td_init (M,rho, rho_0, overlap, RMM5, dim3)

   use dftb_data, only: MTB, MDFTB, rho_aDFTB, rho_bDFTB, rhold_AOTB,          &
                        rhonew_AOTB
   implicit none
   integer, intent(in)       :: M
!carlos: dim3 is to declare the 3th dimension of matrices
   integer, intent(in)       :: dim3
   real*8, allocatable, intent(inout) :: overlap(:,:)
   real*8 , intent(in)  :: RMM5(M*(M+1)/2)
#ifdef TD_SIMPLE
   complex*8, intent(in)     :: rho_0(M,M,dim3)
   complex*8, intent(out)  :: rho(MDFTB,MDFTB,dim3)
#else
   complex*16, intent(in)    :: rho_0(M,M,dim3)
   complex*16, intent(out) :: rho(MDFTB,MDFTB,dim3)
#endif
   integer  :: ii,jj

   if (allocated(overlap)) deallocate(overlap)
   allocate(overlap(M,M), rhold_AOTB(MDFTB,MDFTB,dim3),                       &
            rhonew_AOTB(MDFTB,MDFTB,dim3))

   rhold_AOTB=0.0d0
   rhonew_AOTB=0.0d0

   call spunpack('L', M, RMM5, overlap)

      do jj=1, MDFTB
      do ii=1, MDFTB
         if (ii==jj) then
            rho(ii,jj,1)=cmplx(rho_aDFTB(ii,jj), 0.0D0)
         else
            rho(ii,jj,1)=cmplx(rho_aDFTB(ii,jj)/2.0d0, 0.0D0)
         end if
      end do
      end do

!carlos:Open shell option
   if(dim3==2) then
      do jj=1, MDFTB
      do ii=1, MDFTB
         if (ii==jj) then
            rho(ii,jj,2)=cmplx(rho_bDFTB(ii,jj), 0.0D0)
         else
            rho(ii,jj,2)=cmplx(rho_bDFTB(ii,jj)/2.0d0, 0.0D0)
         end if
      end do
      end do
   end if

   do jj= 1, M
   do ii= 1, M
      rho(ii+MTB,jj+MTB,:)=rho_0(ii,jj,:)
   end do
   end do

end subroutine dftb_td_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i%%%%%%%%%%%%%%%%%%%%%%%%

subroutine getXY_DFTB(M_in,x_in,y_in,xmat,ymat)

   use dftb_data, only: MTB, MDFTB

   implicit none
   integer, intent(in)  :: M_in
   real*8 , intent(in)  :: x_in(M_in,M_in)
   real*8 , intent(in)  :: y_in(M_in,M_in)
   real*8 , intent(out) :: xmat(MDFTB,MDFTB)
   real*8 , intent(out) :: ymat(MDFTB,MDFTB)
   integer              :: ii,jj


   xmat=0.0d0
   ymat=0.0d0

   do ii=1, MTB
      xmat(ii,ii)=1.0d0
      xmat(MTB+M_in+ii,MTB+M_in+ii)=1.0d0

      ymat(ii,ii)=1.0d0
      ymat(MTB+M_in+ii,MTB+M_in+ii)=1.0d0
   end do

   do jj=1, M_in
   do ii=1, M_in
      xmat(MTB+ii, MTB+jj)=x_in(ii,jj)
      ymat(MTB+ii, MTB+jj)=y_in(ii,jj)
   end do
   end do
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine read_rhoDFTB(M, MM, RMM1, rhoalpha, rhobeta, OPEN)
   use dftb_data, only: MTB, MDFTB, rho_aDFTB, rho_bDFTB

   implicit none
   integer, intent(in)    :: M, MM
   real*8, intent(inout)  :: RMM1(MM)
   real*8, intent(inout)  :: rhoalpha(MM), rhobeta(MM)
   logical, intent(in)    :: OPEN
   integer                :: ii, jj, kk


  open(unit=1070,file='rhoTB.in')

   if (.not.OPEN) then
      do ii=1,MDFTB
      do jj=1,MDFTB
         read(1070,*) rho_aDFTB(ii,jj)
      end do
      end do

      do jj=1,M
      do kk=jj,M
               RMM1(kk+(2*M-jj)*(jj-1)/2)=rho_aDFTB(jj+MTB,kk+MTB)
      enddo
      enddo

   else if(OPEN) then

      do ii=1,MDFTB
      do jj=1,MDFTB
         read(1070,*) rho_aDFTB(ii,jj), rho_bDFTB(ii,jj)
      end do
      end do

      do jj=1,M
      do kk=jj,M
               rhoalpha(kk+(2*M-jj)*(jj-1)/2)=rho_aDFTB(jj+MTB,kk+MTB)
               rhobeta(kk+(2*M-jj)*(jj-1)/2)=rho_aDFTB(jj+MTB,kk+MTB)
      enddo
      enddo

      RMM1 = rhoalpha+rhobeta

      write(*,*) 'RHOTB readed'

   end if

    close(1070)
end subroutine read_rhoDFTB
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine construct_rhoDFTB(M, rho, rho_0 ,rho_DFTB, TBload, niter,OPEN)
!carlos: corrected the exces of electrons in open shell
   use dftb_data , only: MTB, MDFTB

   implicit none
   integer, intent(in) :: M
   integer, intent(in) :: niter
   real*8, intent(in)  :: rho_0(M,M)
   real*8, intent(in)  :: rho_DFTB(MDFTB,MDFTB)
   logical, intent(in) :: TBload
   logical, intent(in) :: OPEN
   real*8, intent(out) :: rho(MDFTB,MDFTB)
   real*8  :: pop
   integer :: ii, jj

   pop = 1.0d0
   if (OPEN) pop=0.5d0

   if ((TBload.and.niter==1).or.(niter/=1)) then

      do ii=1,MDFTB
      do jj=ii+1,MDFTB
         rho(ii,jj)=rho_DFTB(ii,jj)/2
         rho(jj,ii)=rho(ii,jj)
      end do
      end do

   else if(.not.TBload.and.(niter==1)) then

      rho=0.0D0
      do ii=1, MTB
         rho(ii,ii)=pop
         rho(MTB+M+ii,MTB+M+ii)=pop
      end do
         rho(MTB+1:MTB+M, MTB+1:MTB+M)=rho_0
   end if

end subroutine construct_rhoDFTB
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine find_TB_neighbors(M_in,Nuc,natom)

   use dftb_data, only:IendA_TB, IendB_TB, at_coupled, coup_atoms_A, coup_atoms_B

   implicit none

   integer              ::  ii, jj, kk
   integer, intent(in)  ::  M_in, natom
   integer, intent(in)  ::  Nuc(M_in)


   do jj=1, at_coupled
      kk=0
   do ii=1, M_in
      if (coup_atoms_A(jj)==Nuc(ii)) then
         kk=kk+1
         IendA_TB(jj,kk)=ii
      end if
   end do
   end do

   do jj=1, at_coupled
      kk=0
   do ii=1, M_in
      if (coup_atoms_B(jj)==Nuc(ii)) then
         kk=kk+1
         IendB_TB(jj,kk)=ii
      end if
   end do
   end do

end subroutine find_TB_neighbors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine build_chimera_DFTB (M_in,fock_in, fock_DFTB, natom, nshell, ncont)

   use dftb_data, only:MDFTB, MTB, IendA_TB, IendB_TB, end_bTB, alfaTB, betaTB, gammaTB,  &
                       Vbias_TB, gammaW, Vbias_TB, dftb_transport, at_coupled

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   integer, intent(in)  :: nshell (0:3)
   integer, intent(in)  :: ncont(M_in)
   real*8, intent(in)   :: fock_in (M_in, M_in)
   real*8, intent(out)  :: fock_DFTB (MDFTB, MDFTB)
   integer              :: ii, jj, kk, ns, np, l1, l2, link
   real*8               :: V_aux
   real*8               :: f_err(MTB)
   real*8               :: pi=4.0d0*datan(1.0d0)
   real*8               :: d_err

   l1=0
   l2=0
   jj=0
   kk=0
   ns=nshell(0)
   np=nshell(1)
   V_aux= 0.0d0
   f_err= 0.0d0
   d_err= 0.0d0

   ! d_err=2.0d0*pi/dble(MTB-1)

!charly:this most be change, it's messy:
   V_aux=Vbias_TB
   if (dftb_transport==0) V_aux=0.0d0

   ! do ii=1, MTB
   !    f_err(ii) = (1.0d0-erf(-pi+dble(ii-1)*d_err))
   ! end do


   fock_DFTB(:,:) = 0.0D0

   do jj=1, at_coupled
   do ii=1, end_bTB

      fock_DFTB(IendA_TB(jj,ii)+MTB,MTB)=gammaW(jj,ii)*gammaTB
      fock_DFTB(MTB,IendA_TB(jj,ii)+MTB)=gammaW(jj,ii)*gammaTB
      fock_DFTB(IendB_TB(jj,ii)+MTB,MTB+M_in+1)=gammaW(jj,ii)*gammaTB
      fock_DFTB(MTB+M_in+1,IendB_TB(jj,ii)+MTB)=gammaW(jj,ii)*gammaTB

   end do
   end do

   do ii=1,MTB
      fock_DFTB(ii,ii) = alfaTB-V_aux/2.0d0
      fock_DFTB(MTB+M_in+ii, MTB+M_in+ii)= alfaTB+V_aux/2.0d0

      if (ii<MTB) then !-end_bTB) then

         fock_DFTB(ii,ii+1) = betaTB
         fock_DFTB(ii+1,ii) = betaTB
         fock_DFTB(2*MTB+M_in-ii, 2*MTB+M_in-ii+1)= betaTB
         fock_DFTB(2*MTB+M_in-ii+1, 2*MTB+M_in-ii)= betaTB

      end if
   end do


   fock_DFTB(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)


end subroutine build_chimera_DFTB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine extract_rhoDFT (M_in, rho_in, rho_out)

   use dftb_data, only: MDFTB, MTB

   implicit none
   integer, intent(in)  :: M_in
   real*8, intent(in)   :: rho_in(MDFTB,MDFTB)
   real*8, intent(out)  :: rho_out(M_in,M_in)
   integer              :: ii, jj

   rho_out=0.0D0

   do jj=1, M_in
   do ii=1, M_in
      rho_out(ii,jj)=rho_in(MTB+ii,MTB+jj)
   end do
   end do

end subroutine extract_rhoDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine chimeraDFTB_evol(M_in,fock_in, fock_DFTB, natom, nshell,ncont, istep)

   use dftb_data, only:MDFTB, MTB, IendA_TB, IendB_TB, end_bTB, alfaTB, betaTB, gammaTB,   &
                       Vbias_TB, start_tdtb, end_tdtb, gammaW, dftb_transport, at_coupled

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   integer, intent(in)  :: ncont(M_in)
   integer, intent(in)  :: istep
   integer, intent(in)  :: nshell (0:3)
   real*8, intent(in)   :: fock_in (M_in, M_in)
   real*8, intent(out)  :: fock_DFTB (MDFTB, MDFTB) !temporal dimensions
   real*8               :: lambda, t_step, f_t
   integer              :: ii, jj, kk, ns, np, l1, l2, link
   real*8               :: f_err(MTB)
   real*8               :: pi=4.0d0*datan(1.0d0)
   real*8               :: d_err


   l1=0
   l2=0
   jj=0
   kk=0
   ns=nshell(0)
   np=nshell(1)
   f_err= 0.0d0
   d_err= 0.0d0

  ! d_err=2.0d0*pi/dble(MTB-1)

   ! do ii=1, MTB
   !    f_err(ii) = (1.0d0-erf(-pi+dble(ii-1)*d_err))
   ! end do

   print*, "istep, start, end", istep, start_tdtb, end_tdtb

   lambda=1.0d0/real(end_tdtb-start_tdtb)

!charly: this is also messy
   if(dftb_transport==0) then
      if (istep < start_tdtb) then
         f_t=0.0D0
      else if(istep >= start_tdtb .and. istep < end_tdtb) then
         t_step=real(istep-start_tdtb)
         f_t=(-Cos(pi*lambda*t_step)+1.0D0)/2.0D0
      else if(istep >= end_tdtb) then
         f_t=1.0D0
      end if
   else
      f_t=1.0D0
   end if

   ! if (istep < start_tdtb) then
   !    f_t=1.0D0
   ! else if(istep >= start_tdtb .and. istep < end_tdtb) then
   !    t_step=real(istep-start_tdtb)
   !    f_t=(Cos(pi*lambda*t_step)+1.0D0)/2.0D0
   ! else if(istep >= end_tdtb) then
   !    f_t=0.0D0
   ! end if

   print*, "factor V", f_t

   fock_DFTB(:,:) = 0.0D0

   do jj=1, at_coupled
   do ii=1, end_bTB

      fock_DFTB(IendA_TB(jj,ii)+MTB,MTB)=gammaW(jj,ii)*gammaTB
      fock_DFTB(MTB,IendA_TB(jj,ii)+MTB)=gammaW(jj,ii)*gammaTB
      fock_DFTB(IendB_TB(jj,ii)+MTB,MTB+M_in+1)=gammaW(jj,ii)*gammaTB
      fock_DFTB(MTB+M_in+1,IendB_TB(jj,ii)+MTB)=gammaW(jj,ii)*gammaTB

   end do
   end do


   do ii=1,MTB
      fock_DFTB(ii,ii) = alfaTB-(Vbias_TB/2.0d0)*f_t!*f_err(ii)
      fock_DFTB(MTB+M_in+ii, MTB+M_in+ii)= alfaTB+(Vbias_TB/2.0d0)*f_t!*f_err(MTB+1-ii)

      if (ii<MTB) then ! -end_bTB) then
         fock_DFTB(ii,ii+1) = betaTB
         fock_DFTB(ii+1,ii) = betaTB
         fock_DFTB(2*MTB+M_in-ii, 2*MTB+M_in-ii+1)= betaTB
         fock_DFTB(2*MTB+M_in-ii+1, 2*MTB+M_in-ii)= betaTB

      end if
   end do


   fock_DFTB(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)


end subroutine chimeraDFTB_evol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TB_current (M_in,rhold,rhonew, overlap, TB_A, TB_B, TB_M)
   use dftb_data, only:MDFTB, MTB

   implicit none
   integer, intent(in)  :: M_in
   real*8, intent(in)   :: overlap(M_in, M_in)
#ifdef TD_SIMPLE
   complex*8, intent(in)   :: rhold(MDFTB,MDFTB)
   complex*8, intent(in)   :: rhonew(MDFTB,MDFTB)
#else
   complex*16, intent(in)   :: rhold(MDFTB,MDFTB)
   complex*16, intent(in)   :: rhonew(MDFTB,MDFTB)
#endif
   real*8, intent(out)  :: TB_A, TB_B, TB_M
   integer              :: ii, jj
   real*8, allocatable  :: delta_rho(:,:)
   real*8               :: qe

   allocate(delta_rho(MDFTB,MDFTB))

   delta_rho=real(rhonew)-real(rhold)


   TB_A=0.0D0
   TB_B=0.0D0
   TB_M=0.0D0

   do ii=1, MTB

      TB_A=delta_rho(ii,ii) + TB_A
      TB_B=delta_rho(MTB+M_in+ii, MTB+M_in+ii) + TB_B

   end do

   do ii=1,M_in
   do jj=1,M_in

      qe = delta_rho(MTB+ii,MTB+jj) * overlap(ii, jj)
      TB_M = qe + TB_M

   enddo
   enddo

end subroutine TB_current

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine write_rhoDFTB(M_in, OPEN)
   use dftb_data , only: rho_aDFTB, rho_bDFTB, dftb_transport, rhofirst_TB

   implicit none
   integer, intent(in) :: M_in
   logical, intent(in) :: OPEN
   integer :: ii, jj

   open(unit=1070,file='rhoTB.out')
   if(dftb_transport==1) open(unit=1001,file='rhofirstTB')

   if (OPEN) then
      do jj=1,M_in
      do ii=1,M_in
         write(1070,*) rho_aDFTB(ii,jj), rho_bDFTB(ii,jj)
         if(dftb_transport==1) write(1001,*) rho_aDFTB(ii,jj), rho_bDFTB(ii,jj)
      end do
      end do
   else
      do jj=1,M_in
      do ii=1,M_in
         write(1070,*) rho_aDFTB(ii,jj)
         if(dftb_transport==1) write(1001,*) rho_aDFTB(ii,jj)
      end do
      end do
   end if

   close(1070)

   write(*,*) 'RHOTB wrtted'

   if (dftb_transport==1) then
      rhofirst_TB(:,:,1) = rho_aDFTB(:,:)
      if (OPEN) rhofirst_TB(:,:,2) = rho_bDFTB(:,:)
      close(1001)
      write(*,*) 'RHOFIRST_TB writed'
   end if

end subroutine write_rhoDFTB


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_output(M, dim3, rho_aux, overlap, istep, Iz, natom, Nuc, OPEN)
   use dftb_data, only: rhold_AOTB, rhonew_AOTB, MTB, MDFTB, dftb_transport
   implicit none

   logical, intent(in) :: OPEN
   integer, intent(in) :: M, istep, dim3
   integer, intent(in) :: natom
   integer, intent(in) :: Nuc(M)
   integer, intent(in) :: Iz(natom)
   real*8, intent(in)  :: overlap(M,M)
#ifdef TD_SIMPLE
   complex*8, intent(in)  :: rho_aux(MDFTB, MDFTB,dim3)
#else
   complex*16, intent(in) :: rho_aux(MDFTB, MDFTB,dim3)
#endif
   integer  :: n, ii, jj
   real*8   :: I_TB_A(dim3), I_TB_B(dim3), I_TB_M(dim3)
   real*8   :: chargeA_TB, chargeB_TB, chargeM_TB
   real*8   :: orb_charge, tot_orb_charge
   real*8   :: qe(natom)
   real*8   :: rhoscratch(M,M)

   if (istep==1) then
!charly: the current is temp off
      if(dftb_transport==0)open(unit=10101,file='currentTB')
      open(unit=20202,file='mullikenTB')

   else
      if(dftb_transport==0) then
      if (OPEN) then
         call TB_current(M,rhold_AOTB(:,:,1),rhonew_AOTB(:,:,1), overlap,      &
                         I_TB_A(1), I_TB_B(1), I_TB_M(1))
         call TB_current(M,rhold_AOTB(:,:,2),rhonew_AOTB(:,:,2), overlap,      &
                         I_TB_A(2), I_TB_B(2), I_TB_M(2))

         write(10101,*) "A", I_TB_A(1) + I_TB_A(2)
         write(10101,*) "B", I_TB_B(1) + I_TB_B(2)
         write(10101,*) "M", I_TB_M(1) + I_TB_M(2)
      else
         call TB_current(M,rhold_AOTB(:,:,1),rhonew_AOTB(:,:,1), overlap,      &
                         I_TB_A(1), I_TB_B(1), I_TB_M(1))

         write(10101,*) "A", I_TB_A(1)
         write(10101,*) "B", I_TB_B(1)
         write(10101,*) "M", I_TB_M(1)
      end if
      end if


      chargeA_TB=MTB
      chargeB_TB=MTB
      do ii=1, MTB
         chargeA_TB=chargeA_TB-real(rho_aux(ii,ii,1))
         chargeB_TB=chargeB_TB-real(rho_aux(MTB+M+ii,MTB+M+ii,1))
      end do

      if (OPEN) then
         do ii=1, MTB
            chargeA_TB=chargeA_TB-real(rho_aux(ii,ii,2))
            chargeB_TB=chargeB_TB-real(rho_aux(MTB+M+ii,MTB+M+ii,2))
         end do
      end if

      chargeM_TB=0.0D0
      do n=1,natom
         qe(n)=Iz(n)
      enddo

      rhoscratch=real(rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,1))

      call mulliken_calc(natom, M, rhoscratch, overlap, Nuc, qe)

      if(OPEN) then
         rhoscratch=real(rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,2))

         call mulliken_calc(natom, M, rhoscratch, overlap, Nuc, qe)

      end if

      do n=1,natom
         chargeM_TB= chargeM_TB + qe(n)
      enddo

      write(20202,*) "A", chargeA_TB
      write(20202,*) "B", chargeB_TB
      write(20202,*) "M", chargeM_TB

   endif

end subroutine dftb_output
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine transport_TB(M, natom, dim3, overlap, rho_aux ,devPtrY,Nuc,istep, OPEN)
   use dftb_data, only: MTB, MDFTB, rhofirst_TB, driving_rateTB
   use cublasmath    , only: basechange_cublas

   implicit none
   logical, intent(in)      :: OPEN
   integer, intent(in)      :: M, natom, dim3, istep
   integer, intent(in)      :: Nuc(M)
   integer*8,intent(in)     :: devPtrY
   real*8 , intent(in)      :: overlap(M,M)
   complex*8, intent(inout) :: rho_aux(MDFTB,MDFTB,dim3)
   real*8     :: rho_real(MDFTB,MDFTB,dim3)
   real*8     :: scratchgamma
   real*8     :: qe(natom)
   real*8     :: currentA, currentB, currentMol
   complex*8  :: rhoscratch(MDFTB,MDFTB,dim3)
   integer    :: ii, jj

   rhoscratch   = 0.0d0
   rho_real     = 0.0d0
   scratchgamma = 0.0d0
   qe           = 0.0d0
   currentA     = 0.0d0
   currentB     = 0.0d0
   currentMol   = 0.0d0

   if (istep>200.and.istep<=1200) then
      scratchgamma= driving_rateTB * dexp(-0.0001D0 *                              &
                    (dble(istep-1200)**2))
   else if(istep<=200) then
      scratchgamma = 0.0d0
   else if(istep>1200) then
      scratchgamma = driving_rateTB
   end if
!carlos: The driving term is calculated

   do jj=1,MTB
   do ii=1,MTB
      rhoscratch(ii,jj,1) = rho_aux(ii,jj,1) - rhofirst_TB(ii,jj,1)
      rhoscratch(ii+M+MTB,jj+M+MTB,1) = rho_aux(ii+M+MTB,jj+M+MTB,1) -          &
                                        rhofirst_TB(ii+M+MTB,jj+M+MTB,1)
      if(OPEN) then
         rhoscratch(ii,jj,2) = rho_aux(ii,jj,2) - rhofirst_TB(ii,jj,2)
         rhoscratch(ii+M+MTB,jj+M+MTB,2) = rho_aux(ii+M+MTB,jj+M+MTB,2) -       &
                                        rhofirst_TB(ii+M+MTB,jj+M+MTB,2)
      end if
   end do
   end do

   do jj=MTB+1,MTB+M
   do ii=1,MTB
      rhoscratch(ii,jj,1) = 0.5d0*(rho_aux(ii,jj,1) - rhofirst_TB(ii,jj,1))
      if(OPEN) rhoscratch(ii,jj,2) = 0.5d0*(rho_aux(ii,jj,2) -                  &
                                            rhofirst_TB(ii,jj,2))
      rhoscratch(jj,ii,:) = rhoscratch(ii,jj,:)
   end do
   end do

   do jj=MTB+M+1,MDFTB
   do ii=1,MTB
      rhoscratch(ii,jj,1) = rho_aux(ii,jj,1) - rhofirst_TB(ii,jj,1)
      if(OPEN) rhoscratch(ii,jj,2) = rho_aux(ii,jj,2) - rhofirst_TB(ii,jj,2)
      rhoscratch(jj,ii,:) = rhoscratch(ii,jj,:)
   end do
   end do

   do jj=MTB+M+1,MDFTB
   do ii=MTB+1,MTB+M
      rhoscratch(ii,jj,1) = 0.5d0*(rho_aux(ii,jj,1)-rhofirst_TB(ii,jj,1))
      if(OPEN) rhoscratch(ii,jj,2) = 0.5d0*(rho_aux(ii,jj,2) -                  &
                                           rhofirst_TB(ii,jj,2))
      rhoscratch(jj,ii,:) = rhoscratch(ii,jj,:)
   end do
   end do

   rhoscratch = scratchgamma*rhoscratch

   rho_real = real(rhoscratch)

   call mulliken_calc(natom, M, rho_real(MTB+1:MTB+M,MTB+1:MTB+M,1), overlap,   &
                      Nuc, qe)
   if (OPEN) call mulliken_calc(natom, M, rho_real(MTB+1:MTB+M,MTB+1:MTB+M,2),  &
                                overlap,Nuc, qe)

   do ii=1,natom
      currentMol=currentMol+qe(ii)
   end do

   do ii=1,MTB
      currentA  = currentA + rho_real(ii,ii,1)
      if (OPEN) currentA  = currentA + rho_real(ii,ii,2)
      currentB  = currentB + rho_real(ii+MTB+M,ii+MTB+M,1)
      if (OPEN) currentB  = currentB + rho_real(ii+MTB+M,ii+MTB+M,2)
   end do

   open(unit=30303,file='current_DLVN_TB')

   write(30303,*) "A", currentA
   write(30303,*) "B", currentB
   write(30303,*) "M", currentMol

   rho_aux(:,:,1) = basechange_cublas(MDFTB, rhoscratch(:,:,1), devPtrY, 'dir')
   if (OPEN) rho_aux(:,:,2) = basechange_cublas(MDFTB, rhoscratch(:,:,2),        &
                                                   devPtrY, 'dir')

end subroutine transport_TB
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module dftb_subs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
