module armonic_subs
   implicit none

contains

!###############################################################################
subroutine grad_calc(E, rho_aop, fock_aop,rho_bop, fock_bop )

   use garcha_mod      , only : r, natom, rqm
   use basis_data      , only : M
   use armonic_data    , only : armonic_calc, hess_norder, delta_h, mass_w,    &
                                armonic_freq, armonic_vec
   use typedef_operator, only : operator

   implicit none
   type(operator), intent(inout)           :: rho_aop, fock_aop
   type(operator), intent(inout), optional :: rho_bop, fock_bop
   real(kind=8)  , intent(inout)           :: E

   real(kind=8)              :: dhau
   real(kind=8)              :: r_init(natom,3)
   real(kind=8)              :: r_new(natom,3)
   real(kind=8)              :: grad0(3*natom)
   real(kind=8), allocatable :: grad(:,:,:)
   real(kind=8), allocatable :: Df_Dr(:,:,:)
   real(kind=8), allocatable :: dxyzqm(:,:)
   real(kind=8), allocatable :: Dfock_a(:,:,:)
   real(kind=8), allocatable :: hess_mat(:,:)
   integer                   :: ii, jj, kk, rr, ss
   ! real(kind=8), allocatable :: DEner(:,:) !TEMPORAL ARRAY FOR TESTS
   ! real(kind=8), allocatable :: f_posta(:,:) !TEMPORAL ARRAY FOR TESTS

   if (armonic_calc == 0) return

   call init_grad_calc(natom)

   dhau = delta_h / 0.529177D0

!Storin input coordinates
   r_init = r * 0.529177D0 ! r is in a.u. and r_init is in Amstrongs

   allocate(dxyzqm(3, natom), hess_mat(3*natom, 3*natom))

   if (armonic_calc == 2) allocate(Dfock_a(M,M,3*natom))

   if (hess_norder==1) then
      allocate (grad(3*natom,-1:1,3*natom))
      if (armonic_calc == 2) allocate(Df_Dr(M,M,-1:1))
   else
      allocate (grad(3*natom,-2:2,3*natom))
      if (armonic_calc == 2) allocate(Df_Dr(M,M,-2:2))
   end if
   grad=0.0D0
   Df_Dr   = 0.0d0
   Dfock_a = 0.0d0

! Calculating energy and gradient at the initial geometry

   call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
   call dft_get_qm_forces(dxyzqm)

   do ii =1, natom
   do jj=1, 3
      grad0(3*(ii-1)+jj) = dxyzqm(jj,ii)
   end do
   end do

!  Calculating gradients in displaced positions along the
!  QM coordinates.

   do ii = 1, natom
   do jj=1, 3
      kk=3*(ii-1)+jj
!     Forward displacement
!     --------------------
      r_new = r_init
      r_new(ii,jj) = r_init(ii,jj) + delta_h * mass_w(ii)
      r = r_new / 0.529177D0
!carlos: this equality is not necessary with the correct dimensions, and is
!        repeated in all the parts of the loop
      rqm = r
      call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
      call dft_get_qm_forces(dxyzqm)

      do rr=1,natom
      do ss=1,3
         grad(kk, 1, 3*(rr-1)+ss) = dxyzqm(ss,rr)
      end do
      end do

      if (armonic_calc == 2) call fock_aop%Gets_data_AO(Df_Dr(:,:,1))

!     Backward displacement
!     --------------------
      r_new = r_init
      r_new(ii,jj) = r_init(ii,jj) - delta_h * mass_w(ii)
      r = r_new / 0.529177D0
      rqm = r
      call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
      call dft_get_qm_forces(dxyzqm)

      do rr=1,natom
      do ss=1,3
         grad(kk, 1, 3*(rr-1)+ss) = dxyzqm(ss,rr)
      end do
      end do

      if (armonic_calc == 2) call fock_aop%Gets_data_AO(Df_Dr(:,:,-1))

      if(hess_norder == 2) then
!        Forward 2x displacement
!        -----------------------
         r_new = r_init
         r_new(ii,jj) = r_init(ii,jj) + 2.0d0 * delta_h * mass_w(ii)
         r = r_new / 0.529177D0
         rqm = r
         call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
         call dft_get_qm_forces(dxyzqm)

         do rr=1,natom
         do ss=1,3
            grad(kk, 1, 3*(rr-1)+ss) = dxyzqm(ss,rr)
         end do
         end do

         if (armonic_calc == 2) call fock_aop%Gets_data_AO(Df_Dr(:,:,2))

!        Backward 2x displacement
!        ------------------------
         r_new = r_init
         r_new(ii,jj) = r_init(ii,jj) - 2.0d0 * delta_h * mass_w(ii)
         r = r_new / 0.529177D0
         rqm = r
         call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
         call dft_get_qm_forces(dxyzqm)

         do rr=1,natom
         do ss=1,3
            grad(kk, 1, 3*(rr-1)+ss) = dxyzqm(ss,rr)
         end do
         end do
         if (armonic_calc == 2) call fock_aop%Gets_data_AO(Df_Dr(:,:,-2))

      endif

      if (armonic_calc == 2) then
         if (hess_norder == 2) then
            Dfock_a(:,:,kk) = (8.0d0*Df_Dr(M,M,1) - 8.0d0*Df_Dr(M,M,-1) +      &
                               Df_Dr(M,M,-2) - Df_Dr(M,M,2))                   &
                              /(12.0d0 * dhau * mass_w(ii))
         else
            Dfock_a(:,:,kk) = (Df_Dr(M,M,1) - Df_Dr(M,M,-1))*0.5d0             &
                              /(dhau * mass_w(ii))
         end if
      end if

   end do
   end do

!Construction and diagonalization of the hessian matrix

call build_hessian(natom, dhau, grad, hess_mat, mass_w)
call diagonalize_hessian(natom, hess_mat, armonic_freq, armonic_vec)

!charly: checking
   ! write(777,*) "te cabieron las fuerzas?", delta_h
   ! do ii = 1, natom
   ! do jj =1, 3
   !    kk=3*(ii-1)+jj
   !    write(777,*)grad0(kk) , DEner(kk,0)
   ! end do
   ! end do

end subroutine grad_calc

!###############################################################################

subroutine init_grad_calc(natom)
   use armonic_data    , only : move_atom, atom_mass, mass_w, emass_d,         &
                                armonic_vec, armonic_freq

   implicit none
   integer, intent(in) :: natom
   integer             :: ii


   allocate(move_atom(natom), atom_mass(natom), mass_w(natom),                 &
            armonic_freq(3*natom), armonic_vec(3*natom,3*natom))

!carlos: the name of the file is temporal, it should be dynamic

   move_atom = .false.
   atom_mass = 0.0d0

!carlos: the armonic.in file contains the mass information of each atom and also
!        the boolean information to allow or not the movement of the atom. The
!        mass is expresed in a.m.u.
   open(unit=10101, file = 'armonic.in')

   do ii = 1, natom
      read(10101,*) atom_mass(ii), move_atom(ii)
   end do

   close(10101)

!Converting the mass from amu to au
   ! atom_mass = atom_mass / emass_d

!Storing inverse square of the masses to weight the cartesian coordinates.
   do ii=1, natom
      mass_w(ii) = 1.0d0 / dsqrt(atom_mass(ii))
   end do

end subroutine init_grad_calc

!###############################################################################

subroutine build_hessian(natom, dhau, grad, hess_mat, mass_w)

   use armonic_data, only : hess_norder

   implicit none
   integer     , intent(in)  :: natom
   real(kind=8), intent(in)  :: grad(3*natom, -hess_norder : hess_norder,      &
                                3*natom)
   real(kind=8), intent(in)  :: mass_w(natom)
   real(kind=8), intent(in)  :: dhau
   real(kind=8), intent(out) :: hess_mat(3*natom, 3*natom)
   real(kind=8)              :: tmp1, tmp2, hhi, hhj
   integer                   :: ii, jj ,kk, ll

   if (hess_norder == 1) then
      kk = 1
      jj = 1
      do ii=1, 3*natom
      do jj=ii, 3*natom
         tmp1 =0.0d0
         tmp2 =0.0d0
         hhi  = dhau * mass_w(kk)
         hhj  = dhau * mass_w(ll)
         if (ii==jj) then
            tmp1=(grad(ii,1,jj)-grad(ii,-1,jj))*0.5d0/hhi!Finite difference.
            hess_mat(ii,jj)=tmp1 * mass_w(kk) * mass_w(ll)   !Mass weighting hessian
         else
            tmp1=(grad(ii,1,jj)-grad(ii,-1,jj))*0.5d0/hhi
            tmp2=(grad(jj,1,ii)-grad(jj,-1,ii))*0.5d0/hhj
            hess_mat(ii,jj)=(tmp1+tmp2)*0.5d0*mass_w(kk)*mass_w(ll)
            hess_mat(jj,ii)=hess_mat(ii,jj)
         end if
         if (mod(jj,3) == 0) ll = ll + 1
      end do
         if (mod(ii,3) == 0) then
            kk = kk + 1
            ll = kk
         end if
      end do

   else
      kk = 1
      jj = 1
      do ii=1, 3*natom
      do jj=ii,3*natom
         tmp1=0.0d0
         tmp2=0.0d0
         hhi  = dhau * mass_w(kk)
         hhj  = dhau * mass_w(ll)
         if (ii==jj) then
            tmp1=(grad(ii,-2,jj)+8d0*grad(ii,1,jj)-8d0*grad(ii,-1,jj)-          &
                  grad(ii,2,jj))/(12d0*hhi)            ! Finite difference.
            hess_mat(ii,jj)=tmp1 * mass_w(kk) * mass_w(ll) ! Mass weighting hessian
         else
            tmp1=(grad(ii,-2,jj)+8d0*grad(ii,1,jj)-8d0*grad(ii,-1,jj)-         &
                  grad(ii,2,jj))/(12d0*hhi)           ! Finite difference.
            tmp2=(grad(jj,-2,ii)+8d0*grad(jj,1,ii)-8d0*grad(jj,-1,ii)-         &
                  grad(jj,2,ii))/(12d0*hhj)           ! Finite difference.
            hess_mat(ii,jj)=(tmp1+tmp2)*0.5d0*mass_w(kk)*mass_w(ll)
            hess_mat(jj,ii)=hess_mat(ii,jj)
         end if
         if (mod(jj,3) == 0) ll = ll + 1
      end do
         if (mod(ii,3) == 0) then
            kk = kk + 1
            ll = kk
         end if
      end do

   end if

end subroutine build_hessian

!###############################################################################

subroutine diagonalize_hessian(natom, hess_mat, armonic_freq, armonic_vec)

   implicit none
   integer, intent(in)         :: natom
   real(kind=8), intent(in)    :: hess_mat(3*natom,3*natom)
   real(kind=8), intent(out)   :: armonic_freq(3*natom)
   real(kind=8), intent(out)   :: armonic_vec(3*natom,3*natom)


   integer     , allocatable    :: IWORK1(:), IWORK2(:)
   real(kind=8), allocatable    :: WORK1(:), WORK2(:)
   integer                      :: LWORK, LIWORK, INFO
   integer                      :: ncoords

   ncoords = 3*natom

   armonic_vec = hess_mat

   allocate ( WORK1(1000), IWORK1(1000) )
   LWORK=-1
   call dsyevd('V','U',ncoords,armonic_vec,ncoords,armonic_freq,WORK1,LWORK,   &
               IWORK1,LWORK,INFO)
   LWORK=WORK1(1)
   LIWORK=IWORK1(1)
   if(allocated(WORK2)) deallocate (WORK2,IWORK2)
   allocate (WORK2(LWORK),IWORK2(LIWORK))
   call dsyevd('V','U',ncoords,armonic_vec,ncoords,armonic_freq,WORK2,LWORK,   &
               IWORK2,LIWORK,INFO)

   deallocate( WORK1, IWORK1, WORK2, IWORK2 )

end subroutine diagonalize_hessian

!###############################################################################

end module armonic_subs
