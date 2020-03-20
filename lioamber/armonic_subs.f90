module armonic_subs
   implicit none

contains

!###############################################################################
subroutine grad_calc(E, rho_aop, fock_aop,rho_bop, fock_bop )

   use garcha_mod      , only : r, natom, rqm, Smat
   use basis_data      , only : M, Nuc !nuc is temporary
   use armonic_data    , only : armonic_calc, hess_norder, delta_h, mass_w,    &
                                armonic_freq, armonic_vec, freq2cm
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
   real(kind=8), allocatable :: DS_Dr(:,:,:)
   real(kind=8), allocatable :: dxyzqm(:,:)
   real(kind=8), allocatable :: Dfock_a(:,:,:)
   real(kind=8), allocatable :: DS_mat(:,:,:)
   real(kind=8), allocatable :: hess_mat(:,:)
   real(kind=8), allocatable :: S_inv(:,:)
   real(kind=8), allocatable :: fock0(:,:)
   integer                   :: ii, jj, kk, rr, ss
   real(kind=8)              :: freq(3*natom)
   ! real(kind=8), allocatable :: DEner(:,:) !TEMPORAL ARRAY FOR TESTS
   ! real(kind=8), allocatable :: f_posta(:,:) !TEMPORAL ARRAY FOR TESTS

   if (armonic_calc == 0) return

   call init_grad_calc(natom)

   dhau = delta_h / 0.529177D0

!Storin input coordinates
   r_init = r * 0.529177D0 ! r is in a.u. and r_init is in Amstrongs

   allocate(dxyzqm(3, natom), hess_mat(3*natom, 3*natom))

   if (armonic_calc == 2) allocate(Dfock_a(M,M,3*natom),DS_mat(M,M,3*natom),   &
                                   S_inv(M,M), fock0(M,M))

   if (hess_norder==1) then
      allocate (grad(3*natom,-1:1,3*natom))
      if (armonic_calc == 2) allocate(Df_Dr(M,M,-1:1), DS_Dr(M,M,-1:1))
   else
      allocate (grad(3*natom,-2:2,3*natom))
      if (armonic_calc == 2) allocate(Df_Dr(M,M,-2:2), DS_Dr(M,M,-2:2))
   end if

   grad     = 0.0D0
   Df_Dr    = 0.0d0
   Dfock_a  = 0.0d0
   S_inv    = 0.0d0

! Calculating energy and gradient at the initial geometry

   call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
   call dft_get_qm_forces(dxyzqm)
   if (armonic_calc == 2) then
      call invert_mat(Smat, S_inv, M)
      call fock_aop%Gets_data_AO(fock0)
   end if
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

      if (armonic_calc == 2) then
         call fock_aop%Gets_data_AO(Df_Dr(:,:,1))
         DS_Dr(:,:,1) = Smat
      end if

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
         grad(kk, -1, 3*(rr-1)+ss) = dxyzqm(ss,rr)
      end do
      end do

      if (armonic_calc == 2) then
         call fock_aop%Gets_data_AO(Df_Dr(:,:,-1))
         DS_Dr(:,:,-1) = Smat
      end if

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
            grad(kk, 2, 3*(rr-1)+ss) = dxyzqm(ss,rr)
         end do
         end do

         if (armonic_calc == 2) then
            call fock_aop%Gets_data_AO(Df_Dr(:,:,2))
            DS_Dr(:,:,2) = Smat
         end if
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
            grad(kk, -2, 3*(rr-1)+ss) = dxyzqm(ss,rr)
         end do
         end do
         if (armonic_calc == 2) then
            call fock_aop%Gets_data_AO(Df_Dr(:,:,-2))
            DS_Dr(:,:,-2) = Smat
         end if
      endif

      if (armonic_calc == 2) then
         if (hess_norder == 2) then
            Dfock_a(:,:,kk) = (8.0d0*Df_Dr(:,:,1) - 8.0d0*Df_Dr(:,:,-1) +      &
                               Df_Dr(:,:,-2) - Df_Dr(:,:,2))                   &
                              /(12.0d0 * dhau * mass_w(ii))
            DS_mat(:,:,kk) = (8.0d0*DS_Dr(:,:,1) - 8.0d0*DS_Dr(:,:,-1) +      &
                               DS_Dr(:,:,-2) - DS_Dr(:,:,2))                   &
                              /(12.0d0 * dhau * mass_w(ii))
         else
            Dfock_a(:,:,kk) = (Df_Dr(:,:,1) - Df_Dr(:,:,-1))*0.5d0             &
                              /(dhau * mass_w(ii))
            DS_mat(:,:,kk) = (DS_Dr(:,:,1) - DS_Dr(:,:,-1))*0.5d0             &
                              /(dhau * mass_w(ii))
         end if
      end if

   end do
   end do

!Construction and diagonalization of the hessian matrix
   call build_hessian(natom, dhau, grad, hess_mat, mass_w)
   call diagonalize_hessian(natom, hess_mat, armonic_freq, armonic_vec)

!Freeing memory as this arrays are no more neccesary

   if (allocated(DS_Dr)) deallocate(DS_Dr)
   if (allocated(Df_Dr)) deallocate(Df_Dr)
   if (allocated(grad)) deallocate(grad)

!Changing units and printing frequencies:

   do ii = 1, 3*natom
      freq(ii) = sign(1.0d0,armonic_freq(ii)) * freq2cm *                &
                 sqrt(abs(armonic_freq(ii)))
   end do

end subroutine grad_calc

!###############################################################################

subroutine init_grad_calc(natom)
   use armonic_data    , only : move_atom, atom_mass, mass_w, emass_d,         &
                                armonic_vec, armonic_freq

   implicit none
   integer, intent(in) :: natom
   logical             :: file_exists
   integer             :: ii


   allocate(move_atom(natom), atom_mass(natom), mass_w(natom),                 &
            armonic_freq(3*natom), armonic_vec(3*natom,3*natom))

!carlos: the name of the file is temporal, it should be dynamic

   move_atom = .false.
   atom_mass = 0.0d0

!carlos: the armonic.in file contains the mass information of each atom and also
!        the boolean information to allow or not the movement of the atom. The
!        mass is expresed in a.m.u.

   inquire(file='armonic.in', exist=file_exists)
   if (file_exists) then
      open(unit=10101, file = 'armonic.in')
   else
      write(*,*) "File 'armonic.in' not found"
      write(*,*) "Closing program"
      stop
   end if

   do ii = 1, natom
      read(10101,*) atom_mass(ii), move_atom(ii)
   end do

   close(10101)

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
      ll = 1
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
         end if
         ll = kk
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
         ll  = kk
         if (mod(ii,3) == 0) then
            kk = kk + 1
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
!scaling
   armonic_vec = hess_mat * 1000d0

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
!descaling
   armonic_freq = armonic_freq/1000d0
   armonic_vec  = armonic_vec/1000d0

   deallocate( WORK1, IWORK1, WORK2, IWORK2 )

end subroutine diagonalize_hessian

!###############################################################################
!Invertion subroutine
subroutine invert_mat(A, Ainv,M)

  implicit none
  integer, intent(in):: M
  real*8, intent(in) :: A(M,M)
  real*8, intent(inout) :: Ainv(M,M)

  real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end subroutine invert_mat

!###############################################################################

subroutine build_elc_phon_coupling(natom, M, Dfock_a, DS_mat, S_inv, fock0,    &
                                   armonic_vec, Nuc)
   implicit none
   integer     , intent(in)    :: natom
   integer     , intent(in)    :: M
   integer     , intent(in)    :: Nuc(M)
   real(kind=8), intent(inout) :: Dfock_a(M,M, 3*natom)
   real(kind=8), intent(in)    :: DS_mat(M,M, 3*natom)
   real(kind=8), intent(in)    :: armonic_vec(3*natom,3*natom)
   real(kind=8), intent(in)    :: S_inv(M,M)
   real(kind=8), intent(in)    :: fock0(M,M)

   integer      :: ii, jj, kk, ll, aa, rr
   real(kind=8) :: Df_armonic(M,M, 3*natom)
   real(kind=8) :: Df_coupling(M,M,3*natom)

   Df_armonic  = 0.0d0
   Df_coupling = 0.0d0

   aa =1
   do rr = 1, 3*natom
   do jj = 1, M
   do ii = 1, M
!carlos: el segundo loop solo barrerlo si es necesario
      Df_coupling(ii,jj,rr) = Dfock_a(ii,jj,rr)
!we are assuming that dSij/Ra is equal to <i'|j> + <i|j'> = 0 =>
! <i'|j> = - <i|j'> = 0
      if (Nuc(ii) == aa .and. Nuc(jj) == aa ) then
         do kk = 1, M
         do ll = 1, M
            Df_coupling(ii,jj,rr) = Df_coupling(ii,jj,rr) -                    &
                                    DS_mat(ii,kk)*S_inv(kk,ll)*fock0(ll,jj) -  &
                                    fock0(ii,kk)*S_inv(kk,ll)*DS_mat(ll,jj)
         end do
         end do
      else if(Nuc(ii) == aa) then
         do kk = 1, M
         do ll = 1, M
            Df_coupling(ii,jj,rr) = Df_coupling(ii,jj,rr) -                    &
                                    DS_mat(ii,kk)*S_inv(kk,ll)*fock0(ll,jj)
         end do
         end do
      else if(Nuc(jj) == aa) then
         do kk = 1, M
         do ll = 1, M
            Df_coupling(ii,jj,rr) = Df_coupling(ii,jj,rr) -                    &
                                    fock0(ii,kk)*S_inv(kk,ll)*DS_mat(ll,jj)
         end do
         end do
      end if
   end do
   end do
      if (mod(aa,3)) aa = aa + 1
   end do

end subroutine build_elc_phon_coupling()

!###############################################################################
end module armonic_subs
