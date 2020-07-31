#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module vib_KE_subs
   implicit none

contains
!###############################################################################
subroutine init_vib_calc(natom, M)
   use vib_KE_data    , only : move_atom, atom_mass, mass_w, emass_d,         &
                               armonic_vec, armonic_freq, tot_at_move

   implicit none
   integer, intent(in) :: natom
   integer, intent(in) :: M
   logical             :: vec_move(natom)
   logical             :: file_exists
   integer             :: ii, jj


   allocate(atom_mass(natom), mass_w(natom))

!carlos: the name of the file is temporal, it should be dynamic

   move_atom = .false.
   atom_mass = 0.0d0

!carlos: the vibration.in file contains the mass information of each atom and also
!        the boolean information to allow or not the movement of the atom. The
!        mass is expresed in a.m.u.

   inquire(file='vibration.in', exist=file_exists)
   if (file_exists) then
      open(unit=10101, file = 'vibration.in')
   else
      write(*,*) "File 'vibration.in' not found"
      write(*,*) "Closing program"
      stop
   end if

   tot_at_move = 0

   do ii = 1, natom
      read(10101,*) atom_mass(ii), vec_move(ii)
      if (vec_move(ii)) tot_at_move = tot_at_move + 1
   end do

   close(10101)

   allocate(move_atom(tot_at_move), armonic_freq(3*tot_at_move),               &
            armonic_vec(3*tot_at_move,3*tot_at_move), dS_dq())

   jj = 0
   do ii=1, natom
      if (vec_move(ii)) then
         jj = jj + 1
         move_atom(jj) = ii
      end if
   end do
!Storing inverse square of the masses to weight the cartesian coordinates.
   do ii=1, natom
      mass_w(ii) = 1.0d0 / dsqrt(atom_mass(ii))
   end do

end subroutine init_vib_calc

!###############################################################################
subroutine vibrational_calc(E, rho_aop, fock_aop,rho_bop, fock_bop )

   use garcha_mod      , only : r, natom, rqm, Smat
   use basis_data      , only : M, Nuc !nuc is temporary
   use vib_KE_data     , only : vib_calc, ke_calc, hess_norder, delta_h,       &
                                mass_w, armonic_freq, armonic_vec, freq2cm,    &
                                tot_at_move, move_atom
   use typedef_operator, only : operator

   implicit none
   type(operator), intent(inout)           :: rho_aop, fock_aop
   type(operator), intent(inout), optional :: rho_bop, fock_bop
   LIODBLE  , intent(inout)           :: E

   LIODBLE              :: dhau
   LIODBLE              :: r_init(natom,3)
   LIODBLE              :: r_new(natom,3)
   LIODBLE, allocatable :: grad0(:)
   LIODBLE, allocatable :: grad(:,:,:)
   LIODBLE, allocatable :: Df_Dr(:,:,:)
   LIODBLE, allocatable :: DS_dr(:,:,:)
   LIODBLE, allocatable :: Sinv(:,:)
   LIODBLE, allocatable :: dxyzqm(:,:)
   LIODBLE, allocatable :: Dfock_a(:,:,:)
   LIODBLE, allocatable :: hess_mat(:,:)
   LIODBLE, allocatable :: fock0(:,:)
   integer              :: ii, jj, kk, rr, ss, iat
   LIODBLE, allocatable :: freq(:)
   ! LIODBLE, allocatable :: DEner(:,:) !TEMPORAL ARRAY FOR TESTS
   ! LIODBLE, allocatable :: f_posta(:,:) !TEMPORAL ARRAY FOR TESTS

   if (.not.vib_calc) return

   call init_vib_calc(natom)

   dhau = delta_h / 0.529177D0

!Storin input coordinates
   r_init = r * 0.529177D0 ! r is in a.u. and r_init is in Amstrongs

   allocate(dxyzqm(3, natom), hess_mat(3*tot_at_move, 3*tot_at_move),          &
            freq(3*tot_at_move), grad0(3*tot_at_move))

   if (ke_calc == 1) allocate(Dfock_a(M,M,3*tot_at_move),fock0(M,M),           &
                              DS_dr(M,M,3*natom), Sinv(M,M))

   if (hess_norder==1) then
      allocate (grad(3*tot_at_move,-1:1,3*tot_at_move))
      if (ke_calc == 2) allocate(Df_Dr(M,M,-1:1))
   else
      allocate (grad(3*tot_at_move,-2:2,3*tot_at_move))
      if (ke_calc == 1) allocate(Df_Dr(M,M,-2:2))
   end if

   grad     = 0.0D0

! Calculating energy and gradient at the initial geometry

   call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
   call dft_get_qm_forces(dxyzqm)
   if (ke_calc == 1) then

      Df_Dr    = 0.0d0
      Dfock_a  = 0.0d0
      Sinv     = 0.0d0

      call fock_aop%Gets_data_AO(fock0)

      do ii = 1, natom
      do jj = 1, 3

         call calc_dSdR(DS_dr(M,M,jj+(ii-1)*3), r_init, ii, jj, M, natom)

      end do
      end do

   end if
   do ii =1, tot_at_move
   do jj=1, 3
      grad0(3*(ii-1)+jj) = dxyzqm(jj,move_atom(ii))
   end do
   end do


!  Calculating gradients in displaced positions along the
!  QM coordinates.

   do ii = 1, tot_at_move
   do jj=1, 3
      kk=3*(ii-1)+jj
      iat = move_atom(ii)

!     Forward displacement
!     --------------------
      r_new = r_init
      r_new(iat,jj) = r_init(iat,jj) + delta_h * mass_w(iat)
      r = r_new / 0.529177D0
!carlos: this equality is not necessary with the correct dimensions, and is
!        repeated in all the parts of the loop
      rqm = r
      call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
      call dft_get_qm_forces(dxyzqm)

      do rr=1,tot_at_move
      do ss=1,3
         grad(kk, 1, 3*(rr-1)+ss) = dxyzqm(ss,move_atom(rr))
      end do
      end do

      if (ke_calc == 1) call fock_aop%Gets_data_AO(Df_Dr(:,:,1))

!     Backward displacement
!     --------------------
      r_new = r_init
      r_new(iat,jj) = r_init(iat,jj) - delta_h * mass_w(iat)
      r = r_new / 0.529177D0
      rqm = r
      call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
      call dft_get_qm_forces(dxyzqm)

      do rr=1,tot_at_move
      do ss=1,3
         grad(kk, -1, 3*(rr-1)+ss) = dxyzqm(ss,move_atom(rr))
      end do
      end do

      if (ke_calc == 2) call fock_aop%Gets_data_AO(Df_Dr(:,:,-1))

      if(hess_norder == 2) then
!        Forward 2x displacement
!        -----------------------
         r_new = r_init
         r_new(iat,jj) = r_init(iat,jj) + 2.0d0 * delta_h * mass_w(iat)
         r = r_new / 0.529177D0
         rqm = r
         call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
         call dft_get_qm_forces(dxyzqm)

         do rr=1,tot_at_move
         do ss=1,3
            grad(kk, 2, 3*(rr-1)+ss) = dxyzqm(ss,move_atom(rr))
         end do
         end do

         if (ke_calc == 1) call fock_aop%Gets_data_AO(Df_Dr(:,:,2))
!        Backward 2x displacement
!        ------------------------
         r_new = r_init
         r_new(iat,jj) = r_init(iat,jj) - 2.0d0 * delta_h * mass_w(iat)
         r = r_new / 0.529177D0
         rqm = r
         call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
         call dft_get_qm_forces(dxyzqm)

         do rr=1,tot_at_move
         do ss=1,3
            grad(kk, -2, 3*(rr-1)+ss) = dxyzqm(ss,move_atom(rr))
         end do
         end do
         if (ke_calc == 1) call fock_aop%Gets_data_AO(Df_Dr(:,:,-2))
      endif

      if (ke_calc == 1) then
         if (hess_norder == 2) then
            Dfock_a(:,:,kk) = (8.0d0*Df_Dr(:,:,1) - 8.0d0*Df_Dr(:,:,-1) +      &
                               Df_Dr(:,:,-2) - Df_Dr(:,:,2))                   &
                              /(12.0d0 * dhau * mass_w(iat))
         else
            Dfock_a(:,:,kk) = (Df_Dr(:,:,1) - Df_Dr(:,:,-1))*0.5d0             &
                              /(dhau * mass_w(iat))
         end if
      end if

   end do
   end do

!Construction and diagonalization of the hessian matrix
   call build_hessian(natom, dhau, grad, hess_mat, mass_w)
   call diagonalize_hessian(tot_at_move, hess_mat, armonic_freq, armonic_vec)

!Freeing memory as this arrays are no more neccesary

   if (allocated(Df_Dr)) deallocate(Df_Dr)
   if (allocated(grad)) deallocate(grad)

!Changing units and printing frequencies:

   do ii = 1, 3*tot_at_move
      freq(ii) = sign(1.0d0,armonic_freq(ii)) * freq2cm *                &
                 sqrt(abs(armonic_freq(ii)))
   end do

end subroutine vibrational_calc
!###############################################################################
subroutine build_hessian(natom, dhau, grad, hess_mat, mass_w)

   use vib_KE_data, only : hess_norder, tot_at_move, move_atom

   implicit none
   integer     , intent(in)  :: natom
   real(kind=8), intent(in)  :: grad(3*tot_at_move, -hess_norder : hess_norder,&
                                3*natom)
   real(kind=8), intent(in)  :: mass_w(natom)
   real(kind=8), intent(in)  :: dhau
   real(kind=8), intent(out) :: hess_mat(3*tot_at_move, 3*tot_at_move)
   real(kind=8)              :: tmp1, tmp2, hhi, hhj
   integer                   :: ii, jj ,kk, ll

   if (hess_norder == 1) then
      kk = 1
      ll = 1
      do ii=1, 3*tot_at_move
      do jj=ii, 3*tot_at_move
         tmp1 =0.0d0
         tmp2 =0.0d0
         hhi  = dhau * mass_w(move_atom(kk))
         hhj  = dhau * mass_w(move_atom(ll))
         if (ii==jj) then
            tmp1=(grad(ii,1,jj)-grad(ii,-1,jj))*0.5d0/hhi!Finite difference.
            hess_mat(ii,jj)=tmp1 * mass_w(move_atom(kk)) * mass_w(move_atom(ll))
            !Mass weighting hessian
         else
            tmp1=(grad(ii,1,jj)-grad(ii,-1,jj))*0.5d0/hhi
            tmp2=(grad(jj,1,ii)-grad(jj,-1,ii))*0.5d0/hhj
            hess_mat(ii,jj)=(tmp1+tmp2)*0.5d0*mass_w(move_atom(kk)) *          &
                            mass_w(move_atom(ll))
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
      do ii=1, 3*tot_at_move
      do jj=ii,3*tot_at_move
         tmp1=0.0d0
         tmp2=0.0d0
         hhi  = dhau * mass_w(move_atom(kk))
         hhj  = dhau * mass_w(move_atom(ll))
         if (ii==jj) then
            tmp1=(grad(ii,-2,jj)+8d0*grad(ii,1,jj)-8d0*grad(ii,-1,jj)-          &
                  grad(ii,2,jj))/(12d0*hhi)            ! Finite difference.
            hess_mat(ii,jj)=tmp1 * mass_w(move_atom(kk)) * mass_w(move_atom(ll))
                                                       ! Mass weighting hessian
         else
            tmp1=(grad(ii,-2,jj)+8d0*grad(ii,1,jj)-8d0*grad(ii,-1,jj)-         &
                  grad(ii,2,jj))/(12d0*hhi)           ! Finite difference.
            tmp2=(grad(jj,-2,ii)+8d0*grad(jj,1,ii)-8d0*grad(jj,-1,ii)-         &
                  grad(jj,2,ii))/(12d0*hhj)           ! Finite difference.
            hess_mat(ii,jj)=(tmp1+tmp2)*0.5d0*mass_w(move_atom(kk))*           &
                            mass_w(move_atom(ll))
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


!###############################################################################
function gigj_int(ea, eb, ra, rb, iang, jang)

   implicit none
   LIODBLE             :: gigj_int
   integer, intent(in) :: iang(3), jang(3)
   LIODBLE, intent(in) :: ea, eb
   LIODBLE, intent(in) :: ra(3), rb(3)
   LIODBLE             :: Kab(3)
   LIODBLE             :: rp(3)
   LIODBLE             :: ep
   LIODBLE             :: Rab(3)
   LIODBLE             :: emu
   LIODBLE             :: int_omega
   LIODBLE             :: sqep
   LIODBLE             :: aux1
   LIODBLE, parameter  :: sqpi =  1.7724538509055160272981674833411
   integer             :: ii

   ep   = ea + eb
   sqep = dsqrt(ep)
   rp   = (ea*ra+eb*rb)/ep
   Rab  = ra - rb
   emu  = ea*eb/(ea+eb)
   Kab  = 0.0d0
   int_omega = 1.0d0

   do ii = 1, 3

      Kab(ii) = dexp(-emu * Rab(ii)**2.0d0)
      if (iang(ii)==0 .and. jang(ii)==0) then
         int_omega = int_omega * Kab(ii) * sqpi / sqep
      else if (iang(ii)==0 .and. jang(ii)==1) then
         int_omega = int_omega * Kab(ii) * sqpi / sqep * (rp(ii) - rb(ii))
      else if (iang(ii)==1 .and. jang(ii)==0) then
         int_omega = int_omega * Kab(ii) * sqpi / sqep * (rp(ii) - ra(ii))
      else if (iang(ii)==1 .and. jang(ii)==1) then
         aux1 = 1.0d0 + 2.0d0 * ep * (ra(ii) - rp(ii)) *(rb(ii) - rp(ii))
         int_omega = int_omega * Kab(ii) * sqpi / (2.0d0*sqep**3.0d0) * aux1
      else if (iang(ii)==0 .and. jang(ii)==2) then
         aux1 = 1.0d0 + 2.0d0 * ep * (rb(ii) - rp(ii))**2.0d0
         int_omega = int_omega * Kab(ii) * sqpi / (2.0d0*sqep**3.0d0) * aux1
      else if (iang(ii)==2 .and. jang(ii)==0) then
         aux1 = 1.0d0 + 2.0d0 * ep * (ra(ii) - rp(ii))**2.0d0
         int_omega = int_omega * Kab(ii) * sqpi / (2.0d0*sqep**3.0d0) * aux1
      else if (iang(ii)==1 .and. jang(ii)==2) then
         aux1 = 2.0d0 * ep * rp(ii) *(rb(ii) - rp(ii)**2.0d0) -                &
                ra(ii) * (1.0d0 + 2.0d0 * ep * (rb(ii)-rp(ii))**2.0d0) +       &
                3.0d0 * rp(ii) - 2.0d0 * rb(ii)
         int_omega = int_omega * Kab(ii) * sqpi / (2.0d0*sqep**3.0d0) * aux1
      else if (iang(ii)==2 .and. jang(ii)==1) then
         aux1 = 2.0d0 * ep * rp(ii) *(ra(ii) - rp(ii)**2.0d0) -                &
                rb(ii) * (1.0d0 + 2.0d0 * ep * (ra(ii)-rp(ii))**2.0d0) +       &
                3.0d0 * rp(ii) - 2.0d0 * ra(ii)
         int_omega = int_omega * Kab(ii) * sqpi / (2.0d0*sqep**3.0d0) * aux1
      else if (iang(ii)==2 .and. jang(ii)==2) then
         aux1 = 3.0d0 + 4.0d0 * ep**2.0d0 *                                    &
                ((ra(ii)-rp(ii))*(rb(ii)-rp(ii)))**2.0d0 + 2.0d0 * ep *        &
                (ra(ii)**2.0d0 + 4.0d0*ra(ii)*rb(ii)+rb(ii)**2.0d0 -           &
                 6.0d0 * (ra(ii)+rb(ii))*rp(ii) + 6.0d0 * rp(ii)**2.0d0)
         int_omega = int_omega * Kab(ii) * sqpi / (4.0d0*sqep**5.0d0) * aux1
      end if

   end do

   gigj_int = int_omega
end function gigj_int

subroutine calc_dSdR(dSdR_mat, rn, atom, axis, M_in, ntatom)

   use basis_data   , only: Nuc, a, c, ncont, NORM, M, nshell

   implicit none
   integer, intent(in)  :: M_in
   integer, intent(in)  :: atom
   integer, intent(in)  :: ntatom
   integer, intent(in)  :: axis
   LIODBLE, intent(in)  :: rn(ntatom, 3)
   LIODBLE, intent(out) :: dSdR_mat(M_in, M_in)
   integer              :: ifunct, jfunct, ns, np, nci, ncj, ati, atj, i_ind,  &
                           j_ind
   integer              :: iang(3), jang(3)
   integer              :: l1, l2
   LIODBLE              :: ccoef
   LIODBLE              :: aux1

   ns  = nshell(0); np = nshell(1)
   dSdR_mat = 0.0d0

   !this loop evaluate each element (i|j')

   l1 = 0 ; l2 = 0
   do ifunct = 1, ns+np
   do jfunct = 1, ns+np
      if (Nuc(jfunct) == atom) then

         if (ifunct>np .and. l1 < 3) then ; l1 = l1 + 1
         else ; l1 = 0
         end if
         if (jfunct>np .and. l2 < 3) then ; l2 = l2 + 1
         else ; l2 = 0
         end if

         ati = Nuc(ifunct)
         atj = Nuc(jfunct)
         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            iang(:) = 0
            jang(:) = 0
            if (l1 /= 0) iang(l1) = 1
            if (l2 /= 0) jang(l2) = 1
            jang(axis) = jang(axis) + 1

            ccoef   = c(ifunct,nci) * 2.0d0*a(jfunct,ncj)*c(jfunct,ncj)

            aux1    = gigj_int(a(ifunct,nci), a(jfunct,ncj), rn(ati,:),     &
                               rn(atj,:), iang, jang)
            dSdR_mat(i_ind,jfunct) = dSdR_mat(i_ind,jfunct) + ccoef * aux1

            if (axis == l2) then
               jang(axis) = jang(axis) - 2
               ccoef      = - c(ifunct,nci) * c(jfunct,ncj)
               aux1       = gigj_int(a(ifunct,nci), a(jfunct,ncj), rn(ati,:),  &
                                     rn(atj,:), iang, jang)
               dSdR_mat(i_ind,jfunct) = dSdR_mat(i_ind,jfunct) + ccoef * aux1
            end if

         end do
         end do
      end if
   end do
   end do

end subroutine calc_dSdR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function inv_mat(A) result(Ainv)
    implicit none
    LIODBLE,intent(in):: A(:,:)
    LIODBLE           :: Ainv(size(A,1),size(A,2))
    LIODBLE           :: work(size(A,1))            ! work array for LAPACK
    integer           :: n,info,ipiv(size(A,1))     ! pivot indices

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
end function inv_mat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

end module vib_KE_subs
