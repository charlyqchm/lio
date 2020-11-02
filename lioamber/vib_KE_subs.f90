#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module vib_KE_subs
   implicit none

contains
!###############################################################################
subroutine init_vib_calc(natom, M)
   use vib_KE_data    , only : move_atom, atom_mass, mass_w, armonic_vec,      &
                               armonic_freq, nat_move, ke_coef, ke_eorb,       &
                               ke_calc

   implicit none
   integer, intent(in) :: natom
   integer, intent(in) :: M
   logical             :: vec_move(natom)
   logical             :: file_exists
   integer             :: ii, jj


   allocate(atom_mass(natom), mass_w(natom))
   if (ke_calc==1) allocate(ke_coef(M,M), ke_eorb(M))

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

   nat_move = 0

   do ii = 1, natom
      read(10101,*) atom_mass(ii), vec_move(ii)
      if (vec_move(ii)) nat_move = nat_move + 1
   end do

   close(10101)

   allocate(move_atom(nat_move), armonic_freq(3*nat_move),               &
            armonic_vec(3*nat_move,3*nat_move))

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

   use garcha_mod      , only : r, natom, rqm, Smat, MO_coef_at, Eorbs
   use basis_data      , only : M
   use vib_KE_data     , only : vib_calc, ke_calc, hess_norder, delta_h,       &
                                mass_w, armonic_freq, armonic_vec, freq2cm,    &
                                nat_move, move_atom, ke_coef, ke_eorb
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
   LIODBLE, allocatable :: hess_mat(:,:), hessp_mat(:,:)
   LIODBLE, allocatable :: fock0(:,:)
   integer              :: ii, jj, kk, rr, ss, iat
   LIODBLE, allocatable :: freq(:)
   ! LIODBLE, allocatable :: DEner(:,:) !TEMPORAL ARRAY FOR TESTS
   ! LIODBLE, allocatable :: f_posta(:,:) !TEMPORAL ARRAY FOR TESTS

   if (.not.vib_calc) return

   call init_vib_calc(natom, M)

   dhau = delta_h / 0.529177D0

!Storin input coordinates
   r_init = r * 0.529177D0 ! r is in a.u. and r_init is in Amstrongs

   allocate(dxyzqm(3, natom), hess_mat(3*nat_move, 3*nat_move),                &
            freq(3*nat_move), grad0(3*nat_move),                               &
            hessp_mat(3*nat_move, 3*nat_move),Dfock_a(1,1,1),                  &
            fock0(1,1), DS_dr(1,1,1), Sinv(1,1), Df_Dr(1,1,1),                 &
            grad(3*nat_move,-1:1,3*nat_move))

   if (ke_calc == 1) then
      deallocate(Dfock_a,fock0,DS_dr,Sinv,Df_Dr)
      allocate(Dfock_a(M,M,3*nat_move),fock0(M,M),DS_dr(M,M,3*nat_move),       &
               Sinv(M,M))
      if (hess_norder==1) then ; allocate(Df_Dr(M,M,-1:1))
      else ; allocate(Df_Dr(M,M,-2:2))
      end if
   end if

   if (hess_norder==2) then
      deallocate(grad)
      allocate (grad(3*nat_move,-2:2,3*nat_move))
   end if

   grad     = 0.0D0

! Calculating energy and gradient at the initial geometry

   call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
   call dft_get_qm_forces(dxyzqm)

!Initialization of matrices to store fock derivatives and dS/dq calculation
   if (ke_calc == 1) then

      Df_Dr    = 0.0d0
      Dfock_a  = 0.0d0
      Sinv     = 0.0d0

      call fock_aop%Gets_data_AO(fock0)

      ke_coef = MO_coef_at
      ke_eorb = Eorbs

      do ii = 1, nat_move
      do jj = 1, 3
         call calc_dSdR(DS_dr(:,:,jj+(ii-1)*3), r, move_atom(ii), jj, M,  &
                        natom)
         DS_dr(:,:,jj+(ii-1)*3) = DS_dr(:,:,jj+(ii-1)*3) / mass_w(move_atom(ii))
      end do
      end do

      Sinv = inv_mat(Smat)

      do ii= 1, M
      do jj= 1, M
      do kk= 1, 3*nat_move
         write(888,*) DS_dr(ii,jj,kk)
      end do
      end do
      end do

   end if

   do ii =1, nat_move
   do jj=1, 3
      grad0(3*(ii-1)+jj) = dxyzqm(jj,move_atom(ii))
   end do
   end do

!  Calculating gradients in displaced positions along the
!  QM coordinates.

   do ii = 1, nat_move
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

      do rr=1,nat_move
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

      do rr=1,nat_move
      do ss=1,3
         grad(kk, -1, 3*(rr-1)+ss) = dxyzqm(ss,move_atom(rr))
      end do
      end do

      if (ke_calc == 1) call fock_aop%Gets_data_AO(Df_Dr(:,:,-1))

      if(hess_norder == 2) then
!        Forward 2x displacement
!        -----------------------
         r_new = r_init
         r_new(iat,jj) = r_init(iat,jj) + 2.0d0 * delta_h * mass_w(iat)
         r = r_new / 0.529177D0
         rqm = r
         call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
         call dft_get_qm_forces(dxyzqm)

         do rr=1,nat_move
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

         do rr=1,nat_move
         do ss=1,3
            grad(kk, -2, 3*(rr-1)+ss) = dxyzqm(ss,move_atom(rr))
         end do
         end do
         if (ke_calc == 1) call fock_aop%Gets_data_AO(Df_Dr(:,:,-2))
      endif
      !Building fock derivatives wei by the mass
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
   call eckart(hess_mat,r_init,mass_w,move_atom,natom,nat_move,hessp_mat)
   hess_mat = hessp_mat
   call diagonalize_hessian(nat_move, hess_mat, armonic_freq, armonic_vec)

!Calculating and storing dH/dQ only for ke_calc = 1
   call calc_dHdQ(Dfock_a, fock0, Sinv, DS_dr, armonic_vec, M, nat_move)

!Changing units and printing frequencies:

   write(*,201)
   write(*,202)
   write(*,203)
   do ii = 1, 3*nat_move
      freq(ii) = sign(1.0d0,armonic_freq(ii)) * freq2cm *                &
                 sqrt(abs(armonic_freq(ii)))
      write(*,204) ii, freq(ii)
   end do
   write(*,205)

!Freeing memory as this arrays are no more neccesary

   if (allocated(Df_Dr)) deallocate(Df_Dr)
   if (allocated(grad))  deallocate(grad)
   if (allocated(grad0)) deallocate(grad0)
   if (allocated(DS_dr)) deallocate(DS_dr)
   if (allocated(Sinv)) deallocate(Sinv)
   if (allocated(dxyzqm)) deallocate(dxyzqm)
   if (allocated(Dfock_a)) deallocate(Dfock_a)
   if (allocated(hess_mat)) deallocate(hess_mat)
   if (allocated(hessp_mat)) deallocate(hessp_mat)
   if (allocated(fock0)) deallocate(fock0)
   if (allocated(freq)) deallocate(freq)
   call vib_ke_finalize()

201 FORMAT(2x,"--------------------------------")
202 FORMAT(2x,"| VIBRATIONAL MODES            |")
203 FORMAT(2x,"--------------------------------")
204 FORMAT(2x,I5,"   ",F14.8," cm-1")
205 FORMAT(2x,"--------------------------------")

end subroutine vibrational_calc
!###############################################################################
subroutine eckart(hess,rc,mass_w,move_atom,natom,nat_move,hessp)

   use vib_KE_data, only: emass_d

   implicit none
!    ------------------------------------------------------------------
   integer,intent(in)  :: natom
   integer,intent(in)  :: nat_move
   integer,intent(in)  :: move_atom(nat_move)
   LIODBLE,intent(in)  :: mass_w(natom) ! Atomic numbers of QM atoms.
   LIODBLE,intent(in)  :: hess(nat_move*3,nat_move*3)
   LIODBLE,intent(in)  :: rc(3,natom)
   LIODBLE,intent(out) :: hessp(nat_move*3,nat_move*3)
   !    Internal variables
   integer             :: ii,jj,iat
   integer             :: ndf
   LIODBLE             :: Mass(nat_move*3)
   LIODBLE             :: X0(3,nat_move)
   LIODBLE             :: Itsr(3,3)
   LIODBLE             :: innp(3,3)
   LIODBLE             :: outp(3,3)
   LIODBLE             :: ai(3)
   LIODBLE             :: totM
   LIODBLE             :: P(nat_move*3,nat_move*3)
   LIODBLE             :: tmp(nat_move*3,nat_move*3)
   LIODBLE             :: cmass(3)
   LIODBLE             :: TRmat(nat_move*3,6)
   LIODBLE             :: unitv(3,3)
   LIODBLE             :: SV(6)
   LIODBLE             :: massau

   LIODBLE, external    :: ddot

   LIODBLE, allocatable :: VT(:,:),U(:,:)

   LIODBLE, allocatable :: WORK(:),IWORK(:)
   integer              :: LWORK,INFO
   !    ------------------------------------------------------------------

   ndf   = nat_move*3
   cmass = 0.0d0

   unitv=reshape((/ 1d0, 0d0, 0d0, &
               &   0d0, 1d0, 0d0, &
               &   0d0, 0d0, 1d0 /),(/3,3/))

   !    Converting initial coordinates to AU
   do ii=1, nat_move
   do jj=1, 3
      X0(jj,ii)=rc(move_atom(ii),jj)/0.5291771D0
   end do
   end do
   !    Computing total mass and atomic mass array and center of mass.
   totM=0.0d0
   do ii=1,nat_move
      do jj=1,3
          massau = (dsqrt(emass_d)/mass_w(move_atom(ii)))**2.0d0
          Mass(3*(ii-1)+jj) = dsqrt(massau)
          cmass(jj)=cmass(jj)+massau*X0(jj,ii)
   !             version with atomic masses in AU
   !             Mass(3*(i-1)+j) = sqrt_atomic_masses_au(k)
   !             cmass(j)=cmass(j)+atomic_masses_au(k)*X0(j,i)
      end do
      totM=totM+massau
   end do
   cmass = cmass/totM

!    Translating to center of mass and mass-weighting.
   do ii=1,nat_move
     do jj=1,3
        X0(jj,ii)=(X0(jj,ii)-cmass(jj))*Mass(3*(ii-1)+jj)
     end do
   end do


   !    Moment of inertia tensor.
   !    Itsr = Sum_i [ai**T.ai - ai.ai**T]
   Itsr=0.0d0
   do ii=1,nat_move
     ai=X0(:,ii)
     innp=0.0d0
     outp=0.0d0
     innp(1,1) = ddot(3,ai,1,ai,1)
     innp(2,2) = innp(1,1)
     innp(3,3) = innp(1,1)
     call dger(3,3,1d0,ai,1,ai,1,outp,3)
     Itsr=Itsr+innp-outp
   end do

   !    Symmetrizing inertia tensor.
   do ii=1,2
     do jj=ii+1,3
       Itsr(jj,ii)=Itsr(ii,jj)
     end do
   end do

   !    WE NOW COMPUTE THE PROJECTOR MATRIX P
   !    THE FORMULA CAN BE FOUND IN
   !    Szalay, JCP 140:234107, (2014), eq 9, and eq. 17
   !    Note that although eq. 12 is only an orthonormalized version of
   !    eq. 9. Also, we only implemented Eckart conditions. Sayvetz
   !    conditions are left for another time.

   !    Also, the following paper may be of interest.
   !    Miller, Handy and Adams, JCP 72:99,(1980) eq 4.11

   !    The imposition of Eckart conditions is made by projection.

   !                  PHess = P.Hess.P

   !    Initializing TRmat
   Trmat=0d0

   !    Computing translation and rotation vectors.
   do ii=1,3
     do iat=1,nat_move
        TRmat(3*(iat-1)+ii,ii)=Mass(3*(iat-1)+ii)/Sqrt(totM)  ! Translation
        TRmat(3*(iat-1)+ii,4:6) = crossp(X0(:,iat),unitv(ii,:)) ! Rotation
     end do
   end do

   !    Orthonormalization of translation and rotation vectors by singular value decomposition method.
   LWORK=-1
   allocate(VT(6,6),U(0,0),WORK(1000),IWORK(8*6))
   call dgesdd('O',ndf,6,TRmat,ndf,SV,U,ndf,VT,6,WORK,LWORK,IWORK,INFO)
   LWORK=int(WORK(1))
   deallocate(WORK)
   allocate(WORK(LWORK))
   call dgesdd('O',ndf,6,TRmat,ndf,SV,U,ndf,VT,6,WORK,LWORK,IWORK,INFO)
   deallocate(WORK,IWORK,U,VT)
   if (INFO /= 0) STOP ('Error during singular value decomposition in eckart subroutine')

   !    Building projection matrix by external product of TRmat.
   call dgemm('n','t',ndf,ndf,6,1d0,TRmat,ndf,TRmat,ndf,0d0,P,ndf)

   P=-P
   do ii=1,ndf
     P(ii,ii) = 1.0d0 + P(ii,ii)
   end do

   !    PROJECTION
   !    HESSP = (1-P).HESS.(1-P)
   hessp=0.0d0
   call dsymm('L','U',ndf,ndf,1d0,hess,ndf,P,ndf,0d0,tmp,ndf)
   call dgemm('N','N',ndf,ndf,ndf,1d0,P,ndf,tmp,ndf,0d0,hessp,ndf)
   do ii=1,ndf
     do jj=1,ndf
        if (abs(hessp(ii,jj)) < 1d-10) hessp(ii,jj)=0d0
     end do
   end do

end subroutine

!###############################################################################
subroutine build_hessian(natom, dhau, grad, hess_mat, mass_w)

   use vib_KE_data, only : hess_norder, nat_move, move_atom

   implicit none
   integer     , intent(in)  :: natom
   real(kind=8), intent(in)  :: grad(3*nat_move, -hess_norder : hess_norder,&
                                3*natom)
   real(kind=8), intent(in)  :: mass_w(natom)
   real(kind=8), intent(in)  :: dhau
   real(kind=8), intent(out) :: hess_mat(3*nat_move, 3*nat_move)
   real(kind=8)              :: tmp1, tmp2, hhi, hhj
   integer                   :: ii, jj ,kk, ll

   if (hess_norder == 1) then
      kk = 1
      ll = 1
      do ii=1, 3*nat_move
      do jj=ii, 3*nat_move
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
      ll = 1
      do ii=1, 3*nat_move
      do jj=ii,3*nat_move
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
   LWORK=int(WORK1(1))
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

!###############################################################################
subroutine calc_dSdR(dSdR_mat, rn, atom, axis, M_in, ntatom)

   use basis_data   , only: Nuc, a, c, ncont, nshell

   implicit none
   integer, intent(in)  :: M_in
   integer, intent(in)  :: atom
   integer, intent(in)  :: ntatom
   integer, intent(in)  :: axis
   LIODBLE, intent(in)  :: rn(ntatom, 3)
   LIODBLE, intent(out) :: dSdR_mat(M_in, M_in)
   integer              :: ifunct, jfunct, ns, np, nci, ncj, ati, atj
   integer              :: iang(3), jang(3)
   integer              :: l1, l2
   LIODBLE              :: ccoef
   LIODBLE              :: aux1

   ns  = nshell(0); np = nshell(1)
   dSdR_mat = 0.0d0

   !this loop evaluate each element (i|j')

   l1 = 0
   do ifunct = 1, ns+np
      if (ifunct>ns .and. l1 < 3) then ; l1 = l1 + 1
      else if (ifunct>ns) then; l1 = 1
      end if
   l2 = 0
   do jfunct = 1, ns+np
      if (Nuc(jfunct) == atom) then
         if (jfunct>ns .and. l2 < 3) then ; l2 = l2 + 1
         else if (jfunct>ns) then ; l2 = 1
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
            dSdR_mat(ifunct,jfunct) = dSdR_mat(ifunct,jfunct) + ccoef * aux1

            if (axis == l2) then
               jang(axis) = jang(axis) - 2
               ccoef      = - c(ifunct,nci) * c(jfunct,ncj)
               aux1       = gigj_int(a(ifunct,nci), a(jfunct,ncj), rn(ati,:),  &
                                     rn(atj,:), iang, jang)
               dSdR_mat(ifunct,jfunct) = dSdR_mat(ifunct,jfunct) + ccoef * aux1
            end if

         end do
         end do
      end if
   end do
   end do

end subroutine calc_dSdR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine calc_dHdQ(Dfock, fock0, Sinv, dSdR, Lmat, M, nat_move)
   use vib_KE_data, only: ke_calc, XCmat_ke, ke_coef, n_vib, armonic_freq,     &
                          ke_eorb

   integer, intent(in) :: M
   integer, intent(in) :: nat_move
   LIODBLE, intent(in) :: Dfock(M,M,3*nat_move)
   LIODBLE, intent(in) :: dSdR(M,M,3*nat_move)
   LIODBLE, intent(in) :: fock0(M,M)
   LIODBLE, intent(in) :: Sinv(M,M)
   LIODBLE, intent(in) :: Lmat(3*nat_move,3*nat_move)
   LIODBLE             :: dHdQ(M,M,3*nat_move)
   LIODBLE             :: dH_mat(M,M,3*nat_move)
   LIODBLE             :: freq(3*nat_move)
   integer             :: ii, jj, kk, ll, rr, ind_i

   if (ke_calc /=1) return

   dH_mat = 0.0d0
   dHdQ   = 0.0d0
   call XCmat_ke%init(M, ke_coef)

   do rr=1, nat_move*3
   do ii=1, M
   do jj=1, M
      dH_mat(ii,jj,rr) = Dfock(ii,jj,rr)
      write(777,*) dH_mat(ii,jj,rr)
      do kk=1, M
      do ll=1, M
         dH_mat(ii,jj,rr) = dH_mat(ii,jj,rr) -                                 &
                            dSdR(kk,ii,rr) * Sinv(kk,ll) * fock0(ll,jj) -      &
                            fock0(ii,kk) * Sinv(kk,ll) * dSdR(ll,jj,rr)
      end do
      end do
   end do
   end do
   end do

   do jj=1, nat_move*3
   do ii=1, nat_move*3
      dHdQ(:,:,jj) = dHdQ(:,:,jj) + dH_mat(:,:,ii)*Lmat(ii,jj)
   end do
   end do

   do ii=1, nat_move*3
      call XCmat_ke%change_base(dHdQ(:,:,ii), 'dir')
   end do

   n_vib = 0
   do ii = 1, 3*nat_move
      freq(ii) = sign(1.0d0,armonic_freq(ii)) * dsqrt(dabs(armonic_freq(ii)))* &
                 0.023421782790063343d0
      if (freq(ii) > 1.0d-05) n_vib = n_vib + 1
   end do

!All the data for electron-phonon calculations is stored in elecphon.out
   open(unit=10101, file = 'elecphon.out')

   write(10101,*) n_vib
   ind_i = 3*nat_move - n_vib + 1
   do ii = ind_i, 3*nat_move
      write(10101,*) freq(ii)
   end do

   do kk= ind_i, 3*nat_move
   do jj=1, M
   do ii=1, M
      write(10101,*) dHdQ(ii,jj,kk)
      write(999,*) dHdQ(ii,jj,kk)
   end do
   end do
   end do

   do ii=1, M
      write(10101,*) ke_eorb(ii)
   end do

   do jj=1, M
   do ii=1, M
      write(10101,*) ke_coef(ii,jj)
   end do
   end do

   close(10101)

   write(*,301)

301 FORMAT(2x,"File elecphon.out has been succesfully created")
end subroutine calc_dHdQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine init_KE_evol(M)
   use vib_KE_data, only: armonic_freq, n_vib, ke_eorb, ke_calc, YCinv_ke,     &
                          phon_pop, ke_coef, phon_pop, PFW_vec, phon_temp,     &
                          ke_lmin, ke_lmax, XCmatC_ke, gamma_mat

   integer, intent(in)  :: M
   logical              :: file_exists
   integer              :: ii, jj, kk, len_PFW
   LIODBLE              :: coef_inv(M,M)
   LIODBLE, allocatable :: dHdQ(:,:,:)
   TDCOMPLEX            :: coef_cmplx(M,M)
   TDCOMPLEX            :: liocmplx


   if (.not.(ke_calc==2)) return

!Selecting relevant levels
   if(ke_lmin == 0) ke_lmin = 1
   if(ke_lmax == 0) ke_lmax = M

   allocate(ke_eorb(M), ke_coef(M,M), gamma_mat(M,M))

!elecphon.in must contain the number of active vibrations and the coupling terms
   inquire(file='elecphon.in', exist=file_exists)
   if (file_exists) then
      open(unit=10101, file = 'elecphon.in')
   else
      write(*,*) "File 'elecphon.in' not found"
      write(*,*) "Closing program"
      stop
   end if

   read(10101,*) n_vib
   allocate(dHdQ(M,M,n_vib), armonic_freq(n_vib), phon_pop(n_vib))
   do ii=1, n_vib
      read(10101,*) armonic_freq(ii)
   end do

   do kk=1, n_vib
   do jj=1, M
   do ii=1, M
      read(10101,*) dHdQ(ii,jj,kk)
   end do
   end do
   end do

   do ii=1, M
      read(10101,*) ke_eorb(ii)
   end do

   do jj=1, M
   do ii=1, M
      read(10101,*) ke_coef(ii,jj)
   end do
   end do

   close(10101)
!YCinv stores the amtrix YC' to proyect \rho in the molecular orbital base.
   coef_inv = inv_mat(ke_coef)
   coef_inv = transpose(coef_inv)

   do jj=1, M
   do ii=1, M
      coef_cmplx(ii,jj) = liocmplx(coef_inv(ii,jj),0.0d0)
   end do
   end do

   call YCinv_ke%init(M,coef_cmplx)

   do jj=1, M
   do ii=1, M
      coef_cmplx(ii,jj) = liocmplx(ke_coef(ii,jj),0.0d0)
   end do
   end do

   call XCmatC_ke%init(M, coef_cmplx)

   call neglect_terms(dHdQ, M)
   call create_phon_bath(phon_pop, armonic_freq ,phon_temp, n_vib)

   deallocate(dHdQ, ke_coef)

   len_PFW = size(PFW_vec)

   write(*,101)
   write(*,102)
   write(*,103) n_vib
   write(*,104)
   do ii=1,n_vib
      write(*,105) ii, armonic_freq(ii)
   end do
   write(*, 106) len_PFW, M*M*n_vib
   write(*,101)


101 FORMAT(2x,"-----------------------------------------------------")
102 FORMAT(2x,"ELECTRON-PHONON COUPLING SUCCESSFULLY INITIATED")
103 FORMAT(2x,"Number of active vibrations =", I5)
104 FORMAT(2x,"Vibrational mode - Frequencie:")
105 FORMAT(2x, I5," - ",F14.8," a.u.")
106 FORMAT(2x,"Active electron-phone couplings:",I5," of ",I5)

end subroutine init_KE_evol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine neglect_terms(dHdQ, M)
   use vib_KE_data, only: armonic_freq, n_vib, ke_eorb, PFW_vec, ke_sigma,     &
                          ke_ind, ke_tol, ke_ka, ke_lmax, ke_lmin, ke_degen,   &
                          gamma_mat

   integer, intent(in) :: M
   LIODBLE, intent(in) :: dHdQ(M,M,n_vib)
   integer             :: count
   integer             :: ii, jj, kk, ll, M2
   LIODBLE             :: aux1, aux2, wdos
   LIODBLE             :: Ea, Eb, wj, Fj
   LIODBLE, parameter  :: pi =  3.141592653589793

   M2        = M*M
   count     = 0
   gamma_mat = 0.0d0

   do ii=1, M-1
   do jj=ii+1, M
      aux1 = 0.0d0
      aux2 = 0.0d0
      wdos = 0.0d0
   do kk=1, n_vib
      Ea   = ke_eorb(ii)
      Eb   = ke_eorb(jj)
      wj   = armonic_freq(kk)
      Fj   = dHdQ(ii,jj,kk)**2.0d0

      aux1 = aux1  + pi * Fj * dirac_delta(wj,Eb-Ea,ke_tol)/wj
      aux2 = aux2  + dirac_delta(wj,Eb-Ea,ke_tol)
      wdos = wdos + dirac_delta(Eb-Ea,wj,ke_sigma)
      ! if((aux1>ke_tol.or.aux2>ke_tol).and.(ii/=jj).and.                        &
      !    (dabs(Ea-Eb)>ke_degen)) then
      ! if((ii>=ke_lmin.and.jj>=ke_lmin).and.(ii<=ke_lmax.and.jj<=ke_lmax)) then
      !    count = count + 1
      ! end if
      ! end if

   end do
      if (dabs(Ea-Eb)>ke_degen) then
         gamma_mat(ii,jj) = ke_ka * wdos * aux1 / aux2
         gamma_mat(jj,ii) = gamma_mat(ii,jj)
      end if
   end do
   end do
   ! if (count > 0) then
   !    allocate(PFW_vec(count), ke_ind(count))
   ! else
   !    write(*,*) "PFW_vec can't be created, please encrease the broadening"
   !    write(*,*) "of the Dirac delta or the threshold."
   !    stop
   ! end if

   write(*,*) "ELECTRON-PHONON TERMS ii jj gamma"

   do ii=1, M-1
   do jj=ii+1, M
      write(*,*) ii, jj, gamma_mat(ii,jj)
   end do
   end do

   ! ll = 1
   ! do ii=1, M
   ! do jj=1, M
   ! do kk=1, n_vib
   !    Ea   = ke_eorb(ii)
   !    Eb   = ke_eorb(jj)
   !    wj   = armonic_freq(kk)
   !    Fj   = dHdQ(ii,jj,kk)**2.0d0
   !    aux1 = pi * Fj * dirac_delta(Ea,Eb,wj,ke_sigma)/wj
   !    aux2 = pi * Fj * dirac_delta(Ea,Eb,-wj,ke_sigma)/wj
   !
   !
   !    if((aux1>ke_tol.or.aux2>ke_tol).and.(ii/=jj).and.                        &
   !       (dabs(Ea-Eb)>ke_degen)) then
   !    if((ii>=ke_lmin.and.jj>=ke_lmin).and.(ii<=ke_lmax.and.jj<=ke_lmax)) then
   !       PFW_vec(ll) = ke_ka * pi * Fj/wj
   !       ke_ind(ll)   = (ii-1)+M*(jj-1)+M2*(kk-1)
   !       write(*,*) ii, jj, kk, PFW_vec(ll)
   !       ll = ll + 1
   !    end if
   !    end if
   !
   ! end do
   ! end do
   ! end do



end subroutine neglect_terms
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ke_rho_evolve(rho_at, M, istep)
   use  vib_KE_data, only: PFW_vec, ke_sigma, phon_pop, ke_eorb,        &
                           armonic_freq, YCinv_ke, ke_calc, ke_ind, ke_sigma,  &
                           YCinv_ke, XCmatC_ke, phon_temp, gamma_mat

   implicit none
   integer  , intent(in)    :: M, istep
   TDCOMPLEX, intent(inout) :: rho_at(M,M)
   TDCOMPLEX                :: rho_ke(M,M)
   integer                  :: ii, jj, kk, ll, M2
   integer                  :: len_PFW
   LIODBLE                  :: eta(M)
   LIODBLE                  :: lambda(M)
   LIODBLE                  :: Ea, Eb, wj, rhoa, rhob, Nj, exp1, exp2
   TDCOMPLEX                :: rho_OM(M,M)
   TDCOMPLEX                :: traza, liocmplx

   if (.not.(ke_calc==2)) return

   M2      = M*M
   rho_OM  = rho_at
   rho_ke  = liocmplx(0.0d0,0.0d0)
   len_PFW = size(PFW_vec)
   eta     = 0.0d0
   lambda  = 0.0d0

   call YCinv_ke%change_base(rho_OM, 'dir')

   traza = 0.0d0
   do ii=1,M
      traza = traza+rho_OM(ii,ii)
   end do
      write(*,*) "Traza=", real(traza)

   !write(*,*) "HOMO ", real(rho_OM(21,21))
   !write(*,*) "LUMO ", real(rho_OM(22,22))
   !write(*,*) "LUMO+1 ", real(rho_OM(23,23))



   ! if (mod(istep, 100) == 0) then
   !    do ii=1,M
   !       write(777,*) real(rho_OM(ii,ii))
   !    end do
   ! end if

   ! do ll=1, len_PFW
   !    kk = int(ke_ind(ll)/M2)+1
   !    jj = int((ke_ind(ll)-M2*(kk-1))/M)+1
   !    ii = ke_ind(ll)-M2*(kk-1)-M*(jj-1)+1
   !
   !    Ea   = ke_eorb(ii)
   !    Eb   = ke_eorb(jj)
   !    wj   = armonic_freq(kk)
   !    exp1 = dirac_delta(Ea,Eb,wj,ke_sigma)
   !    exp2 = dirac_delta(Ea,Eb,-wj,ke_sigma)
   !    Nj   = phon_pop(kk)
   !    rhoa = dble(rho_OM(ii,ii))/2.0d0
   !    rhob = dble(rho_OM(jj,jj))/2.0d0
   !
   !    eta(ii)   = eta(ii) + PFW_vec(ll)*((Nj+rhob)*exp1+(Nj-rhob+1)*exp2)
   !    lambda(ii)= lambda(ii) + PFW_vec(ll)*rhob*((Nj+1)*exp1 + Nj*exp2)
   ! end do

   do ii=1, M
   do jj=1, M
      Ea   = ke_eorb(ii)
      Eb   = ke_eorb(jj)
      wj   = dabs(Eb-Ea)
      Nj   = 1.0d0 / (dexp(wj/phon_temp)-1.0d0)
      rhoa = dble(rho_OM(ii,ii))/2.0d0
      rhob = dble(rho_OM(jj,jj))/2.0d0

      if (Eb>Ea) then
         eta(ii)    = eta(ii)    + gamma_mat(ii,jj) * (Nj+rhob)
         lambda(ii) = lambda(ii) + gamma_mat(ii,jj) * rhob*(Nj+1)
      else if(Eb<Ea) then
         eta(ii)    = eta(ii)    + gamma_mat(ii,jj) * (Nj-rhob+1)
         lambda(ii) = lambda(ii) + gamma_mat(ii,jj) * rhob*Nj
      end if

   end do
   end do

   do jj=1, M
   do ii=1, M
      if (ii==jj) then
         rho_ke(ii,ii) = real(-eta(ii), COMPLEX_SIZE/2) * rho_OM(ii,ii) +       &
                         real(2.0d0*lambda(ii), COMPLEX_SIZE/2)
      else
         rho_ke(ii,jj) = real(-0.5d0*(eta(ii)+eta(jj)), COMPLEX_SIZE/2) *       &
                         rho_OM(ii,jj)
      end if
   end do
   end do

   traza = 0.0d0
   do ii=1,M
      traza = traza+rho_ke(ii,ii)
   end do
      write(*,*) "KE traza=", real(traza)

   call XCmatC_ke%change_base(rho_ke, 'inv')
   rho_at = rho_ke

end subroutine ke_rho_evolve
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine create_phon_bath(phon_pop, freq ,phon_temp, n_vib)

   implicit none
   integer, intent(in)  :: n_vib
   LIODBLE, intent(in)  :: phon_temp
   LIODBLE, intent(in)  :: freq(n_vib)
   LIODBLE, intent(out) :: phon_pop(n_vib)
   integer              :: ii

   do ii=1, n_vib
      phon_pop(ii) = 1.0d0 / (dexp(freq(ii)/phon_temp)-1.0d0)
      ! write(*,*) "Phonon", ii, phon_pop(ii)
   end do

end subroutine create_phon_bath
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine vib_ke_finalize()
   use vib_KE_data, only: move_atom, atom_mass, mass_w, armonic_freq,          &
                          armonic_vec, ke_coef, ke_eorb, PFW_vec, ke_ind,      &
                          phon_pop, vib_calc, XCmat_ke, YCinv_ke, XCmatC_ke
   implicit none
   if (.not.vib_calc) return

   if (allocated(move_atom))    deallocate(move_atom)
   if (allocated(atom_mass))    deallocate(atom_mass)
   if (allocated(mass_w))       deallocate(mass_w)
   if (allocated(armonic_freq)) deallocate(armonic_freq)
   if (allocated(armonic_vec))  deallocate(armonic_vec)
   if (allocated(ke_coef))      deallocate(ke_coef)
   if (allocated(ke_eorb))      deallocate(ke_eorb)
   if (allocated(PFW_vec))      deallocate(PFW_vec)
   if (allocated(ke_ind))       deallocate(ke_ind)
   if (allocated(phon_pop))     deallocate(phon_pop)
   call XCmat_ke%destroy()
   call YCinv_ke%destroy()
   call XCmatC_ke%destroy()

end subroutine vib_ke_finalize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function dirac_delta(xi,mu,sigma) result(dd)
   implicit none
   LIODBLE, intent(in) :: xi, mu, sigma
   LIODBLE             :: norm, arg
   LIODBLE             :: dd
   LIODBLE, parameter  :: pi =  3.141592653589793

   ! norm = 1.0d0 / dsqrt(2.0d0 * pi * sigma**2.0d0)
   ! arg  = - (xi-mu)**2.0d0/(2.0d0 * sigma**2.0d0)
   ! dd   = norm * dexp(arg)
   dd = 1/pi * sigma/((xi-mu)**2+sigma**2)
end function dirac_delta

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
function crossp(v1,v2)
!       ---------------------------------------------------------------------
!       Computes cross product v1 x v2 = vr, with all vectors dimension 3
!       ---------------------------------------------------------------------
   LIODBLE  :: v1(3), v2(3)
   LIODBLE  :: crossp(3)

   crossp(1) = v1(2)*v2(3)-v1(3)*v2(2)
   crossp(2) = v1(3)*v2(1)-v1(1)*v2(3)
   crossp(3) = v1(1)*v2(2)-v1(2)*v2(1)

   return
end function crossp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

end module vib_KE_subs
