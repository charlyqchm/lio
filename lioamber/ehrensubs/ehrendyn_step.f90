!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!carlos: I added nb_ext as a new variale for the extra basis, if TB is off
!        nb_ext=nbasis
subroutine ehrendyn_step( step_back, propagator_id, time_value, time_step,     &
                        & nbasis, nb_ext, natoms, nucpos, nucvel, nucfor_ds,   &
                        & Sinv, Uinv, Linv, dens_oldi, dens_midi, dens_newo,   &
                        & fock_mid, dipmom, energy )
!------------------------------------------------------------------------------!
!
! All matrices must be in OM
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   logical, intent(inout)    :: step_back
   integer, intent(in)       :: propagator_id
   real*8 , intent(in)       :: time_value
   real*8 , intent(in)       :: time_step
   integer, intent(in)       :: nbasis
   integer, intent(in)       :: nb_ext
   integer, intent(in)       :: natoms

   real*8, intent(in)        :: nucpos(3, natoms)
   real*8, intent(in)        :: nucvel(3, natoms)
   real*8, intent(inout)     :: nucfor_ds(3, natoms)

   real*8, intent(in)        :: Sinv(nbasis, nbasis)
   real*8, intent(in)        :: Uinv(nb_ext, nb_ext)
   real*8, intent(in)        :: Linv(nb_ext, nb_ext)

   complex*16, intent(in)    :: dens_oldi(nb_ext, nb_ext)
   complex*16, intent(in)    :: dens_midi(nb_ext, nb_ext)
   complex*16, intent(out)   :: dens_newo(nb_ext, nb_ext)

   real*8    , intent(inout) :: fock_mid(nb_ext, nb_ext)
   real*8    , intent(inout) :: dipmom(3)
   real*8    , intent(inout) :: energy

   complex*16, allocatable   :: dens_old(:,:)
   complex*16, allocatable   :: dens_mid(:,:)
   real*8    , allocatable   :: fock0(:,:)
   real*8                    :: dipmom0(3)
   real*8                    :: energy0

   integer                   :: substep, substeps, nn
   real*8                    :: time, dt
   real*8                    :: elec_field(3)
   real*8    , allocatable   :: Bmat(:,:)
   real*8    , allocatable   :: Dmat(:,:)
   complex*16, allocatable   :: Tmat(:,:)
   complex*16, allocatable   :: dens_mao(:,:)

   integer, parameter :: propagator_id_verlet = 1
   integer, parameter :: propagator_id_magnus = 2

   allocate( Bmat(nbasis,nbasis), Dmat(nbasis,nbasis), Tmat(nb_ext,nb_ext) )
   allocate( dens_old(nb_ext,nb_ext), dens_mid(nb_ext,nb_ext) )
   allocate( dens_mao(nb_ext,nb_ext) )
   allocate( fock0(nb_ext,nb_ext) )

   dens_old   = dens_oldi
   dens_mid   = dens_midi
   fock0      = fock_mid
   dipmom0(:) = 0.0d0
   energy0    = energy

   time = time_value
   if (step_back) then
      substeps = 20
      dt = (-time_step) / ( (2.0d0)*(substeps) )
      dens_old = dens_mid
   else
      substeps = 0
      dt = time_step
   endif


   do substep = 0, substeps
      fock_mid = fock0
      dipmom   = dipmom0
      energy   = energy0
      nucfor_ds(:,:) = 0.0d0

!     Preparing matrices (received in AO and returned in ON)
      call ehrendyn_prep( Nbasis, nb_ext, Natoms, time, .true., Uinv, Linv,    &
                          Sinv, dens_mid, fock_mid, Bmat, Dmat, Tmat,          &
                          nucpos, nucvel, nucfor_ds, dipmom, energy )

!     Density Propagation (works in ON)
      if (propagator_id==propagator_id_verlet) then
!         print*,'This is verlet!'
         call ehrenaux_verlet( nb_ext, dt, Tmat, dens_old, dens_mid, dens_newo )

      else if (propagator_id==propagator_id_magnus) then
!         print*,'This is magnus!'
         call ehrenaux_magnus( nbasis, 20, dt, Tmat, dens_mid, dens_newo )

      else
         print*,'Unidentified substep!'; stop

      endif

      if (step_back) then
         step_back = .false.
         dt = (-2.0d0) * dt
         dens_old = dens_newo
      else
         dens_old = dens_mid
         dens_mid = dens_newo
      endif

   enddo

   deallocate( Bmat, Dmat, Tmat )
   deallocate( dens_old, dens_mid )
   deallocate( dens_mao )
   deallocate( fock0 )

end subroutine ehrendyn_step
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
