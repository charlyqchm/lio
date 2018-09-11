module DOS_data
   implicit none
   logical               :: dos_calc     = .true.
   logical               :: pdos_calc    = .true.
   logical               :: pdos_allb    = .false.
   integer               :: dos_nsteps   = 50000
   real*8                :: dos_sigma    = 0.001d0
   real*8                :: dos_Eref     = -0.166063225427645
   integer               :: pdos_nbases  = 0            !Option for just 1 atom
   integer               :: pdos_natoms  = 0
   integer, allocatable  :: pdos_nuc(:)
   integer, allocatable  :: pdos_base(:)
   real*8,  allocatable  :: pdos(:)
   real*8,  allocatable  :: pdos_b(:,:)


end module DOS_data
