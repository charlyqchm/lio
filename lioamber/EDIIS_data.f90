module ediis_data

   implicit none

   integer             :: nediis = 70
   integer             :: step_nediis=0
   real*8              :: damp_ediis = 0.8d0
   real*8, allocatable :: fock_ediis_mat(:,:,:,:)
   real*8, allocatable :: rho_ediis_mat(:,:,:,:)
   real*8, allocatable :: BMAT(:,:,:)
   real*8, allocatable :: EDIIS_E(:)
   real*8, allocatable :: fock_damp(:,:,:)
end module ediis_data
