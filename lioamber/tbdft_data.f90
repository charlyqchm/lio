module tbdft_data
   implicit none

   logical      :: tbdft_calc =.false.      ! Logical indicator for tbdft calculation
   integer      :: MTB    = 0               ! Size of the two tight-binding subatrices
   integer      :: n_atTB = 0
   integer      :: MTBDFT = 0               ! Size of the DFT-TB matrix
   integer      :: start_tdtb=0             ! Initial time step for evolution of diagonal TB terms
   integer      :: end_tdtb=0               ! Final time step for evolution of diagonal TB terms
   integer      :: end_bTB                  ! Index matrix size
   integer      :: n_biasTB                 ! Number of electrodes
   integer      :: n_atperbias              ! Number of atoms per bias
   integer    , allocatable :: Iend_TB(:,:) ! Index matrix
   integer    , allocatable :: linkTB(:,:)  ! Link atoms, separated by bias
   integer    , allocatable :: basTB(:)     ! Coupling basis in the LIO order
   real(kind=8) :: alfaTB                   ! Fermi Energy
   real(kind=8) :: betaTB                   ! Offdiagonal tight binding param
   real(kind=8) :: gammaTB                  ! DFT-TB terms
   real(kind=8) :: Vbias_TB                 ! Bias potential
   real(kind=8) :: chargeA_TB
   real(kind=8) :: chargeB_TB
   real(kind=8) :: chargeM_TB
   real(kind=8)   , allocatable :: rhoa_TBDFT(:,:)     ! Matrix to store rho TBDFT for TD
   real(kind=8)   , allocatable :: rhob_TBDFT(:,:)     ! Matrix to store rho TBDFT for TD
   real(kind=8)   , allocatable :: chimerafock (:,:,:) ! Allocated in the central code
   real(kind=8)   , allocatable :: gammaW(:)         ! gamma weight, per atom
   real(kind=8)   , allocatable :: VbiasTB(:)          ! Bias potential for each
                                                       ! electrode
#ifdef TD_SIMPLE
   complex(kind=4), allocatable :: rhold_AOTB(:,:,:)   ! rho in AO to calculate charges
   complex(kind=4), allocatable :: rhonew_AOTB(:,:,:)  ! rho in AO to calculate charges
#else
   complex(kind=8), allocatable :: rhold_AOTB(:,:,:)   ! rho in AO to calculate charges
   complex(kind=8), allocatable :: rhonew_AOTB(:,:,:)  ! rho in AO to calculate charges
#endif

!DLVN-TB variables:
   integer                      :: tbdft_transport = 0
   real(kind=8)                 :: driving_rateTB = 0.00d0
   complex(kind=4), allocatable :: rhofirst_TB(:,:,:)

!Carlos Temporal variables for OM base change:

   real*8, allocatable  :: rho_OM(:,:), rho_OM2(:,:)
   real*8, allocatable  :: auto_vec(:,:), auto_vec_t(:,:)
   real*8, allocatable  :: auto_inv(:,:), auto_t_inv(:,:)
   ! real*8, allocatable  :: Umat(:,:), Umat_t(:,:)
   integer, allocatable :: tr_ind(:) !indice de transporte
   real*8, allocatable  :: auto_val(:)
   real*8, allocatable  :: temp_TB(:)

end module tbdft_data
