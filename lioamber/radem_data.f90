#include "complex_type.fh"
module radem_data
   use typedef_operator, only: operator
   use typedef_cumat   , only: cumat_x, cumat_r

   implicit none

   logical                     :: radiative_calc = .true.
   real(kind=8)                :: preA_radi = 2.592668788247536d-7!2.592668788247536d-7
                                 !2/3 \hbar \alpha/c^2
   real(kind=8)  , allocatable :: d2dip_radi(:,:)
   type(cumat_x)               :: Ymat_radi
   type(cumat_r)               :: Xmat_radi
   type(operator), allocatable :: dip_radi_op(:)
   type(operator), allocatable :: ddip_radi_op(:,:)
   type(operator), allocatable :: fock_radi_op(:)

!charly:variables cochinas:
   real(kind=8), allocatable :: coef_r(:,:), coefT_r(:,:)
   real(kind=8), allocatable :: cinv_r(:,:), cinvT_r(:,:)
   real(kind=8), allocatable :: dip_vec(:)
   integer                   :: npoints = 20

end module radem_data
