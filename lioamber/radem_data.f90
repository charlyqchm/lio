#include "complex_type.fh"
module radem_data
   use typedef_operator, only: operator
   use typedef_cumat   , only: cumat_x

   implicit none

   logical                     :: radiative_calc = .false.
   real(kind=8)                :: preA_radi = 2.592668788247536d-7!2.592668788247536d-7
                                 !2/3 \hbar \alpha/c^2
   real(kind=8)  , allocatable :: d2dip_radi(:,:)
   type(cumat_x)               :: Ymat_radi
   type(operator), allocatable :: dip_radi_op(:)
   type(operator), allocatable :: ddip_radi_op(:,:)
   type(operator), allocatable :: fock_radi_op(:)

end module radem_data
