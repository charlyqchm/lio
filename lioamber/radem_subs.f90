#include "complex_type.fh"
module radem_subs
   implicit none
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine radi_emission_init(M_in, M, open_shell, r, d, natom, ntatom)
   use radem_data, only: radiative_calc, dip_radi_op, ddip_radi_op,            &
                         fock_radi_op, d2dip_radi, Xmat_radi, k_radi
   use faint_cpu , only: intfld

   implicit none
   logical     , intent(in)  :: open_shell
   integer     , intent(in)  :: M_in
   integer     , intent(in)  :: M
   integer     , intent(in)  :: natom
   integer     , intent(in)  :: ntatom
   real(kind=8), intent(in)  :: r(ntatom,3)
   real(kind=8), intent(in)  :: d(natom,natom)
   integer                   :: ii
   integer                   :: MM
   real(kind=8)              :: vec_aux(3)
   real(kind=8), allocatable :: dip_array(:)
   real(kind=8), allocatable :: dip_mat_aux(:,:)

   if (.not. radiative_calc) return

   MM = M*(M+1)/2

   open(unit=2002, file="accelerator.in")
   read(2002,*) k_radi
   close(2002)

   allocate(dip_radi_op(3),dip_array(MM), dip_mat_aux(M, M))

   if (.not.open_shell) then
      allocate(fock_radi_op(1), ddip_radi_op(3,1), d2dip_radi(3,1))
   else
      allocate(fock_radi_op(2), ddip_radi_op(3,2), d2dip_radi(3,2))
   end if

   do ii=1, 3
      dip_array    = 0.0d0
      dip_mat_aux  = 0.0d0
      vec_aux      = 0.0d0
      vec_aux(ii)  = 1.0d0
      call intfld(dip_array, dip_array, r, d, natom, ntatom, .false., 1.0d0,   &
                  vec_aux(1), vec_aux(2), vec_aux(3))
      call spunpack('L', M, dip_array, dip_mat_aux)
!carlols: dip_mat is stored as the orthonormal part of the operator to be used
!         with commutator later.
      call dip_radi_op(ii)%Sets_data_AO(dip_mat_aux)
      call dip_radi_op(ii)%BChange_AOtoON(Xmat_radi, M_in)

   end do

end subroutine radi_emission_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine radi_fock_calculation(fock, rho_aux, M, dim3, t_step)
   use radem_data, only: radiative_calc, dip_radi_op, ddip_radi_op,            &
                         d2dip_radi, preA_radi, fock_radi_op, Ymat_radi, k_radi
   implicit none

   integer     , intent(in)    :: M
   integer     , intent(in)    :: dim3
   integer     , intent(in)    :: t_step
   real(kind=8), intent(in)    :: fock(M,M)
   TDCOMPLEX   , intent(inout) :: rho_aux(M,M)
   integer                     :: ii, jj, kk
   real(kind=8)                :: aux1
   real(kind=8)                :: aux_mat1(M,M)
   real(kind=8)                :: aux_mat2(M,M)
   real(kind=8)                :: aux_mat3(M,M)
   TDCOMPLEX                   :: aux_mat4(M,M)

   if (.not.radiative_calc) return

   aux_mat3 = 0.0d0
   aux_mat4 = 0.0d0

   if (t_step > 5000) then
      do ii=1,3
         call dip_radi_op(ii)%Commut_data_r(fock, aux_mat1, M)
         call ddip_radi_op(ii,dim3)%Sets_data_ON(aux_mat1)
         call ddip_radi_op(ii,dim3)%Commut_data_r(fock, aux_mat2, M)
         d2dip_radi(ii,dim3) = 0.0d0
         do jj=1, M
         do kk=1, M
            aux1 = dble(rho_aux(jj,kk)*aux_mat2(kk,jj))
            d2dip_radi(ii,dim3) = d2dip_radi(ii,dim3) - aux1
         end do
         end do

         aux_mat3 = aux_mat3 + k_radi * preA_radi * d2dip_radi(ii,dim3) * aux_mat1
      end do

      call fock_radi_op(dim3)%Sets_data_ON(aux_mat3)
      call fock_radi_op(dim3)%Commut_data_c(rho_aux, aux_mat4, M)

   end if

   rho_aux = aux_mat4


end subroutine radi_fock_calculation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

end module radem_subs
