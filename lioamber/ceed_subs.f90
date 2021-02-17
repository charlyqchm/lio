#include "datatypes/datatypes.fh"
module ceed_subs
   implicit none
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ceed_init(M, open_shell, r, d, natom, ntatom, propagator, coef_at)
! This subroutine initialize the variables for CEED calculations
   use ceed_data, only: ceed_calc, dip_ceed_op, d2ip_ceed_op,                  &
                        fock_ceed_op, d2dip_ceed, Xmat_ceed, d2mu_vec,         &
                        Xmat_ceed, YCinv_ceed, Xtrans_ceed, C_ON_mat_ceed,     &
                        CR_ON_mat_ceed
   use faint_cpu , only: intfld


   implicit none
   logical, intent(in)  :: open_shell
   integer, intent(in)  :: M
   integer, intent(in)  :: propagator
   integer, intent(in)  :: natom
   integer, intent(in)  :: ntatom
   LIODBLE, intent(in)  :: r(ntatom,3)
   LIODBLE, intent(in)  :: d(natom,natom)
   LIODBLE, intent(in)  :: coef_at(M,M)
   integer              :: ii,jj
   integer              :: MM
   LIODBLE              :: vec_aux(3)
   LIODBLE              :: aux_mat1(M,M)
   TDCOMPLEX            :: aux_Cmat1(M,M), aux_Cmat2(M,M)
   LIODBLE, allocatable :: dip_array(:)
   LIODBLE, allocatable :: dip_mat_aux(:,:)
   TDCOMPLEX            :: liocmplx

   if (ceed_calc==0) return

   ! if (propagator /= 1) then
   !    write(*,*) "CEED calculation can only be performed with Verlet propagator"
   !    write(*,*) "Stopping LIO"
   !    stop
   ! end if

   MM = M*(M+1)/2

   allocate(dip_ceed_op(3), dip_array(MM), dip_mat_aux(M, M))

   if (.not.open_shell) then
      allocate(fock_ceed_op(1), d2ip_ceed_op(3,1), d2dip_ceed(3,1))
      if (ceed_calc==2) allocate(d2mu_vec(M, 1))
   else
      allocate(fock_ceed_op(2), d2ip_ceed_op(3,2), d2dip_ceed(3,2))
      if (ceed_calc==2) allocate(d2mu_vec(M, 2))
   end if

   if (ceed_calc == 2) then
      call Xmat_ceed%init(M, coef_at)
      aux_mat1 = inv_mat(coef_at)
      aux_mat1 = transpose(aux_mat1)

      do jj=1, M
      do ii=1, M
         aux_Cmat1(ii,jj) = liocmplx(aux_mat1(ii,jj),0.0d0)
      end do
      end do
      call YCinv_ceed%init(M,aux_Cmat1)
      call Xtrans_ceed%multiply(aux_Cmat2, aux_Cmat1)
      call C_ON_mat_ceed%init(M, aux_Cmat2)

      do jj=1, M
      do ii=1, M
         aux_mat1(ii,jj) = dble(aux_Cmat2(ii,jj))
      end do
      end do

      call CR_ON_mat_ceed%init(M,aux_mat1)

   end if

   do ii=1, 3
      dip_array    = 0.0d0
      dip_mat_aux  = 0.0d0
      vec_aux      = 0.0d0
      vec_aux(ii)  = 1.0d0
      call intfld(dip_array, dip_array, r, d, natom, ntatom, .false., 1.0d0,   &
                  vec_aux(1), vec_aux(2), vec_aux(3))
      call spunpack('L', M, dip_array, dip_mat_aux)
      call dip_ceed_op(ii)%Sets_data_AO(dip_mat_aux)
      call dip_ceed_op(ii)%BChange_AOtoON(Xmat_ceed, M)
   end do

   deallocate(dip_array, dip_mat_aux)
   call Xtrans_ceed%destroy()

end subroutine ceed_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ceed_fock_calculation(fock_aop, rho_aop, M, t_step, dim3,           &
                                 open_shell, rho_aux, fock_bop, rho_bop)
!This subroutine includes the CEED term in the evolution of the density matrix.
!At the end of the subroutine rho_aux store the commutator -i[H_CEED,\rho],
!were:
!              H_CEED =  -i A_ceed d^2<\mu>/dt^2 . [\mu , H]
!d^2<\mu>/dt^2 is approximated as:
!               d^2<\mu>/dt^2 = - Tr(\rho . [[\mu,H],H])
!and:
!             A_ceed = 2/3 \alpha/c^2 = 2.592668788247536d-7

   use ceed_data, only: ceed_calc, dip_ceed_op, d2ip_ceed_op, d2dip_ceed,     &
                         A_ceed, fock_ceed_op, k_ceed, ceed_td_step, d2mu_vec, &
                         Xmat_ceed, YCinv_ceed, C_ON_mat_ceed, CR_ON_mat_ceed
   use typedef_operator, only: operator
   implicit none

   type(operator), intent(inout)           :: rho_aop, fock_aop
   type(operator), intent(inout), optional :: rho_bop, fock_bop
   logical   , intent(in)    :: open_shell
   integer   , intent(in)    :: M
   integer   , intent(in)    :: dim3
   integer   , intent(in)    :: t_step
   TDCOMPLEX , intent(inout) :: rho_aux(M,M,dim3)
   LIODBLE                   :: fock(M,M,dim3)
   integer                :: ii, jj, kk, ss
   LIODBLE                :: aux1, aux2
   LIODBLE                :: aux_mat1(M,M,dim3)
   LIODBLE                :: aux_mat2(M,M,dim3)
   LIODBLE                :: aux_mat3(M,M,dim3)
   TDCOMPLEX              :: aux_mat4(M,M,dim3)
   LIODBLE                :: aux_mat5(M,M,dim3)
   TDCOMPLEX              :: liocmplx

   if (ceed_calc==0) return


   if(ceed_calc == 1) then
      call fock_aop%Gets_data_ON(fock(:,:,1))
      call rho_aop%Gets_dataC_ON(rho_aux(:,:,1))
      if (open_shell) then
         call fock_bop%Gets_data_ON(fock(:,:,2))
         call rho_bop%Gets_dataC_ON(rho_aux(:,:,2))
      end if
   else if (ceed_calc == 2) then
      call fock_aop%Gets_data_ON(fock(:,:,1))
      call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
      call CR_ON_mat_ceed%change_base(fock(:,:,1), 'dir')
      call YCinv_ceed%change_base(rho_aux(:,:,1), 'dir')
      if (open_shell) then
         call fock_bop%Gets_data_ON(fock(:,:,2))
         call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))
         call CR_ON_mat_ceed%change_base(fock(:,:,2), 'dir')
         call YCinv_ceed%change_base(rho_aux(:,:,2), 'dir')
      end if
   end if

   aux_mat3 = 0.0d0
   aux_mat4 = liocmplx(0.0d0,0.0d0)

   if (t_step > ceed_td_step) then
      do ii = 1,3
         if(ceed_calc==2) d2mu_vec = 0.0d0
      do ss = 1,dim3

         call dip_ceed_op(ii)%Commut_data_r(fock(:,:,ss), aux_mat1(:,:,ss), M)
         call d2ip_ceed_op(ii,ss)%Sets_data_ON(aux_mat1(:,:,ss))
         call d2ip_ceed_op(ii,ss)%Commut_data_r(fock(:,:,ss),aux_mat2(:,:,ss),M)
         d2dip_ceed(ii,ss) = 0.0d0
         do jj=1, M
         do kk=1, M
            aux1 = dble(rho_aux(jj,kk,ss)*aux_mat2(kk,jj,ss))
            d2dip_ceed(ii,ss) = d2dip_ceed(ii,ss) - aux1
            if(ceed_calc==2) d2mu_vec(jj,ss) = d2mu_vec(jj,ss) - aux1
         end do
         end do
!###############################################################################
!testeando dipolo
         ! aux2 = 0.0d0
         !
         ! if (ii==1) then
         !    call dip_ceed_op(ii)%Gets_data_ON(aux_mat5(:,:,ss))
         !    do jj=1, M
         !    do kk=1, M
         !       aux2 = aux2 + dble(rho_aux(jj,kk,ss) * aux_mat5(kk,jj,ss))
         !    end do
         !    end do
         !    if (ss==1) write(777,*) aux2
         !    if (ss==2) write(888,*) aux2
         ! end if

!###############################################################################

      end do
         if (ceed_calc==1) then
            aux_mat1(:,:,1) = d2dip_ceed(ii,1) * aux_mat1(:,:,1)
            if (open_shell) aux_mat1(:,:,2) = d2dip_ceed(ii,2) * aux_mat1(:,:,2)
            ! aux1     = d2dip_ceed(ii,1)
            ! if (open_shell) aux1 = aux1 + d2dip_ceed(ii,2)
            ! aux_mat1 = aux1 * aux_mat1
         else if (ceed_calc==2) then
            do jj=1,M
            do kk=1,M
               aux1 = d2mu_vec(jj,1) + d2mu_vec(kk,1)
               if (open_shell) aux1 = aux1 + d2mu_vec(jj,2) + d2mu_vec(kk,2)
               aux_mat1(jj,kk, :) = aux1 * aux_mat1(jj,kk, :)
            end do
            end do
         end if
         aux_mat3 = aux_mat3 + k_ceed * A_ceed * aux_mat1
      end do

      do ss = 1,dim3
         call fock_ceed_op(ss)%Sets_data_ON(aux_mat3(:,:,ss))
         call fock_ceed_op(ss)%Commut_data_c(rho_aux(:,:,ss),                  &
                                             aux_mat4(:,:,ss), M)
      end do

   end if

   if (ceed_calc==2) call C_ON_mat_ceed%change_base(aux_mat4(:,:,1), 'inv')

   rho_aux = aux_mat4


end subroutine ceed_fock_calculation

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
subroutine ceed_finalize()
   use ceed_data, only: ceed_calc, dip_ceed_op, d2ip_ceed_op,                  &
                        fock_ceed_op, d2dip_ceed, Xmat_ceed, C_ON_mat_ceed,    &
                        YCinv_ceed
   implicit none

   if (ceed_calc==0) return
   if (allocated(dip_ceed_op) ) deallocate(dip_ceed_op)
   if (allocated(d2ip_ceed_op)) deallocate(d2ip_ceed_op)
   if (allocated(fock_ceed_op)) deallocate(fock_ceed_op)
   if (allocated(d2dip_ceed)  ) deallocate(d2dip_ceed)
   call Xmat_ceed%destroy()
   call YCinv_ceed%destroy()
   call C_ON_mat_ceed%destroy()

end subroutine ceed_finalize
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module ceed_subs
