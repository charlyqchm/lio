#include "datatypes/datatypes.fh"
module ceed_subs
   implicit none
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ceed_saving_coef(coef, M)
   use ceed_data, only: ceed_calc, ceed_coef, ceed_coefT
   implicit none
   integer, intent(in) :: M
   LIODBLE, intent(in) :: coef(M,M)

   if (.not.ceed_calc) return

   allocate(ceed_coef(M,M), ceed_coefT(M,M))
   ceed_coef  = coef

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ceed_init(M, open_shell, r, d, natom, ntatom, propagator)
! This subroutine initialize the variables for CEED calculations
   use ceed_data, only: ceed_calc, dip_ceed_op, d2ip_ceed_op,                  &
                         fock_ceed_op, d2dip_ceed, Xmat_ceed, ceed_coef,       &
                         ceed_coefT, Ytrans_ceed
   use faint_cpu , only: intfld

   implicit none
   logical, intent(in)  :: open_shell
   integer, intent(in)  :: M
   integer, intent(in)  :: propagator
   integer, intent(in)  :: natom
   integer, intent(in)  :: ntatom
   LIODBLE, intent(in)  :: r(ntatom,3)
   LIODBLE, intent(in)  :: d(natom,natom)
   integer              :: ii, jj
   integer              :: MM
   LIODBLE              :: vec_aux(3)
   LIODBLE, allocatable :: dip_array(:)
   LIODBLE, allocatable :: dip_mat_aux(:,:)
   LIODBLE              :: aux_mat(M,M)

   if (.not. ceed_calc) return

   if (propagator /= 1) then
      write(*,*) "CEED calculation can only be performed with Verlet propagator"
      write(*,*) "Stopping LIO"
      stop
   end if

   call Ytrans_ceed%multiply(aux_mat, ceed_coef)
   ceed_coef = aux_mat
   ceed_coefT= transpose(ceed_coef)

   MM = M*(M+1)/2

   allocate(dip_ceed_op(3), dip_array(MM), dip_mat_aux(M, M))

   if (.not.open_shell) then
      allocate(fock_ceed_op(1), d2ip_ceed_op(3,1), d2dip_ceed(3,1))
   else
      allocate(fock_ceed_op(2), d2ip_ceed_op(3,2), d2dip_ceed(3,2))
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
!changing to OM base:
      call dip_ceed_op(ii)%Gets_data_ON(dip_mat_aux)
      dip_mat_aux = matmul(dip_mat_aux, ceed_coef)
      dip_mat_aux = matmul(ceed_coefT , dip_mat_aux)
      call dip_ceed_op(ii)%Sets_data_ON(dip_mat_aux)

   end do

   deallocate(dip_array, dip_mat_aux)
end subroutine ceed_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ceed_fock_calculation(fock_in, rho_in, M, t_step, dim3, open_shell)
!This subroutine includes the CEED term in the evolution of the density matrix.
!At the end of the subroutine rho_aux store the commutator -i[H_CEED,\rho],
!were:
!              H_CEED =  -i A_ceed d^2<\mu>/dt^2 . [\mu , H]
!d^2<\mu>/dt^2 is approximated as:
!               d^2<\mu>/dt^2 = - Tr(\rho . [[\mu,H],H])
!and:
!             A_ceed = 2/3 \alpha/c^2 = 2.592668788247536d-7

   use ceed_data, only: ceed_calc, dip_ceed_op, d2ip_ceed_op, d2dip_ceed,      &
                        A_ceed, fock_ceed_op, k_ceed, ceed_td_step, ceed_coef, &
                        ceed_coefT
   implicit none

   logical   , intent(in)    :: open_shell
   integer   , intent(in)    :: M
   integer   , intent(in)    :: dim3
   integer   , intent(in)    :: t_step
   LIODBLE   , intent(in)    :: fock_in(M,M,dim3)
   TDCOMPLEX , intent(inout) :: rho_in(M,M,dim3)
   integer                :: ii, jj, kk, ss
   LIODBLE                :: aux1, aux2
   LIODBLE                :: fock(M,M,dim3)
   TDCOMPLEX              :: rho_aux(M,M, dim3)
   LIODBLE                :: aux_mat1(M,M,dim3)
   LIODBLE                :: aux_mat2(M,M,dim3)
   LIODBLE                :: aux_mat3(M,M,dim3)
   LIODBLE                :: dmu2_vec(M)
   TDCOMPLEX              :: aux_mat4(M,M,dim3)
   TDCOMPLEX              :: liocmplx

   if (.not.ceed_calc) return

   aux_mat3 = 0.0d0
   aux_mat4 = liocmplx(0.0d0,0.0d0)
   fock(:,:,1)   = matmul(fock_in(:,:,1),ceed_coef)
   fock(:,:,1)    = matmul(ceed_coefT,fock(:,:,1))
   rho_aux(:,:,1)  = matmul(rho_in(:,:,1),ceed_coef)
   rho_aux(:,:,1)  = matmul(ceed_coefT,rho_aux(:,:,1))

   if (t_step > ceed_td_step) then
      do ii = 1,3
      do ss = 1,dim3
         call dip_ceed_op(ii)%Commut_data_r(fock(:,:,ss), aux_mat1(:,:,ss), M)
         call d2ip_ceed_op(ii,ss)%Sets_data_ON(aux_mat1(:,:,ss))
         call d2ip_ceed_op(ii,ss)%Commut_data_r(fock(:,:,ss),aux_mat2(:,:,ss),M)
         d2dip_ceed(ii,ss) = 0.0d0
         dmu2_vec = 0.0d0
         do jj=1, M
         do kk=1, M
            aux1 = dble(rho_aux(jj,kk,ss)*aux_mat2(kk,jj,ss))
            d2dip_ceed(ii,ss) = d2dip_ceed(ii,ss) - aux1
            dmu2_vec(jj) = dmu2_vec(jj) - aux1
         end do
!            if(ii == 1 )write(777,*) " ",jj," ", dmu2_vec(jj)
         end do
      end do
         aux1  = d2dip_ceed(ii,1)
!Modificando cada pedaso del conmutador
         do jj =1, M
         do kk =1, M
            ! if (jj/=2 .and. kk/=2) then
               aux2 = (dmu2_vec(jj)+dmu2_vec(kk))*aux_mat1(jj,kk,1)
               aux_mat1(jj,kk,1) =  aux2
            ! else
               ! aux_mat1(jj,kk,ss) = 0.0d0
            ! end if
         end do
         end do
         ! aux_mat1 = aux1 * aux_mat1

         if (open_shell) aux1 = aux1 + d2dip_ceed(ii,2)
         aux_mat3 = aux_mat3 + k_ceed * A_ceed * aux_mat1
      end do

      do ss = 1,dim3
         call fock_ceed_op(ss)%Sets_data_ON(aux_mat3(:,:,ss))
         call fock_ceed_op(ss)%Commut_data_c(rho_aux(:,:,ss),                  &
                                             aux_mat4(:,:,ss), M)
      end do

   end if

   rho_aux(:,:,1) = matmul(aux_mat4(:,:,1),ceed_coefT)
   rho_aux(:,:,1) = matmul(ceed_coef, rho_aux(:,:,1))

   rho_in = rho_aux
end subroutine ceed_fock_calculation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ceed_finalize()
   use ceed_data, only: ceed_calc, dip_ceed_op, d2ip_ceed_op,            &
                        fock_ceed_op, d2dip_ceed, Xmat_ceed
   implicit none

   if (.not.ceed_calc) return
   if (allocated(dip_ceed_op) ) deallocate(dip_ceed_op)
   if (allocated(d2ip_ceed_op)) deallocate(d2ip_ceed_op)
   if (allocated(fock_ceed_op)) deallocate(fock_ceed_op)
   if (allocated(d2dip_ceed)  ) deallocate(d2dip_ceed)
   call Xmat_ceed%destroy()

end subroutine ceed_finalize
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module ceed_subs
