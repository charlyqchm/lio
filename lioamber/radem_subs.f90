#include "complex_type.fh"
module radem_subs
   implicit none
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine radi_emission_init(M_in, M, open_shell, r, d, natom, ntatom)
   use radem_data, only: radiative_calc, dip_radi_op, ddip_radi_op,            &
                         fock_radi_op, d2dip_radi, coef_r, coefT_r, Xmat_radi, &
                         dip_vec, npoints
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

   allocate(dip_radi_op(3),dip_array(MM), dip_mat_aux(M, M), dip_vec(npoints))

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
      call dip_radi_op(ii)%BChange_AOtoON(Xmat_radi, M)
   end do

end subroutine radi_emission_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine radi_fock_calculation(fock, rho_aux, M, dim3, istep)
   use radem_data, only: radiative_calc, dip_radi_op, ddip_radi_op,            &
                         d2dip_radi, preA_radi, fock_radi_op, Ymat_radi,       &
                         Xmat_radi, coef_r, coefT_r, dip_vec, npoints
   implicit none

   integer     , intent(in)    :: M
   integer     , intent(in)    :: dim3
   integer     , intent(in)    :: istep
   real(kind=8), intent(in)    :: fock(M,M)
   TDCOMPLEX   , intent(inout) :: rho_aux(M,M)
   integer                     :: ii, jj, kk
   real(kind=8)                :: aux1,aux2,aux3
   real(kind=8)                :: aux_mat1(M,M)
   real(kind=8)                :: aux_mat2(M,M)
   real(kind=8)                :: aux_mat3(M,M)
   real(kind=8)                :: aux_mat4(M,M)
   TDCOMPLEX                   :: aux_mat5(M,M)
   TDCOMPLEX                   :: i_cmplx=(0.0,1.0)
   real(kind=8)                :: x_dat(npoints)
   real(kind=8)                :: predict

   if (.not.radiative_calc) return
!   do ii=1, npoints
!      x_dat(ii) = dble(ii)
!   end do
   aux_mat4 = 0.0d0
   aux1 = 0.00d0
   aux2 = 0.0d0
   aux3 = 0.0d0

  if (istep > 100) then
    do ii=1,3
      call dip_radi_op(ii)%Gets_data_ON(aux_mat3)
      call dip_radi_op(ii)%Commut_data_r(fock, aux_mat1, M)
      call ddip_radi_op(ii,dim3)%Sets_data_ON(aux_mat1)
      call ddip_radi_op(ii,dim3)%Commut_data_r(fock, aux_mat2, M)
      d2dip_radi(ii,dim3) = 0.0d0
      do jj=1, M
      do kk=1, M
!         aux1 = aux1 + dble(rho_aux(jj,kk)*aux_mat3(kk,jj))
         aux2 = aux2 + dble(rho_aux(jj,kk)*i_cmplx*aux_mat1(kk,jj))
         aux3 = dble(rho_aux(jj,kk)*aux_mat2(kk,jj))
         d2dip_radi(ii,dim3) = d2dip_radi(ii,dim3) - aux3
      end do
      end do

      aux2=0.0d0
      if (ii==1) then
         do jj=1, M
            aux2 = aux2 + real(rho_aux(jj,jj))
         do kk=1, M
            aux1 = aux1 + dble(rho_aux(jj,kk)*aux_mat3(kk,jj))
         end do            
         end do
         write(666,*) aux1
         write(777,*) aux2
         write(888,*) d2dip_radi(ii,dim3)
      end if
!      write(777,*) aux2
!      write(888,*) d2dip_radi(ii,dim3)

       aux_mat4 = aux_mat4 - 1.0d5*preA_radi * d2dip_radi(ii,dim3) * aux_mat1
    end do
  endif

    call fock_radi_op(dim3)%Sets_data_ON(aux_mat4)
!    call fock_radi_op(dim3)%BChange_AOtoON(Xmat_radi, M)
!    call Ymat_radi%change_base(rho_aux, 'dir')
    call fock_radi_op(dim3)%Commut_data_c(rho_aux, aux_mat5, M)

    rho_aux = aux_mat5

end subroutine radi_fock_calculation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
SUBROUTINE FITTING(maxdeg,npoints, x_data,y_data, y, rms)
! Least-squares fit of a univariate data by a polynomial.
!
   implicit none
   integer :: iounit, info, ipoint, npoints, maxdeg, ideg, jdeg, lwork, nsv
   double precision x, y, x0, rms, dev, svtol
   double precision x_data(npoints), y_data(npoints)
   double precision, allocatable :: matrix(:,:), &
                                    vector(:), work(:), sv(:)
!   character(80) filename

   x0=0.0d0 !expansion point


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   form linear least-squares system   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  write(*,'(/1x,a,1x,i3,1x,a,1x,i3,1x,a)') 'compose', maxdeg+1, 'by', maxdeg+1, 'linear system'
! allocate arrays to keep linear system matrix and right-hand-side vector
   allocate(matrix(0:maxdeg,0:maxdeg), vector(0:maxdeg), stat=info)
    if (info/=0) stop 'stop: failed to allocate matrix and vector'
    matrix = 0.0
    vector = 0.0

! ---- start the most computationally expensive part ---- !
!          need to be optimized and parallelized
! compute linear system matrix and right-hand-side vector
  do ideg=0, maxdeg
    do jdeg=0, ideg
   matrix(jdeg,ideg) = sum( (x_data(1:npoints)-x0)**ideg * (x_data(1:npoints)-x0)**jdeg )
   if (jdeg/=ideg) matrix(ideg,jdeg) = matrix(jdeg,ideg)
   enddo
   enddo
   do ideg=0, maxdeg
     vector(ideg) = sum( (x_data(1:npoints)-x0)**ideg * y_data(1:npoints) )
   enddo
! --------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   solve linear equations system using LAPACK library   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,'(/a)') 'solve linear system by singular value decomposition'
! solve system by singular value decomposition
  lwork = (maxdeg+1)**2*32
  svtol = -1.0d-12
  allocate(work(lwork), sv(0:maxdeg), stat=info)
  if (info/=0) stop 'stop: failed to allocate work and sv arrays'
! use linear equations solver from Lapack/MKL (documentation https://software.intel.com/en-us/node/521114)
  call dgelss( maxdeg+1, maxdeg+1, 1, matrix(0:maxdeg,0:maxdeg), &
         maxdeg+1, vector(0:maxdeg), maxdeg+1, sv(0:maxdeg), svtol, nsv, work, lwork, info )
  if (info/=0) stop 'error: SVD failed'
! compute covariance (degree of correlation between two linear solutions
! or measure of linear dependence between two solutions)
! COVAR = V^T * S^{-1} * V,
!  where "V" is right singular matrix,
!   and "S" is diagonal matrix with singular values on diagonal
! allocate arrays to store covariance matrix and some intermediates
!allocate(covar(0:maxdeg,0:maxdeg), w()..., ... stat=info)
!if (info/=0) ...
! 1. compute w=S^{-1}
!w = 0.0
!forall (i=0:maxdeg) ...
! 2. compute S^{-1} * V
!call dgemm(...)
! 3. compute V^T * S^{-1} * V
!call dgemm(...)
! print values of the polynomial coefficients and diagonal elements of covariance matrix

! compute and print deviations between model polynomial and data points
   y = sum( (/( vector(ideg) * (x_data(npoints)+1.0)**ideg, ideg=0, maxdeg )/) )

   rms = 0.0d0
   do  ideg = 2, maxdeg
    rms = rms+vector(ideg)*(ideg-1)*(ideg)*((x_data(npoints)+1.0)**(ideg-2))
   enddo
! deallocate all allocatable arrays
    deallocate(matrix, vector,  work, sv)
   return
END SUBROUTINE FITTING

end module radem_subs
