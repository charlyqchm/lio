module ediis_subs

   implicit none
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine EDIIS_init(M_in, OP_shell)

   use ediis_data, only: nediis, damp_ediis, fock_ediis_mat, rho_ediis_mat,&
                         BMAT, EDIIS_E, fock_damp

   implicit none
   logical, intent(in) :: OP_shell
   integer, intent(in) :: M_in

   if (OP_shell) then
      allocate(fock_ediis_mat(M_in,M_in,nediis,2),                             &
               rho_ediis_mat(M_in,M_in,nediis,2), BMAT(nediis,nediis,2),       &
               EDIIS_E(nediis), fock_damp(M_in,M_in,2))
   else
      allocate(fock_ediis_mat(M_in,M_in,nediis,1),                             &
               rho_ediis_mat(M_in,M_in,nediis,1), BMAT(nediis,nediis,1),       &
               EDIIS_E(nediis), fock_damp(M_in,M_in,1))
   end if

   BMAT           = 0.0d0
   fock_ediis_mat = 0.0d0
   rho_ediis_mat  = 0.0d0
   EDIIS_E        = 0.0d0
   fock_damp      = 0.0d0

end subroutine EDIIS_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef CUBLAS
subroutine ediis_conver (niter, M_in, energy, devPtrX, devPTrY, OPEN, fock_aop,&
                         rho_aop, fock_bop, rho_bop)
#else
subroutine ediis_conver (niter, M_in, energy ,Xmat, Ymat, OPEN, fock_aop,      &
                         rho_aop, fock_bop, rho_bop)
#endif

   use ediis_data,       only: nediis, damp_ediis, fock_ediis_mat,            &
                               rho_ediis_mat, BMAT, step_nediis, EDIIS_E,     &
                               fock_damp
   use typedef_operator, only: operator

   implicit none

   type(operator), intent(inout)           :: fock_aop, rho_aop
   type(operator), intent(inout), optional :: fock_bop, rho_bop

   logical, intent(in) :: OPEN
   integer, intent(in) :: niter
   integer, intent(in) :: M_in
   real*8, intent(in)  :: energy
   logical             :: acceptable
   integer             :: position
   integer             :: ii, jj, kk
   real*8              :: traza
   real*8              :: norm_factor
   real*8, allocatable :: mat_aux1(:,:,:), mat_aux2(:,:,:)
   real*8, allocatable :: fock_new(:,:,:)
   real*8, allocatable :: BMAT_aux(:,:,:)
   real*8, allocatable :: EDIIS_coef(:,:)

#ifdef  CUBLAS
   integer*8, intent(in) :: devPtrX
   integer*8, intent(in) :: devPtrY
#else
   real*8,  intent(in) :: Xmat(M_in,M_in)
   real*8,  intent(in) :: Ymat(M_in,M_in)
#endif

   step_nediis=niter

   position = min(step_nediis, nediis)

   if (allocated(mat_aux1)) deallocate(mat_aux1)
   if (allocated(mat_aux2)) deallocate(mat_aux2)
   if(allocated(EDIIS_coef)) deallocate(EDIIS_coef)
   if(allocated(fock_new)) deallocate(fock_new)


   if (OPEN) then
      allocate(mat_aux1(M_in,M_in,2),mat_aux2(M_in,M_in,2),                    &
               EDIIS_coef(position,2), fock_new(M_in,M_in,2))
   else
      allocate(mat_aux1(M_in,M_in,1),mat_aux2(M_in,M_in,1),                    &
               EDIIS_coef(position,1),fock_new(M_in,M_in,1) )
   end if


   mat_aux1   = 0.0d0
   mat_aux2   = 0.0d0
   EDIIS_coef = 0.0d0
   fock_new   = 0.0d0
   acceptable = .true.
!-------------------------------------------------------------------------------
!Updating Fock, rho and the energy:
!-------------------------------------------------------------------------------

   if (step_nediis>nediis) then
      do ii=1, nediis
         fock_ediis_mat(:,:,ii,1) = fock_ediis_mat(:,:,ii+1,1)
         rho_ediis_mat(:,:,ii,1) = rho_ediis_mat(:,:,ii+1,1)
         if(OPEN) then
            fock_ediis_mat(:,:,ii,2) = fock_ediis_mat(:,:,ii+1,2)
            rho_ediis_mat(:,:,ii,2) = rho_ediis_mat(:,:,ii+1,2)
         endif
         EDIIS_E(ii) = EDIIS_E(ii+1)
      end do
   end if

#ifdef CUBLAS
    call rho_aop%BChange_AOtoON(devPtrY, M_in, 'r')
    call fock_aop%BChange_AOtoON(devPtrX, M_in, 'r')
    if(OPEN)then
      call rho_bop%BChange_AOtoON(devPtrY, M_in, 'r')
      call fock_bop%BChange_AOtoON(devPtrX, M_in, 'r')
    end if
#else
    call rho_aop%BChange_AOtoON(Ymat, M_in, 'r')
    call fock_aop%BChange_AOtoON(Xmat,M_in, 'r')
    if (OPEN) then
       call rho_bop%BChange_AOtoON(Ymat, M_in, 'r')
       call fock_bop%BChange_AOtoON(Xmat,M_in, 'r')
    endif
#endif

   call fock_aop%Gets_data_ON(mat_aux1(:,:,1))
   call rho_aop%Gets_data_ON(mat_aux2(:,:,1))

   fock_ediis_mat(:,:,position,1) = mat_aux1(:,:,1)
   rho_ediis_mat(:,:,position,1) = mat_aux2(:,:,1)

   if (OPEN) then
      call fock_bop%Gets_data_ON(mat_aux1(:,:,2))
      call rho_bop%Gets_data_ON(mat_aux2(:,:,2))

      fock_ediis_mat(:,:,position,1) = mat_aux1(:,:,2)
      rho_ediis_mat(:,:,position,1) = mat_aux2(:,:,2)
   end if


   EDIIS_E(position) = energy

!-------------------------------------------------------------------------------
!First and second steps
!-------------------------------------------------------------------------------
   if (step_nediis==1) then
      fock_new  = mat_aux1
      fock_damp = fock_new
   else if (step_nediis==2) then !>=2.and.step_nediis<=150) then
      fock_new = (mat_aux1 + 80.0d0*fock_damp)/(1.0d0+80.0d0)
      fock_damp=fock_new
   end if

!-------------------------------------------------------------------------------
!Updating BMAT
!-------------------------------------------------------------------------------

   if (allocated(BMAT_aux)) deallocate(BMAT_aux)
   if(OPEN) then
      allocate(BMAT_aux(position,position,2))
   else
      allocate(BMAT_aux(position,position,1))
   endif

   BMAT_aux = 0.0d0

   if (step_nediis>1.and.step_nediis<=nediis) then

      do jj = 1, position-1
      do ii = 1, position-1
         BMAT_aux(ii,jj,1) = BMAT(ii,jj,1)
         if(OPEN)  BMAT_aux(ii,jj,2) = BMAT(ii,jj,2)
      end do
      end do

   else if (step_nediis>nediis) then

      do jj = 1, position-1
      do ii = 1, position-1
         BMAT_aux(ii,jj,1) = BMAT(ii+1,jj+1,1)
         if(OPEN)  BMAT_aux(ii,jj,2) = BMAT(ii+1,jj+1,2)
      end do
      end do

   end if

   mat_aux1 = 0.0d0
   mat_aux2 = 0.0d0

   do jj=1, position-1
      mat_aux1(:,:,1) = fock_ediis_mat(:,:,jj,1)-fock_ediis_mat(:,:,position,1)
      mat_aux2(:,:,1)= rho_ediis_mat(:,:,jj,1)-rho_ediis_mat(:,:,position,1)

      call matmul_traza(mat_aux1(:,:,1),mat_aux2(:,:,1), M_in, traza)
      BMAT_aux(jj,position,1) = traza

      if (OPEN) then
         mat_aux1(:,:,2) = fock_ediis_mat(:,:,jj,2)-fock_ediis_mat(:,:,position,2)
         mat_aux2(:,:,2)= rho_ediis_mat(:,:,jj,2)-rho_ediis_mat(:,:,position, 2)
         call matmul_traza(mat_aux1(:,:,2),mat_aux2(:,:,2), M_in, traza)
         BMAT_aux(jj,position,2) = traza
      end if

      BMAT_aux(position, jj,:) = BMAT_aux(jj,position,:)

   end do


   BMAT(1:position,1:position,:) = BMAT_aux(1:position,1:position,:)

   do ii=1, position
      EDIIS_coef(ii,:) =  EDIIS_E(ii)
   end do

!-------------------------------------------------------------------------------
!Solving linear equation and getting new fock:
!-------------------------------------------------------------------------------

   if (step_nediis>2) then

!charly
      do ii=1, position
         write(222, '(30E14.6)') BMAT_aux(ii,:,1)
      end do

      write(222,*) "-----------------------"
      write(222,*) (EDIIS_E(ii), ii=1,position)

      call solve_linear_constraints(EDIIS_coef(:,1), EDIIS_E(1:position),       &
                                    BMAT_aux(:,:,1), position)
      if (OPEN) call solve_linear_constraints(EDIIS_coef(:,2),                  &
                                              EDIIS_E(1:position),              &
                                              BMAT_aux(:,:,2), position)

      do ii=1, position
         fock_new(:,:,1) = fock_new(:,:,1) +                                    &
                           EDIIS_coef(ii,1)*fock_ediis_mat(:,:,ii,1)
         if (OPEN) fock_new(:,:,2) = fock_new(:,:,2) +                          &
                                     EDIIS_coef(ii,2)*fock_ediis_mat(:,:,ii,2)
      end do

          write(222,*) "-----------------------"
          write(222, '(30E14.6)') EDIIS_coef(:,1)
          write(222,*) "#################################"
   end if


   call fock_aop%Sets_data_ON(fock_new(:,:,1))
   if (OPEN) call fock_bop%Sets_data_ON(fock_new(:,:,2))

end subroutine ediis_conver

subroutine matmul_traza(mat1,mat2, M_in, traza)

   implicit none
   integer, intent(in) :: M_in
   real*8, intent(out) :: traza
   real*8, intent(in)  :: mat1(M_in,M_in), mat2(M_in, M_in)
   integer             :: ii, jj
   real*8              :: mat3(M_in)

   mat3=0.0d0
   traza=0.0d0

   do ii=1, M_in
   do jj=1, M_in
      mat3(ii) = mat1(ii,jj)*mat2(jj,ii) + mat3(ii)
   end do
   end do

   do ii=1, M_in
      traza=traza+mat3(ii)
   end do
end subroutine matmul_traza

subroutine solve_linear_constraints(coef, Ener, BMAT, ndim)

   implicit none
   integer, intent(in) :: ndim
   real*8, intent(out) :: coef(ndim)
   real*8, intent(in)  :: Ener(ndim)
   real*8, intent(in)  :: BMAT(ndim,ndim)
   logical   :: converged
   logical   :: big_alpha1, big_alpha2
   logical   :: update
   integer   :: ii, jj, conv_steps
   integer   :: yind
   integer   :: zind(ndim-1)
   integer   :: lastindex, newindex
   real*8    :: aux_coef(ndim), new_coef(ndim)
   real*8    :: grad(ndim)
   real*8    :: delta(ndim)
   real*8    :: r_grad(ndim-1)
   real*8    :: alpha1, alpha2, alpha3, alpha_aux
   real*8    :: vec_alpha2(ndim-1)
   real*8    :: result1, result2, result3

   coef         = 1.0d0/dble(ndim)
   new_coef     = 0.0d0
   delta        = 0.0d0
   aux_coef     = 0.0d0
   r_grad       = 0.0d0
   lastindex    = 0
   newindex     = 1
   yind         = 1
   alpha1       = 0.0d0
   alpha2       = 0.0d0
   alpha3       = 0.0d0
   alpha_aux    = 0.0d0
   converged    = .true.
   conv_steps   = 0

   do ii=2, ndim
      zind(ii-1) = ii
   end do

   do while (converged.and.conv_steps<=100000)
      conv_steps=conv_steps+1
!charly
      if (conv_steps>100000) then
         write(222,*) "Too many steps"
      end if

      big_alpha1   = .false.
      big_alpha2   = .false.
      update       = .false.
      converged    = .false.

      call gradient(coef, grad, Ener, BMAT ,ndim)
      call displacement (grad, zind, yind, delta, coef, ndim)

      do ii=1, ndim-1
         if (abs(delta(zind(ii)))>1.0D-8) then
            converged=.true.
            exit
         end if
      end do

      if (converged.eqv..false.) then
        exit
      end if

      if (delta(yind) < 0.0d0) then
         alpha1=-coef(yind)/delta(yind)
      else
         big_alpha1 = .true.
      end if

      do ii=1, ndim-1
         vec_alpha2(ii) = -delta(zind(ii))/coef(zind(ii))
      end do

      alpha2=maxval(vec_alpha2)
      alpha2=1.0d0/alpha2

      if (alpha2<=0.0d0) then
         big_alpha2=.true.
      end if

      call min_alpha(Ener,BMAT, coef,delta, alpha3, ndim)
!charly:
       ! write(666,*) big_alpha1, big_alpha2
       ! write(666,'(3E14.6)') alpha1, alpha2, alpha3
       ! write(666,'(10E14.6)') coef
      ! write(666,'(10E14.6)') delta


      if (big_alpha1.and.big_alpha2) then

         call f_coef(Ener,BMAT, coef+alpha3*delta, result1, ndim)
         call f_coef(Ener,BMAT, coef, result2, ndim)
         ! if (result2<result1) then
         !    print*, "convergence problem1"
         !    stop
         ! end if

      else if (big_alpha1) then

         if (alpha3>alpha2) alpha3=alpha2
         call f_coef(Ener,BMAT, coef+alpha2*delta, result1, ndim)
         call f_coef(Ener,BMAT, coef+alpha3*delta, result2, ndim)
         call f_coef(Ener,BMAT, coef, result3, ndim)
         if (result1<result2) alpha3=alpha2
         ! if (result3<result1.and.result3<result2) then
         !    print*, "convergence problem2"
         !    stop
         ! end if

      else if (big_alpha2) then

         if (alpha3>alpha1) then
            alpha3=alpha1
            update=.true.
         end if
         call f_coef(Ener,BMAT, coef+alpha1*delta, result1, ndim)
         call f_coef(Ener,BMAT, coef+alpha3*delta, result2, ndim)
         call f_coef(Ener,BMAT, coef, result3, ndim)
         if (result1<result2) then
            alpha3=alpha1
            update=.true.
         end if
         ! if (result3<result1.and.result3<result2) then
         !    print*, "convergence problem3"
         !    stop
         ! end if

      else

         if(alpha1<alpha2) then
            alpha_aux=alpha1
         else
            alpha_aux=alpha2
         end if
         if (alpha3>alpha_aux) alpha3=alpha_aux
         call f_coef(Ener,BMAT, coef+alpha_aux*delta, result1, ndim)
         call f_coef(Ener,BMAT, coef+alpha3*delta, result2, ndim)
         call f_coef(Ener,BMAT, coef, result3, ndim)
         if (result1<result2) alpha3=alpha_aux
         ! if (result3<result1.and.result3<result2) then
         !    print*, "convergence problem4"
         !    stop
         ! end if
         if(alpha3>=alpha1) update=.true.
      end if

      new_coef = coef+alpha3*delta
      coef     = new_coef

      if (update) then
         newindex       = newindex+1
         if (newindex==ndim+1) newindex = 2
         lastindex      = yind
         yind           = zind(newindex-1)
         zind(newindex-1) = lastindex
      end if

!chaly: checking the restriction
       result1=0.0d0
       do ii=1, ndim
          result1=coef(ii)+result1
       end do

       ! write(666,*) "result", result1

       if (abs(result1-1.0d0)>1.0d-8) then
          print*,"The restricion is not comply"
          stop
       end if

!charly
!       write(666,'(5F10.6)') coef
!      write(555,'(A4,5E14.4)') "grad", grad
!      write(444,'(A5,5F10.6)') "delta", delta
   end do
!charly
       ! write(666,*) "------------------"
!      write(555,*) "------------------"
!      write(444,*) "------------------"
   return

end subroutine solve_linear_constraints

subroutine gradient(coef, grad, Ener, BMAT ,ndim)

   implicit none
   integer, intent(in) :: ndim
   real*8, intent(in)  :: coef(ndim)
   real*8, intent(in)  :: Ener(ndim)
   real*8, intent(in)  :: BMAT(ndim,ndim)
   real*8, intent(out) :: grad(ndim)
   integer   :: ii, jj

   grad=0.0d0

   do ii=1, ndim
   do jj=1, ndim
      grad(ii) = -BMAT(ii,jj)*coef(jj)+grad(ii)
   end do
      grad(ii) = grad(ii) + Ener(ii)
   end do

end subroutine gradient

subroutine displacement(grad, zind, yind, delta, coef, ndim)

   implicit none
   integer, intent(in) :: ndim
   real*8, intent(in)  :: coef(ndim)
   real*8, intent(in)  :: grad(ndim)
   real*8, intent(out) :: delta(ndim)
   integer, intent(in) :: yind
   integer, intent(in) :: zind(ndim-1)
   integer   :: ii
   real*8    :: r_grad(ndim-1)

   delta=0.0d0

   do ii=1, ndim-1
      r_grad(ii) = grad(zind(ii)) - grad(yind)
      if (r_grad(ii)<0.0d0.or.coef(zind(ii))>0.0d0) then
         delta(zind(ii)) = -r_grad(ii)
      else
         delta(zind(ii)) = 0.0d0
      end if
   end do

   do ii=1, ndim-1
       delta(yind) = - delta(zind(ii)) + delta(yind)
   end do
end subroutine displacement

subroutine min_alpha(Ener,BMAT, coef,delta, alpha, ndim)

   implicit none
   integer, intent(in)  :: ndim
   real*8,  intent(in)  :: Ener(ndim)
   real*8,  intent(in)  :: BMAT(ndim,ndim)
   real*8,  intent(in)  :: coef(ndim)
   real*8,  intent(in)  :: delta(ndim)
   real*8,  intent(out) :: alpha
   integer    :: ii, jj
   real*8     :: num1, num2, den1

   num1 = 0.0d0
   num2 = 0.0d0
   den1 = 0.0d0
   alpha= 0.0d0

   do ii=1, ndim
      num1=num1+Ener(ii)*delta(ii)
   end do

   do jj=1, ndim
   do ii=1, ndim
      num2=num2+BMAT(ii,jj)*delta(ii)*coef(jj)
      den1=den1+BMAT(ii,jj)*delta(ii)*delta(jj)
   end do
   end do

   alpha = (num1-num2)/den1

end subroutine min_alpha

subroutine f_coef(Ener,BMAT, coef, result, ndim)

   implicit none
   integer, intent(in)  :: ndim
   real*8,  intent(in)  :: Ener(ndim)
   real*8,  intent(in)  :: BMAT(ndim,ndim)
   real*8,  intent(in)  :: coef(ndim)
   real*8,  intent(out) :: result
   integer              :: ii,jj
   real*8    :: sum1, sum2

   sum1   = 0.0d0
   sum2   = 0.0d0
   result = 0.0d0

   do ii=1, ndim
      sum1=sum1+Ener(ii)*coef(ii)
   end do

   do jj=1, ndim
   do ii=1, ndim
      sum2=sum2+BMAT(ii,jj)*coef(ii)*coef(jj)
   end do
   end do

   result=sum1-0.5d0*sum2

end subroutine

end module ediis_subs
