module exchange_data

   implicit none
   double precision, allocatable :: Ccont(:,:)
   double precision, allocatable :: K_int(:,:,:,:)
   double precision, allocatable :: K_mat(:,:)

end module exchange_data

module exchange_subs

contains

subroutine exchange_init()

   use exchange_data, only: Ccont, K_int, K_mat
   use basis_data   , only: c_raw, max_c_per_atom, M

   implicit none
   integer :: ii, jj, kk, ll

   if (.not.allocated(K_int)) allocate(K_int(M,M,M,M))
   if (.not.allocated(Ccont)) allocate(Ccont(M, max_c_per_atom))
   if (.not.allocated(K_mat)) allocate(K_mat(M,M))

   Ccont=c_raw
   K_int= 0.0d0

   call g2g_ex_int(Ccont, K_int)

end subroutine exchange_init

subroutine exchange_mat(rho, fock, wf, Ex_exact)

   use exchange_data, only: K_int, K_mat
   use basis_data   , only: M

   implicit none
   double precision, intent(in)    :: rho(M,M)
   double precision, intent(inout) :: fock(M,M)
   double precision, intent(in)    :: wf
   double precision, intent(out), optional :: Ex_exact
   integer :: ii, jj

   K_mat=0.0d0

   call g2g_exchange_mat(rho, K_mat, K_int)

   if (present(Ex_exact)) then
      call Ex_calc(rho, K_mat, Ex_exact, M)
      Ex_exact = -wf*0.5d0*Ex_exact
   else
      fock = fock - wf*K_mat
   end if

end subroutine exchange_mat

subroutine Ex_calc(rho, K_mat, Ex_exact, M)

   implicit none
   integer         , intent(in)  :: M
   double precision, intent(in)  :: rho(M,M), K_mat(M,M)
   double precision, intent(out) :: Ex_exact
   integer :: ii, jj

   Ex_exact = 0.0d0

   do ii=1, M
   do jj=1, M
      Ex_exact = Ex_exact + rho(ii,jj)*K_mat(jj,ii)
   end do
   end do

end subroutine Ex_calc

end module exchange_subs
