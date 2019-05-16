module exchange_data

   implicit none
   logical          :: PBE0_calc = .false.
   double precision :: a_PBE0    = 1.0d0
   double precision, allocatable :: Ccont(:,:)
   double precision, allocatable :: K_int(:,:,:,:)
   double precision, allocatable :: K_mat(:,:)

end module exchange_data

module exchange_subs

contains

subroutine exchange_init()

   use exchange_data, only: PBE0_calc, a_PBE0,Ccont, K_int, K_mat
   use basis_data   , only: c_raw, max_c_per_atom, M

   implicit none

   if (.not.PBE0_calc) return

   a_PBE0 = 0.75d0

   if (.not.allocated(Ccont)) allocate(Ccont(M, max_c_per_atom))
   if (.not.allocated(K_mat)) allocate(K_mat(M,M))

   Ccont=c_raw

end subroutine exchange_init

subroutine exchange_mat(rho, fock, Ex_exact)

   use exchange_data, only: K_mat,Ccont, a_PBE0, PBE0_calc
   use basis_data   , only: M

   implicit none
   double precision, intent(in)    :: rho(M,M)
   double precision, intent(inout) :: fock(M,M)
   double precision, intent(in)    :: wf
   double precision, intent(out), optional :: Ex_exact
   integer :: ii, jj

   if (.not.PBE0_calc) return

   K_mat=0.0d0
!carlos: no estas seguro aun que el segundo 0 dependa del paso o
!        si es una variable que indica inicializacion.

   call g2g_ex_exact(rho,Ccont, K_mat, 1, 0)

!carlos: esto no necesariamente sea util.
   if (present(Ex_exact)) then
      call Ex_calc(rho, K_mat, Ex_exact, M)
      Ex_exact = -(1.0d0-a_PBE0)*0.5d0*Ex_exact
   else
      fock = fock - (1.0d0-a_PBE0)*K_mat
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
