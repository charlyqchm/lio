module exchange_data

   implicit none
   double precision, allocatable :: Ccont(:,:)

end module exchange_data

module exchange_subs

contains

subroutine exchange_exact(rho, Kmat)

   use exchange_data, only: Ccont
   use basis_data   , only: c_raw, max_c_per_atom, M

   implicit none

   double precision, intent(in)  :: rho(M,M)
   double precision, intent(inout) :: Kmat(M,M)
   double precision              :: K_int(M,M,M,M)

   if (.not. allocated(Ccont)) allocate(Ccont(M, max_c_per_atom))

   write(333,*) c_raw
   Ccont=c_raw
   K_int= 0.0d0
   Kmat = 0.0d0

   call g2g_ex_exact(rho, Ccont, Kmat, K_int)

end subroutine exchange_exact

end module exchange_subs
