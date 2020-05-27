#include "datatypes/datatypes.fh"
module DOS_subs
    implicit none
contains

subroutine init_PDOS(M)
    ! This subroutine initialize the variables for PDOS and DOS calculation, and
    ! read the PDOS_dat.in file.
    use DOS_data, only: pdos_calc, pdos_allb , pdos, pdos_nuc,       &
                        pdos_natoms, pdos_nbases, pdos_b, min_level, &
                        dos_nsteps, dos_sigma, dos_Eref, dos_calc,   &
                        trans_calc, trans_coef

    implicit none
    integer, intent(in) :: M

    if (.not.dos_calc) return

    open(unit = 10203, file = "PDOS_dat.in")
    read(10203,*) min_level
    read(10203,*) dos_nsteps
    read(10203,*) dos_sigma
    read(10203,*) dos_Eref

   if (pdos_calc) then

      allocate(pdos(M))

      if (.not. pdos_allb) then
         read(10203,*) pdos_natoms
         allocate(pdos_nuc(pdos_natoms))
         read(10203,*) pdos_nuc
      else
         read(10203,*) pdos_natoms, pdos_nbases
         allocate(pdos_nuc(pdos_natoms), pdos_b(M,pdos_nbases))
         read(10203,*) pdos_nuc
      endif
   endif

   if (trans_calc) allocate(trans_coef(M))

   close(10203)
end subroutine

subroutine build_PDOS(coef_mat, overlap, M, M_total, Nuc)
   ! This subroutine calculate the weights of each atom or basis, over the DOS.
   use DOS_data  , only: pdos_calc, pdos_allb, pdos_nuc, pdos_natoms, &
                         pdos_nbases, pdos, pdos_b
   use tbdft_data, only: MTB

   implicit none
   integer     , intent(in)  :: M, M_total
   integer     , intent(in)  :: Nuc(M)
   LIODBLE, intent(in)  :: coef_mat(M_total,M_total)
   LIODBLE, intent(in)  :: overlap(M,M)
   integer     , allocatable :: index_b(:)
   integer :: ii, jj, kk, ll

   if(.not.pdos_calc) return

   pdos   = 0.0d0
   pdos_b = 0.0d0

   allocate(index_b(pdos_nbases))
      do ll = 1, pdos_natoms
      do kk = 1, M
         if (pdos_nuc(ll) == Nuc(kk)) then
            do jj = 1, M_total
            do ii = 1, M
               pdos(jj) = pdos(jj) + coef_mat(ii+MTB,jj) * coef_mat(kk+MTB,jj) &
                          * overlap(ii,kk)
            enddo
            enddo
         endif
      enddo
      enddo

   if (pdos_allb) then

      do ll=1, pdos_natoms
         kk = 1
         do ii = 1, M
            if (pdos_nuc(ll) == Nuc(ii)) then
               index_b(kk) = ii
               kk = kk+1
            endif
         enddo

         do kk=1, pdos_nbases
         do jj = 1, M_total
         do ii = 1, M
            pdos_b(jj,kk) = pdos_b(jj,kk) + coef_mat(ii+MTB,jj) *            &
                            coef_mat(index_b(kk)+MTB,jj) * overlap(ii,index_b(kk))
         enddo
         enddo
         enddo
      end do
   end if

end subroutine build_PDOS

subroutine write_DOS (M_in, morb_energy)
   ! This subroutine write the DOS and PDOS in the corresponding outputs.
   use DOS_data, only: dos_calc, pdos_calc, pdos, dos_nsteps , dos_Eref,        &
                       dos_sigma, pdos_b, pdos_allb, pdos_nbases, min_level

   implicit none
   integer     , intent(in) :: M_in
   LIODBLE, intent(in) :: morb_energy(M_in)

   LIODBLE      :: min_e, max_e, delta
   LIODBLE      :: pexpf, expf, density
   LIODBLE      :: x0, xx
   integer           :: ii, jj, kk
   character(len=12) :: outfile

   if ( .not. (dos_calc .or. pdos_calc)) return

   pexpf = 1.0d0
   expf  = -1.0d0/(2.0d0*dos_sigma**2)

   min_e = morb_energy(min_level) - &
          (morb_energy(M_in) - morb_energy(min_level)) / (M_in +1 -min_level)
   max_e = morb_energy(M_in) - &
          (morb_energy(M_in) - morb_energy(min_level)) / (M_in +1 -min_level)
   delta = (max_e - min_e) / dble(dos_nsteps)

   open(unit = 10203, file = "DOS.out")
   x0 = min_e - dos_Eref
   xx = 0.0D0

   do ii = 1, dos_nsteps
      density = 0.0d0
      xx = x0 + delta * dble(ii)
      do jj = min_level, M_in
         density = density + pexpf * &
                   dexp(expf * (xx + dos_Eref - morb_energy(jj)) ** 2.0D0)
      enddo
      if (density > 1.0D-16) write(10203, *) xx, density
   enddo
   close(10203)

   if (pdos_calc) then
      open(unit = 10803, file = "PDOS_all.out")
      xx = 0.0D0

      do ii = 1, dos_nsteps
         density = 0.0d0
         xx = x0 + delta * dble(ii)
         do jj = min_level, M_in
            density = density + pexpf * pdos(jj) * &
                      dexp(expf*(xx + dos_Eref - morb_energy(jj)) ** 2.0D0)
         enddo

         if (density > 1.0D-6) write(10803, *) xx, density
      enddo

      close(10803)
   endif

   if (pdos_calc .and. pdos_allb) then
      do kk = 1, pdos_nbases
         if (kk < 10) then
            write(outfile, "(A7,I1,A4)") "PDOS_b0", kk, ".out"
         else
            write(outfile, "(A6,I2,A4)") "PDOS_b", kk, ".out"
         endif

         open(unit = 10803,file = trim(outfile))

         xx = 0.0d0

         do ii = 1, dos_nsteps
            density = 0.0d0
            xx = x0 + delta * dble(ii)
            do jj = min_level, M_in
               density = density + pexpf * pdos_b(jj,kk) * &
                         dexp(expf*(xx + dos_Eref - morb_energy(jj)) ** 2.0D0)
            enddo

            if (density > 1.0D-16) write(10803, *) xx, density
         enddo

         close(10803)
     enddo
   endif

endsubroutine write_DOS

subroutine write_trans_func(coef_mat, morb_energy, overlap, M, M_total,        &
                            fock_op)

   use DOS_data        , only: trans_coef, trans_calc, min_level, dos_nsteps, &
                               dos_sigma, dos_Eref
   use tbdft_data      , only: MTB, n_biasTB, tbdft_calc
   use typedef_operator, only: operator

   implicit none
   integer, intent(in)        :: M
   integer, intent(in)        :: M_total
   LIODBLE, intent(in)        :: coef_mat(M_total,M_total)
   LIODBLE, intent(in)        :: morb_energy(M_total)
   LIODBLE, intent(in)        :: overlap(M, M)
   type(operator), intent(in) :: fock_op

   integer  :: ii, jj, kk, ee
   integer  :: Nb
   LIODBLE  :: fockAO(M_total, M_total)
   LIODBLE  :: aux1
   LIODBLE  :: min_e, max_e, delta
   LIODBLE  :: pexpf, expf, T_E
   LIODBLE  :: x0, xx

   if(.not.trans_calc) return

   Nb = MTB/n_biasTB

   trans_coef = 0.0d0

   call fock_op%Gets_data_AO(fockAO)

   do ee = 1, M_total
      aux1 = 0.0d0
      do ii = 1, M
      do jj = 1, M
      do kk = 1, Nb
          aux1 = aux1 + fockAO(MTB+ii,kk)*coef_mat(kk,ee)*coef_mat(MTB+jj,ee)* &
                 overlap(jj,ii)
      enddo
      enddo
      enddo
      trans_coef(ee) = abs(aux1)
   enddo

   pexpf = 1.0d0
   expf  = -1.0d0/(2.0d0*dos_sigma**2)

   min_e = morb_energy(min_level) - &
          (morb_energy(M_total) - morb_energy(min_level))/(M_total+1-min_level)
   max_e = morb_energy(M_total) - &
          (morb_energy(M_total) - morb_energy(min_level))/(M_total+1 -min_level)
   delta = (max_e - min_e) / dble(dos_nsteps)

   open(unit = 10203, file = "Transmission.out")
   x0 = min_e - dos_Eref
   xx = 0.0D0

   do ii = 1, dos_nsteps
      T_E = 0.0d0
      xx = x0 + delta * dble(ii)
      do jj = min_level, M_total
         T_E = T_E + trans_coef(jj) * pexpf * &
                   dexp(expf * (xx + dos_Eref - morb_energy(jj)) ** 2.0D0)
      enddo
      if (T_E > 1.0D-16) write(10203, *) xx, T_E
   enddo
   close(10203)

endsubroutine write_trans_func

end module DOS_subs
