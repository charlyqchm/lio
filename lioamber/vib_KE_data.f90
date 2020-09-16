#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module vib_KE_data
   use typedef_cumat   , only: cumat_r
   use typedef_cumat   , only: cumat_x

   implicit none

   logical              :: vib_calc = .false.
   integer              :: ke_calc  = 0
   integer, allocatable :: move_atom(:)
   integer              :: nat_move
   integer              :: hess_norder  = 1
   integer              :: n_vib
   integer              :: ke_start_t   = 100
   LIODBLE              :: delta_h      = 0.01d0 !Amstrongs
   LIODBLE              :: ke_sigma     = 0.1d0
   LIODBLE              :: phon_temp    = 0.00001d0!0.0d0
   integer, allocatable :: ke_ind(:)
   LIODBLE, allocatable :: atom_mass(:)
   LIODBLE, allocatable :: mass_w(:)
   LIODBLE, allocatable :: armonic_freq(:)
   LIODBLE, allocatable :: armonic_vec(:,:)
   LIODBLE, allocatable :: ke_coef(:,:)
   LIODBLE, allocatable :: ke_eorb(:)
   LIODBLE, allocatable :: PFW_vec(:)
   LIODBLE, allocatable :: phon_pop(:)
   LIODBLE, parameter   :: emass_d = 5.485799111d-4 !Electron mass in
                                                         !Daltons
   LIODBLE, parameter   :: freq2cm = 5152.792213813936d0!5141.1182009496200901D0 ! freq to cm-1
   type(cumat_r)        :: XCmat_ke
   type(cumat_x)        :: YCinv_ke

end module vib_KE_data
