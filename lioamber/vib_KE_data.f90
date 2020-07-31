#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module vib_KE_data
   implicit none

   logical              :: vib_calc = .true.
   integer              :: ke_calc  = 0
   integer, allocatable :: move_atom(:)
   integer              :: tot_at_move
   integer              :: hess_norder  = 1
   LIODBLE              :: delta_h      = 0.01d0 !Amstrongs
   LIODBLE, allocatable :: atom_mass(:)
   LIODBLE, allocatable :: mass_w(:)
   LIODBLE, allocatable :: armonic_freq(:)
   LIODBLE, allocatable :: armonic_vec(:,:)
   LIODBLE, parameter   :: emass_d = 5.485799111d-4 !Electron mass in
                                                         !Daltons
   LIODBLE, parameter   :: freq2cm = 5152.792213813936d0!5141.1182009496200901D0 ! freq to cm-1

end module vib_KE_data
