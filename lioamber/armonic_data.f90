module armonic_data
   implicit none
!Important input variables
! * armonic_calc : integer, indicates the different calculations options()

   logical     , allocatable :: move_atom(:)
   integer                   :: armonic_calc = 1
   integer                   :: hess_norder  = 1
   real(kind=8)              :: delta_h      = 0.01d0 !Amstrongs
   real(kind=8), allocatable :: atom_mass(:)
   real(kind=8), allocatable :: mass_w(:)
   real(kind=8), allocatable :: armonic_freq(:)
   real(kind=8), allocatable :: armonic_vec(:,:)
   real(kind=8), parameter   :: emass_d = 5.485799111d-4 !Electron mass in
                                                         !Daltons
   real(kind=8), parameter   :: freq2cm = 5141.1182009496200901D0 ! freq to cm-1


end module armonic_data
