!carlos: for the moment these subroutines only conmut the ON and real type matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Commut_data_rr(this, Bmat, AB_BAmat, Nsize)

#ifdef CUBLAS
   use cublasmath, only : commutator_cublas
#else
   use mathsubs,   only: commutator
#endif

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Nsize
   real*8, intent(in)             :: Bmat(Nsize,Nsize)
   real*8, intent(inout)          :: AB_BAmat(Nsize,Nsize)

   real*8, allocatable :: Amat(:,:)
   allocate(Amat(Nsize,Nsize))

   Amat=this%data_ON

#ifdef CUBLAS
      AB_BAmat = commutator_cublas(Amat, Bmat)
#else
      AB_BAmat = commutator (Amat, Bmat)
#endif

end subroutine Commut_data_rr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Commut_data_rc(this, Bmat, AB_BAmat, Nsize)

#ifdef CUBLAS
   use cublasmath, only : commutator_cublas
#else
   use mathsubs,   only: commutator
#endif

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Nsize

   TDCOMPLEX, intent(in)       :: Bmat(Nsize,Nsize)
   TDCOMPLEX, intent(out)      :: AB_BAmat(Nsize,Nsize)

   real*8, allocatable         :: Amat(:,:)

   allocate(Amat(Nsize,Nsize))

   Amat=this%data_ON

#ifdef CUBLAS
      AB_BAmat = commutator_cublas(Amat, Bmat)
#else
      AB_BAmat = commutator (Amat, Bmat)
#endif

end subroutine Commut_data_rc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine Commut_data_cr(this, Bmat, AB_BAmat, Nsize)

#ifdef CUBLAS
   use cublasmath, only : commutator_cublas
#else
   use mathsubs,   only: commutator
#endif

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Nsize
   real*8, intent(in)             :: Bmat(Nsize,Nsize)
   TDCOMPLEX, intent(inout)       :: AB_BAmat(Nsize,Nsize)

   TDCOMPLEX, allocatable  :: Amat(:,:)

   allocate(Amat(Nsize,Nsize))

   Amat=this%dataC_ON

#ifdef CUBLAS
      AB_BAmat = commutator_cublas(Amat, Bmat)
#else
      AB_BAmat = commutator (Amat, Bmat)
#endif

end subroutine Commut_data_cr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine Commut_data_cc(this, Bmat, AB_BAmat, Nsize)

#ifdef CUBLAS
   use cublasmath, only : commutator_cublas
#else
   use mathsubs,   only: commutator
#endif

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Nsize

   TDCOMPLEX, intent(in)       :: Bmat(Nsize,Nsize)
   TDCOMPLEX, intent(out)      :: AB_BAmat(Nsize,Nsize)

   TDCOMPLEX, allocatable      :: Amat(:,:)

   allocate(Amat(Nsize,Nsize))

   Amat=this%dataC_ON

#ifdef CUBLAS
      AB_BAmat = commutator_cublas(Amat, Bmat)
#else
      AB_BAmat = commutator (Amat, Bmat)
#endif

end subroutine Commut_data_cc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
