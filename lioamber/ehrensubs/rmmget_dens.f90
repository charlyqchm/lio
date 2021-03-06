!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine rmmget_dens_r(DensMao)
  use garcha_mod, only: Pmat_vec
  use basis_data, only: M

  implicit none
  real*4,intent(out)     :: DensMao(M,M)
  integer                :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      DensMao(ii,jj)=Pmat_vec(idx)/2.0d0
      DensMao(jj,ii)=Pmat_vec(idx)/2.0d0
    enddo
    DensMao(jj,jj)=DensMao(jj,jj)*2.0d0
  enddo

  return;end subroutine
!--------------------------------------------------------------------!
  subroutine rmmget_dens_d(DensMao)
  use garcha_mod, only: Pmat_vec
  use basis_data, only: M

  implicit none
  real*8,intent(out)     :: DensMao(M,M)
  integer                :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      DensMao(ii,jj)=Pmat_vec(idx)/2.0d0
      DensMao(jj,ii)=Pmat_vec(idx)/2.0d0
    enddo
    DensMao(jj,jj)=DensMao(jj,jj)*2.0d0
  enddo

  return;end subroutine
!--------------------------------------------------------------------!
  subroutine rmmget_dens_c(DensMao)
  use garcha_mod, only: Pmat_vec
  use basis_data, only: M

  implicit none
  complex*8,intent(out)  :: DensMao(M,M)
  integer                :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      DensMao(ii,jj)=cmplx(Pmat_vec(idx)/2.0d0)
      DensMao(jj,ii)=cmplx(Pmat_vec(idx)/2.0d0)
    enddo
    DensMao(jj,jj)=DensMao(jj,jj)*2.0d0
  enddo

  return;end subroutine
!--------------------------------------------------------------------!
  subroutine rmmget_dens_z(DensMao)
  use garcha_mod, only: Pmat_vec
  use basis_data, only: M

  implicit none
  complex*16,intent(out) :: DensMao(M,M)
  integer                :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      DensMao(ii,jj)=dcmplx(Pmat_vec(idx)/2.0d0)
      DensMao(jj,ii)=dcmplx(Pmat_vec(idx)/2.0d0)
    enddo
    DensMao(jj,jj)=DensMao(jj,jj)*2.0d0
  enddo

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
