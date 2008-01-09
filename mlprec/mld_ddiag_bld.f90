!!$ 
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$	                 Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MLD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$
! File: mld_ddiag_bld.f90
!
! Subroutine: mld_ddiag_bld
! Version:    real
!
!  This routine builds the diagonal preconditioner corresponding to a given
!  sparse matrix A.    
!
!
! Arguments:
!    a       -  type(psb_dspmat_type), input.
!               The sparse matrix structure containing the local part of the
!               matrix A to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor associated to the sparse matrix A.
!    p       -  type(mld_dbaseprc_type), input/output.
!               The 'base preconditioner' data structure containing the local 
!               part of the diagonal preconditioner.
!    info    -  integer, output.
!               Error code.
!  
subroutine mld_ddiag_bld(a,desc_a,p,upd,info)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_ddiag_bld

  Implicit None

! Arguments
  type(psb_dspmat_type), target           :: a
  type(psb_desc_type), intent(in)         :: desc_a
  type(mld_dbaseprc_type),intent(inout)   :: p
  character, intent(in)                   :: upd
  integer, intent(out)                    :: info

! Local variables
  Integer      :: err, n_row, n_col,I,j,k,ictxt,&
       & me,np,mglob,lw, err_act
  integer             :: debug_level, debug_unit
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  name  = 'mld_ddiag_bld'
  info  = 0
  ictxt = psb_cd_get_context(desc_a)
  n_row = psb_cd_get_local_rows(desc_a)
  n_col = psb_cd_get_local_cols(desc_a)
  mglob = psb_cd_get_global_rows(desc_a)
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),' Enter'

  call psb_realloc(n_col,p%d,info)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_realloc')
    goto 9999
  end if

  !
  ! Retrieve the diagonal entries of the matrix A
  !
  call psb_sp_getdiag(a,p%d,info)
  if(info /= 0) then
    info=4010
    ch_err='psb_sp_getdiag'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  !
  ! Copy into p%desc_data the descriptor associated to A
  !
  call psb_cdcpy(desc_a,p%desc_Data,info)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_cdcpy')
    goto 9999
  end if

  !
  ! The i-th diagonal entry of the preconditioner is set to one if the
  ! corresponding entry a_ii of the sparse matrix A is zero; otherwise 
  ! it is set to one/a_ii
  !
  do i=1,n_row
    if (p%d(i) == dzero) then
      p%d(i) = done
    else
      p%d(i) = done/p%d(i)
    endif
  end do

  if (a%pl(1) /= 0) then
    !
    ! Apply the same row permutation as in the sparse matrix A
    !
    call  psb_gelp('n',a%pl,p%d,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_gelp'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  endif

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),'Done'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_ddiag_bld

