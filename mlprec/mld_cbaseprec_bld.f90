!!$ 
!!$ 
!!$                           MLD2P4  version 1.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.2)
!!$  
!!$  (C) Copyright 2008
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
!!$                      Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
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
! File: mld_cbaseprec_bld.f90
!
! Subroutine: mld_cbaseprc_bld
! Version:    complex
!
!  This routine builds a 'base preconditioner' related to a matrix A.
!  In a multilevel framework, it is called by mld_mlprec_bld to build the
!  base preconditioner at each level.
!
!  Details on the base preconditioner to be built are stored in the iprcparm
!  field of the preconditioner data structure (for a description of this
!  data structure see mld_prec_type.f90).
!    
!
! Arguments:
!    a       -  type(psb_cspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix A to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_cbaseprec_type), input/output.
!               The 'base preconditioner' data structure containing the local 
!               part of the preconditioner at the selected level.
!    info    -  integer, output.
!               Error code.              
!    upd     -  character, input, optional.
!               If upd='F' then the base preconditioner is built from
!               scratch; if upd=T' then the matrix to be preconditioned
!               has the same sparsity pattern of a matrix that has been
!               previously preconditioned, hence some information is reused
!               in building the new preconditioner.
!  
subroutine mld_cbaseprc_bld(a,desc_a,p,info,upd)

  use psb_base_mod
  use mld_inner_mod, mld_protect_name => mld_cbaseprc_bld

  Implicit None

  ! Arguments
  type(psb_cspmat_type), target           :: a
  type(psb_desc_type), intent(in), target :: desc_a
  type(mld_cbaseprc_type),intent(inout)   :: p
  integer, intent(out)                    :: info
  character, intent(in), optional         :: upd

  ! Local variables
  Integer      :: err, n_row, n_col,ictxt, me,np,mglob, err_act
  character    :: iupd
  integer             :: debug_level, debug_unit
  character(len=20)   :: name, ch_err

  if (psb_get_errstatus() /= 0) return 
  name = 'mld_cbaseprc_bld'
  info=0
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt   = psb_cd_get_context(desc_a)
  n_row   = psb_cd_get_local_rows(desc_a)
  n_col   = psb_cd_get_local_cols(desc_a)
  mglob   = psb_cd_get_global_rows(desc_a)
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  if (present(upd)) then 
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),'UPD ', upd
    if ((psb_toupper(UPD) ==  'F').or.(psb_toupper(UPD) == 'T')) then
      IUPD=psb_toupper(UPD)
    else
      IUPD='F'
    endif
  else
    IUPD='F'
  endif

  !
  ! Should add check to ensure all procs have the same... 
  !

  call mld_check_def(p%iprcparm(mld_prec_type_),'base_prec',&
       &  mld_diag_,is_legal_base_prec)


  call psb_nullify_desc(p%desc_data)

  select case(p%iprcparm(mld_prec_type_))

  case (mld_noprec_)
    ! No preconditioner 

    ! Do nothing 
    call psb_cdcpy(desc_a,p%desc_data,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_cdcpy'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case (mld_diag_)
    ! Diagonal preconditioner

    call mld_diag_bld(a,desc_a,p,info)
    if(debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ': out of mld_diag_bld'
    if(info /= 0) then
      info=4010
      ch_err='mld_diag_bld'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(mld_bjac_,mld_as_)
    ! Additive Schwarz preconditioners/smoothers

    call mld_check_def(p%iprcparm(mld_n_ovr_),'overlap',&
         &  0,is_legal_n_ovr)
    call mld_check_def(p%iprcparm(mld_sub_restr_),'restriction',&
         &  psb_halo_,is_legal_restrict)
    call mld_check_def(p%iprcparm(mld_sub_prol_),'prolongator',&
         &  psb_none_,is_legal_prolong)
    call mld_check_def(p%iprcparm(mld_sub_ren_),'renumbering',&
         &  mld_renum_none_,is_legal_renum)
    call mld_check_def(p%iprcparm(mld_sub_solve_),'fact',&
         &  mld_ilu_n_,is_legal_ml_fact)

    ! Set parameters for using SuperLU_dist on the local submatrices
    if (p%iprcparm(mld_sub_solve_)==mld_sludist_) then
      p%iprcparm(mld_n_ovr_)         = 0
      p%iprcparm(mld_smooth_sweeps_) = 1
    end if

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ': Calling mld_as_bld'

    ! Build the local part of the base preconditioner/smoother
    call mld_as_bld(a,desc_a,p,iupd,info)
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='mld_as_bld')
      goto 9999
    end if

  case default

    info=4001
    ch_err='Unknown mld_prec_type_'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999

  end select

  p%base_a    => a
  p%base_desc => desc_a
  p%iprcparm(mld_prec_status_) = mld_prec_built_
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Done'
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_cbaseprc_bld

