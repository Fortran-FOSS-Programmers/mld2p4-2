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
! File: mld_sprecset.f90
!
! Subroutine: mld_sprecseti
! Version: real
!
!  This routine sets the integer parameters defining the preconditioner. More
!  precisely, the integer parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set character and real parameters, see mld_sprecsetc and mld_sprecsetr,
!  respectively.
!
!
! Arguments:
!    p       -  type(mld_sprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in MLD2P4 user's guide.
!    val     -  integer, input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in MLD2P4 user's guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!   
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  MLD2P4 developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface mld_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see mld_prec_mod.f90).
!   
subroutine mld_sprecseti(p,what,val,info,ilev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_sprecseti

  implicit none

! Arguments
  type(mld_sprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  integer, intent(in)                    :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

! Local variables
  integer                                :: ilev_, nlev_
  character(len=*), parameter            :: name='mld_precseti'

  info = 0

  if (.not.allocated(p%baseprecv)) then 
    info = 3111
    write(0,*) name,': Error: uninitialized preconditioner, should call MLD_PRECINIT'
    return 
  endif
  nlev_ = size(p%baseprecv)

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    info = -1
    write(0,*) name,': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif
  if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
    info = 3111
    write(0,*) name,&
         &': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
    return 
  endif

  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then 
    
    if (ilev_ == 1) then
      ! 
      ! Rules for fine level are slightly different.
      ! 
      select case(what) 
      case(mld_smoother_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_,&
           & mld_sub_ren_,mld_sub_ovr_,mld_sub_fillin_,mld_smoother_sweeps_)
        p%baseprecv(ilev_)%iprcparm(what)  = val
      case default
        write(0,*) name,': Error: invalid WHAT'
        info = -2
      end select

    else if (ilev_ > 1) then 
      select case(what) 
      case(mld_smoother_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_,&
           & mld_sub_ren_,mld_sub_ovr_,mld_sub_fillin_,&
           & mld_smoother_sweeps_,mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,&
           & mld_smoother_pos_,mld_aggr_omega_alg_,mld_aggr_eig_)
        p%baseprecv(ilev_)%iprcparm(what)  = val
      case(mld_coarse_mat_)
        if (ilev_ /= nlev_ .and. val /= mld_distr_mat_) then 
          write(0,*) name,': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        p%baseprecv(ilev_)%iprcparm(mld_coarse_mat_)  = val
      case(mld_coarse_subsolve_)
        if (ilev_ /= nlev_) then 
          write(0,*) name,': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        p%baseprecv(ilev_)%iprcparm(mld_sub_solve_)  = val
      case(mld_coarse_solve_)
        if (ilev_ /= nlev_) then 
          write(0,*) name,': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if

        if (nlev_ > 1) then 
          p%baseprecv(nlev_)%iprcparm(mld_coarse_solve_)  = val
          p%baseprecv(nlev_)%iprcparm(mld_smoother_type_) = mld_bjac_
          p%baseprecv(nlev_)%iprcparm(mld_coarse_mat_)    = mld_distr_mat_
          select case (val) 
          case(mld_umf_, mld_slu_)
            p%baseprecv(nlev_)%iprcparm(mld_coarse_mat_)  = mld_repl_mat_
            p%baseprecv(nlev_)%iprcparm(mld_sub_solve_)   = val
          case(mld_sludist_)
            p%baseprecv(nlev_)%iprcparm(mld_sub_solve_)   = val
          end select
        endif
      case(mld_coarse_sweeps_)
        if (ilev_ /= nlev_) then 
          write(0,*) name,': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        p%baseprecv(ilev_)%iprcparm(mld_smoother_sweeps_)  = val
      case(mld_coarse_fillin_)
        if (ilev_ /= nlev_) then 
          write(0,*) name,': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        p%baseprecv(ilev_)%iprcparm(mld_sub_fillin_)  = val
      case default
        write(0,*) name,': Error: invalid WHAT'
        info = -2
      end select

    endif

  else if (.not.present(ilev)) then 
    !
    ! ilev not specified: set preconditioner parameters at all the appropriate
    ! levels
    !
    select case(what) 
    case(mld_smoother_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_,&
         & mld_sub_ren_,mld_sub_ovr_,mld_sub_fillin_,&
         & mld_smoother_sweeps_)
      do ilev_=1,max(1,nlev_-1)
        if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
          write(0,*) name,&
               &': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
          info = -1 
          return 
        endif
        p%baseprecv(ilev_)%iprcparm(what)  = val
      end do
    case(mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,&
         & mld_smoother_pos_,mld_aggr_omega_alg_,mld_aggr_eig_)
      do ilev_=2,nlev_
        if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
          write(0,*) name,&
               &': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
          info = -1 
          return 
        endif
        p%baseprecv(ilev_)%iprcparm(what)  = val
      end do
    case(mld_coarse_mat_)
      if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
        write(0,*) name,&
             & ': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
        info = -1 
        return 
      endif
      if (nlev_ > 1) p%baseprecv(nlev_)%iprcparm(mld_coarse_mat_)  = val
    case(mld_coarse_solve_)
      if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
        write(0,*) name,&
             &': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
        info = -1 
        return 
      endif

      if (nlev_ > 1) then 
        p%baseprecv(nlev_)%iprcparm(mld_coarse_solve_)  = val
        p%baseprecv(nlev_)%iprcparm(mld_smoother_type_) = mld_bjac_
        p%baseprecv(nlev_)%iprcparm(mld_coarse_mat_)    = mld_distr_mat_
        select case (val) 
        case(mld_umf_, mld_slu_)
          p%baseprecv(nlev_)%iprcparm(mld_coarse_mat_)  = mld_repl_mat_
          p%baseprecv(nlev_)%iprcparm(mld_sub_solve_)   = val
        case(mld_sludist_)
          p%baseprecv(nlev_)%iprcparm(mld_sub_solve_)   = val
        end select
      endif
    case(mld_coarse_subsolve_)
      if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
        write(0,*) name,&
             &': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
        info = -1 
        return 
      end if
      if (nlev_ > 1) p%baseprecv(nlev_)%iprcparm(mld_sub_solve_)  = val

    case(mld_coarse_sweeps_)
      if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
        write(0,*) name,&
             &': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
        info = -1 
        return 
      endif
      if (nlev_ > 1) p%baseprecv(nlev_)%iprcparm(mld_smoother_sweeps_)  = val
    case(mld_coarse_fillin_)
      if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
        write(0,*) name,&
             &': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
        info = -1 
        return 
      endif
      if (nlev_ > 1) p%baseprecv(nlev_)%iprcparm(mld_sub_fillin_)  = val
    case default
      write(0,*) name,': Error: invalid WHAT'
      info = -2
    end select

  endif

end subroutine mld_sprecseti

!
! Subroutine: mld_sprecsetc
! Version: real
! Contains: mld_stringval
!
!  This routine sets the character parameters defining the preconditioner. More
!  precisely, the character parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set integer and real parameters, see mld_sprecseti and mld_sprecsetr,
!  respectively.
!
!
! Arguments:
!    p       -  type(mld_sprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in MLD2P4 user's guide.
!    string  -  character(len=*), input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in MLD2P4 user's guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  MLD2P4 developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface mld_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see mld_prec_mod.f90).
!   
!   
subroutine mld_sprecsetc(p,what,string,info,ilev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_sprecsetc

  implicit none

  ! Arguments
  type(mld_sprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  character(len=*), intent(in)           :: string
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

! Local variables
  integer                                :: ilev_, nlev_,val
  character(len=*), parameter            :: name='mld_precseti'

  info = 0

  if (.not.allocated(p%baseprecv)) then 
    info = 3111
    return 
  endif
  nlev_ = size(p%baseprecv)

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    write(0,*) name,': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    info = -1
    return
  endif
  if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
    write(0,*) name,': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
    info = 3111
    return 
  endif


  call mld_stringval(string,val,info)
  if (info == 0) call mld_inner_precset(p,what,val,info,ilev=ilev)


end subroutine mld_sprecsetc


!
! Subroutine: mld_sprecsetr
! Version: real
!
!  This routine sets the real parameters defining the preconditioner. More
!  precisely, the real parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set integer and character parameters, see mld_sprecseti and mld_sprecsetc,
!  respectively.
!
! Arguments:
!    p       -  type(mld_sprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in MLD2P4 user's guide.
!    val     -  real(psb_spk_), input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in MLD2P4 user's guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  MLD2P4 developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface mld_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see mld_prec_mod.f90).
!   
!   
subroutine mld_sprecsetr(p,what,val,info,ilev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_sprecsetr

  implicit none

  ! Arguments
  type(mld_sprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  real(psb_spk_), intent(in)           :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

! Local variables
  integer                                :: ilev_,nlev_
  character(len=*), parameter            :: name='mld_precsetd'

  info = 0

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if (.not.allocated(p%baseprecv)) then 
    write(0,*) name,': Error: uninitialized preconditioner, should call MLD_PRECINIT' 
    info = 3111
    return 
  endif
  nlev_ = size(p%baseprecv)

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    write(0,*) name,': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    info = -1
    return
  endif
  if (.not.allocated(p%baseprecv(ilev_)%rprcparm)) then 
    write(0,*) name,': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
    info = 3111
    return 
  endif

  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then 
    
      if (ilev_ == 1) then 
        !
        ! Rules for fine level are slightly different. 
        !
        select case(what) 
        case(mld_sub_iluthrs_)
          p%baseprecv(ilev_)%rprcparm(what)  = val
        case default
          write(0,*) name,': Error: invalid WHAT'
          info = -2
        end select

      else if (ilev_ > 1) then 
        select case(what) 
        case(mld_aggr_omega_val_,mld_aggr_thresh_,mld_sub_iluthrs_)
          p%baseprecv(ilev_)%rprcparm(what)  = val
        case default
          write(0,*) name,': Error: invalid WHAT'
          info = -2
        end select
      endif

  else if (.not.present(ilev)) then 
      !
      ! ilev not specified: set preconditioner parameters at all the appropriate levels
      !

      select case(what) 
      case(mld_sub_iluthrs_)
        do ilev_=1,nlev_
          if (.not.allocated(p%baseprecv(ilev_)%rprcparm)) then 
            write(0,*) name,': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
            info = -1 
            return 
          endif
          p%baseprecv(ilev_)%rprcparm(what)  = val
        end do
      case(mld_coarse_iluthrs_)
        ilev_=nlev_
        if (.not.allocated(p%baseprecv(ilev_)%rprcparm)) then 
          write(0,*) name,': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
          info = -1 
          return 
        endif
        p%baseprecv(ilev_)%rprcparm(mld_sub_iluthrs_)  = val
      case(mld_aggr_omega_val_)
        do ilev_=2,nlev_
          if (.not.allocated(p%baseprecv(ilev_)%rprcparm)) then 
            write(0,*) name,': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
            info = -1 
            return 
          endif
          p%baseprecv(ilev_)%rprcparm(what)  = val
        end do
      case(mld_aggr_thresh_)
        do ilev_=2,nlev_
          if (.not.allocated(p%baseprecv(ilev_)%rprcparm)) then 
            write(0,*) name,': Error: uninitialized preconditioner component, should call MLD_PRECINIT' 
            info = -1 
            return 
          endif
          p%baseprecv(ilev_)%rprcparm(what)  = val
        end do
      case default
        write(0,*) name,': Error: invalid WHAT'
        info = -2
      end select

  endif

end subroutine mld_sprecsetr
