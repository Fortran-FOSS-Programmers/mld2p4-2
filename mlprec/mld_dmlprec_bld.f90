!!$ 
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$                      Alfredo Buttari      CNRS-IRIT, Toulouse
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
! File: mld_dmlprec_bld.f90
!
! Subroutine: mld_dmlprec_bld
! Version:    real
!
!  This routine builds the preconditioner according to the requirements made by
!  the user trough the subroutines mld_precinit and mld_precset.
!  
!  A multilevel preconditioner is regarded as an array of 'one-level' data structures,
!  each containing the part of the preconditioner associated to a certain level,
!  (for more details see the description of mld_Tonelev_type in mld_prec_type.f90).
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level. No transfer operators are associated to level 1.
! 
!
! Arguments:
!    a       -  type(psb_dspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!  
subroutine mld_dmlprec_bld(a,desc_a,p,info)

  use psb_sparse_mod
  use mld_inner_mod, mld_protect_name => mld_dmlprec_bld
  use mld_prec_mod
  use mld_d_jac_smoother
  use mld_d_as_smoother
  use mld_d_diag_solver
  use mld_d_ilu_solver
  use mld_d_umf_solver

  Implicit None

  ! Arguments
  type(psb_dspmat_type),intent(in), target  :: a
  type(psb_desc_type), intent(in), target    :: desc_a
  type(mld_dprec_type),intent(inout),target  :: p
  integer, intent(out)                       :: info
!!$  character, intent(in), optional         :: upd

  ! Local Variables
  type(mld_dprec_type)    :: t_prec
  Integer      :: err,i,k,ictxt, me,np, err_act, iszv, newsz
  integer      :: ipv(mld_ifpsz_), val
  integer      :: int_err(5)
  character    :: upd_
  integer            :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  if (psb_get_errstatus().ne.0) return 
  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'mld_dmlprec_bld'
  info = psb_success_
  int_err(1) = 0
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
  !
  ! For the time being we are commenting out the UPDATE argument
  ! we plan to resurrect it later. 
!!$  if (present(upd)) then 
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),'UPD ', upd
!!$
!!$    if ((psb_toupper(upd).eq.'F').or.(psb_toupper(upd).eq.'T')) then
!!$      upd_=psb_toupper(upd)
!!$    else
!!$      upd_='F'
!!$    endif
!!$  else
!!$    upd_='F'
!!$  endif
  upd_ = 'F'

  if (.not.allocated(p%precv)) then 
    !! Error: should have called mld_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  !
  ! Check to ensure all procs have the same 
  !   
  newsz = -1
  iszv  = size(p%precv)
  call psb_bcast(ictxt,iszv)
  if (iszv /= size(p%precv)) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of precv')
    goto 9999
  end if

  if (iszv <= 1) then
    ! We should only ever get here for multilevel.
    info=psb_err_from_subroutine_
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  if (iszv > 1) then

    !
    ! Build the matrix and the transfer operators corresponding
    ! to the remaining levels
    !
    !
    ! Check on the iprcparm contents: they should be the same
    ! on all processes.
    !
    if (me == psb_root_) ipv(:) = p%precv(1)%iprcparm(:) 
    call psb_bcast(ictxt,ipv) 
    if (any(ipv(:) /=  p%precv(1)%iprcparm(:) )) then
      write(debug_unit,*) me,name,&
           &': Inconsistent arguments among processes, forcing a default'
      p%precv(1)%iprcparm(:) = ipv(:) 
    end if
    !
    ! Finest level first; remember to fix base_a and base_desc
    ! 
    p%precv(1)%base_a    => a
    p%precv(1)%base_desc => desc_a

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,a_err='Base level precbuild.')
      goto 9999
    end if


    do i=2, iszv
      !
      ! Check on the iprcparm contents: they should be the same
      ! on all processes.
      !
      if (me == psb_root_) ipv(:) = p%precv(i)%iprcparm(:) 
      call psb_bcast(ictxt,ipv) 
      if (any(ipv(:) /=  p%precv(i)%iprcparm(:) )) then
        write(debug_unit,*) me,name,&
             &': Inconsistent arguments among processes, resetting.'
        p%precv(i)%iprcparm(:) = ipv(:) 
      end if

      !
      ! Sanity checks on the parameters
      !
      if (i<iszv) then 
        !
        ! A replicated matrix only makes sense at the coarsest level
        !
        call mld_check_def(p%precv(i)%iprcparm(mld_coarse_mat_),'Coarse matrix',&
             &   mld_distr_mat_,is_distr_ml_coarse_mat)

      else if (i == iszv) then 

        call check_coarse_lev(p%precv(i)) 

      end if

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Calling mlprcbld at level  ',i
      !
      ! Build the mapping between levels i-1 and i and the matrix
      ! at level i
      ! 
      if (info == psb_success_) call mld_coarse_bld(p%precv(i-1)%base_a,&
           & p%precv(i-1)%base_desc, p%precv(i),info)

      if (info /= psb_success_) then 
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Init upper level preconditioner')
        goto 9999
      endif

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Return from ',i,' call to mlprcbld ',info      

      if (i>2) then 
        if (all(p%precv(i)%map%naggr == p%precv(i-1)%map%naggr)) then 
          newsz=i-1
        end if
        call psb_bcast(ictxt,newsz)
        if (newsz > 0) exit
      end if
    end do

    if (newsz > 0) then 
      if (me == 0) then 
        write(debug_unit,*) trim(name),&
             &': Warning: aggregates from level ',&
             & newsz
        write(debug_unit,*) trim(name),&
             &':                       to level ',&
             & iszv,' coincide.'
        write(debug_unit,*) trim(name),&
             &': Number of levels actually used :',newsz
        write(debug_unit,*)
      end if
      allocate(t_prec%precv(newsz),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,&
             & a_err='prec reallocation')
        goto 9999
      endif
      do i=1,newsz-1
        call mld_move_alloc(p%precv(i),t_prec%precv(i),info)
      end do
      call mld_move_alloc(p%precv(iszv),t_prec%precv(newsz),info)
      do i=newsz+1, iszv
        call mld_precfree(p%precv(i),info)
      end do
      call mld_move_alloc(t_prec,p,info) 
      ! Ignore errors from transfer
      info = psb_success_
      !
      ! Restart
      iszv = newsz
      ! Fix the pointers, but the level 1 should
      ! be already OK
      do i=2, iszv - 1 
        p%precv(i)%base_a    => p%precv(i)%ac
        p%precv(i)%base_desc => p%precv(i)%desc_ac
        p%precv(i)%map%p_desc_X => p%precv(i-1)%base_desc
        p%precv(i)%map%p_desc_Y => p%precv(i)%base_desc
      end do


      i    = iszv 
      call check_coarse_lev(p%precv(i)) 
      if (info == psb_success_) call mld_coarse_bld(p%precv(i-1)%base_a,&
           & p%precv(i-1)%base_desc, p%precv(i),info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='coarse rebuild')
        goto 9999
      endif
    end if
  end if

  do i=1, iszv
    !
    ! build the base preconditioner at level i
    !
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Calling mlprcbld at level  ',i
    select case(p%precv(i)%prec%iprcparm(mld_sub_solve_))
    case(mld_ilu_n_,mld_milu_n_)      
      call mld_check_def(p%precv(i)%prec%iprcparm(mld_sub_fillin_),&
           & 'Level',0,is_legal_ml_lev)
    case(mld_ilu_t_)                 
      call mld_check_def(p%precv(i)%prec%rprcparm(mld_sub_iluthrs_),&
           & 'Eps',dzero,is_legal_fact_thrs)
    end select
    call mld_check_def(p%precv(i)%iprcparm(mld_smoother_sweeps_),&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)
    call mld_check_def(p%precv(i)%iprcparm(mld_smoother_sweeps_pre_),&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)
    call mld_check_def(p%precv(i)%iprcparm(mld_smoother_sweeps_post_),&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)

    !
    !  Test version for beginning of OO stuff. 
    ! 
    if (allocated(p%precv(i)%sm)) then 
      call p%precv(i)%sm%free(info)
      if (info == psb_success_) deallocate(p%precv(i)%sm,stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_alloc_dealloc_,name,&
             & a_err='One level preconditioner build.')
        goto 9999
      endif
    end if
    select case (p%precv(i)%prec%iprcparm(mld_smoother_type_)) 
    case(mld_bjac_, mld_jac_) 
      allocate(mld_d_jac_smoother_type :: p%precv(i)%sm, stat=info)
    case(mld_as_)
      allocate(mld_d_as_smoother_type  :: p%precv(i)%sm, stat=info)
    case default
      info = -1 
    end select
    if (info /= psb_success_) then 
      write(0,*) ' Smoother allocation error',info,&
           & p%precv(i)%prec%iprcparm(mld_smoother_type_)
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='One level preconditioner build.')
      goto 9999
    endif
    call p%precv(i)%sm%set(mld_sub_restr_,p%precv(i)%prec%iprcparm(mld_sub_restr_),info)
    call p%precv(i)%sm%set(mld_sub_prol_,p%precv(i)%prec%iprcparm(mld_sub_prol_),info)
    call p%precv(i)%sm%set(mld_sub_ovr_,p%precv(i)%prec%iprcparm(mld_sub_ovr_),info)

    select case (p%precv(i)%prec%iprcparm(mld_sub_solve_)) 
    case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 
      allocate(mld_d_ilu_solver_type :: p%precv(i)%sm%sv, stat=info)
      if (info == psb_success_) call  p%precv(i)%sm%sv%set(mld_sub_solve_,&
           & p%precv(i)%prec%iprcparm(mld_sub_solve_),info)
      if (info == psb_success_) call  p%precv(i)%sm%sv%set(mld_sub_fillin_,&
           & p%precv(i)%prec%iprcparm(mld_sub_fillin_),info)
      if (info == psb_success_) call  p%precv(i)%sm%sv%set(mld_sub_iluthrs_,&
           & p%precv(i)%prec%rprcparm(mld_sub_iluthrs_),info)
    case(mld_diag_scale_)
      allocate(mld_d_diag_solver_type :: p%precv(i)%sm%sv, stat=info)
    case(mld_umf_)
      allocate(mld_d_umf_solver_type :: p%precv(i)%sm%sv, stat=info)
    case default
      info = -1 
    end select

    if (info /= psb_success_) then 
      write(0,*) ' Solver allocation error',info,&
           & p%precv(i)%prec%iprcparm(mld_sub_solve_)
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='One level preconditioner build.')
      goto 9999
    endif

    call p%precv(i)%sm%build(p%precv(i)%base_a,p%precv(i)%base_desc,'F',info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='One level preconditioner build.')
      goto 9999
    endif

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Return from ',i,' call to mlprcbld ',info      
  end do


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains

  subroutine check_coarse_lev(prec)
    type(mld_donelev_type) :: prec

    !
    ! At the coarsest level, check mld_coarse_solve_ 
    !
    val = prec%iprcparm(mld_coarse_solve_)  
    select case (val) 
    case(mld_jac_)   
      
      if (prec%prec%iprcparm(mld_sub_solve_) /= mld_diag_scale_) then 
        if (me == 0) write(debug_unit,*)&
             & 'Warning: inconsistent coarse level specification.'
        if (me == 0) write(debug_unit,*)&
             & '         Resetting according to the value specified for mld_coarse_solve_.'
        prec%prec%iprcparm(mld_sub_solve_)    = mld_diag_scale_
      end if
      prec%prec%iprcparm(mld_smoother_type_)  = mld_jac_          
      
    case(mld_bjac_)   
      
      if ((prec%prec%iprcparm(mld_sub_solve_) == mld_diag_scale_).or.&
           & ( prec%prec%iprcparm(mld_smoother_type_) /= mld_bjac_))  then 
        if (me == 0) write(debug_unit,*)&
             & 'Warning: inconsistent coarse level specification.'
        if (me == 0) write(debug_unit,*)&
             & '         Resetting according to the value specified for mld_coarse_solve_.'
!!$#if defined(HAVE_UMF_)
!!$        prec%prec%iprcparm(mld_sub_solve_)       = mld_umf_
!!$#elif defined(HAVE_SLU_)
!!$        prec%prec%iprcparm(mld_sub_solve_)       = mld_slu_
!!$#else 
        prec%prec%iprcparm(mld_sub_solve_)       = mld_ilu_n_
!!$#endif
      end if
      prec%prec%iprcparm(mld_smoother_type_)  = mld_bjac_          
      
    case(mld_umf_, mld_slu_)
      if ((prec%iprcparm(mld_coarse_mat_)  /= mld_repl_mat_).or.&
           & (prec%prec%iprcparm(mld_sub_solve_)  /= val)) then 
        if (me == 0) write(debug_unit,*)&
             & 'Warning: inconsistent coarse level specification.'
        if (me == 0) write(debug_unit,*)&
             & '         Resetting according to the value specified for mld_coarse_solve_.'
        prec%iprcparm(mld_coarse_mat_)         = mld_repl_mat_
        prec%prec%iprcparm(mld_sub_solve_)     = val
        prec%prec%iprcparm(mld_smoother_type_) = mld_bjac_          
      end if
    case(mld_sludist_)
      if ((prec%iprcparm(mld_coarse_mat_)  /= mld_distr_mat_).or.&
           & (prec%prec%iprcparm(mld_sub_solve_)  /= val)) then 
        if (me == 0) write(debug_unit,*)&
             & 'Warning: inconsistent coarse level specification.'
        if (me == 0) write(debug_unit,*)&
             & '         Resetting according to the value specified for mld_coarse_solve_.'
        prec%iprcparm(mld_coarse_mat_)           = mld_distr_mat_
        prec%prec%iprcparm(mld_sub_solve_)       = val
        prec%prec%iprcparm(mld_smoother_type_)   = mld_bjac_          
        prec%prec%iprcparm(mld_smoother_sweeps_) = 1
      end if
    end select
  end subroutine check_coarse_lev

end subroutine mld_dmlprec_bld
