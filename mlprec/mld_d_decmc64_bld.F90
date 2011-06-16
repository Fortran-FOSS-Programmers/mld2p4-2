subroutine mld_d_decmc64_bld(nsteps,theta,a,desc_a,nlaggr,ilaggr,info)

  use psb_base_mod
  use mld_mc64_tools_mod
  use mld_d_inner_mod, mld_protect_name => mld_d_decmc64_bld

  implicit none

  ! Arguments
  integer, intent(in)               :: nsteps
  type(psb_dspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)   :: desc_a
  real(psb_dpk_), intent(in)        :: theta
  integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  integer, intent(out)              :: info

  ! Local variables
  type(psb_dspmat_type)       :: atmp
  type(psb_d_csc_sparse_mat)  :: acsc, acred
  type(psb_d_csr_sparse_mat)  :: acsr
  type(psb_d_coo_sparse_mat)  :: acoo
  integer :: icnt,nlp,k,n,ia,isz,naggr,i,j,m, nz, nagg1, nagg2, nagg3
  integer :: num_agg, num_unagg
  real(psb_dpk_)  :: cpling, tcl
  logical :: recovery
  integer :: debug_level, debug_unit
  integer :: ictxt,np,me,err_act
  integer :: nrow, ncol, nr, nc, nnz
  integer           :: irmax, icmax, idx, mxk2, istep
  real(psb_dpk_)    :: rmax, cmax
  character(len=20) :: name, ch_err
  integer :: icntl(10), infov(10)
  real(psb_dpk_) :: dcntl(10) 
  integer              :: liw, ldw
  integer, allocatable :: mark(:), mark2(:), mark3(:)
  real(psb_dpk_), allocatable :: diag(:)
  type :: match_data_type
    integer              :: nagg=0
    integer, allocatable :: mark(:)
  end type match_data_type
  type(match_data_type), allocatable :: match_data(:)

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  name = 'mld_dec_map_bld'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  !
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt,me,np)
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

#ifdef HAVE_MC64_
  call psb_realloc(ncol,ilaggr,info)
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  allocate(match_data(nsteps),stat=info) 
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/nsteps,0,0,0,0/),&
         & a_err='match_data')
    goto 9999
  end if

  nr = a%get_nrows()
  call a%csclip(atmp,info,jmax=nr)
  call atmp%mv_to(acoo)
  acoo%val=abs(acoo%val)
  call atmp%mv_from(acoo)
  call psb_realloc(nr,diag,info)
  call atmp%get_diag(diag,info)
  call atmp%clip_diag(info)

  do istep = 1, nsteps
    if (istep >1) then 
      call atmp%mv_to(acoo)
      ! First prepare an aggregated matrix
      nz = acoo%get_nzeros()
      nr = acoo%get_nrows()
      ! Put back the diagonal
      call psb_ensure_size(nz+nr,acoo%ia,info)
      call psb_ensure_size(nz+nr,acoo%ja,info)
      call psb_ensure_size(nz+nr,acoo%val,info)
      do i=1, nr
        acoo%ia(nz+i)  = i
        acoo%ja(nz+i)  = i
        acoo%val(nz+i) = diag(i)
      end do
      call acoo%set_sorted(.false.)
      call acoo%set_nzeros(nz+nr)
      nz = acoo%get_nzeros()
      do i=1, nz
        acoo%ia(i) = match_data(istep-1)%mark(acoo%ia(i))
        acoo%ja(i) = match_data(istep-1)%mark(acoo%ja(i))
      end do
      call acoo%set_nrows(match_data(istep-1)%nagg)
      call acoo%set_ncols(match_data(istep-1)%nagg)
      call acoo%set_dupl(psb_dupl_add_)
      call acoo%fix(info)
      call atmp%mv_from(acoo)
    endif
    ! Take the diagonal out of the matrix
    call atmp%get_diag(diag,info)
    call atmp%clip_diag(info)
    if (.true.) then 
      call mld_d_match_new(theta,diag,atmp,match_data(istep)%nagg,&
           & match_data(istep)%mark,info)
    else
      call mld_d_match_old(atmp,match_data(istep)%nagg,&
           & match_data(istep)%mark,info)
    end if

  end do

  !
  ! Go back to initial NR
  !
  nr    = a%get_nrows()
  naggr = match_data(nsteps)%nagg
  do i=1, nr
    ilaggr(i) = i
    do istep=1, nsteps
      ilaggr(i) = match_data(istep)%mark(ilaggr(i))
    end do
  end do
  do i=1, nr
    if ((ilaggr(i)<1).or.(ilaggr(i)>naggr)) then 
      write(0,*) 'Invalid entry ',i,ilaggr(i)
    end if
  end do
  
  call atmp%free()


  allocate(nlaggr(np),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/np,0,0,0,0/),&
         & a_err='integer')
    goto 9999
  end if

  nlaggr(:) = 0
  nlaggr(me+1) = naggr
  call psb_sum(ictxt,nlaggr(1:np))


  call psb_erractionrestore(err_act)
  return
#else 

  write(psb_err_unit,*) 'MC64 matching support not configured, aborting'
  info = psb_err_internal_error_
  call psb_errpush(info,name)

  goto 9999

#endif
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine mld_d_decmc64_bld

