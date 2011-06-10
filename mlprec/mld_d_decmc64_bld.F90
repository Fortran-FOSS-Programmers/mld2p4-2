subroutine mld_d_decmc64_bld(theta,a,desc_a,nlaggr,ilaggr,info)

  use psb_base_mod
  use mld_mc64_tools_mod
  use mld_d_inner_mod, mld_protect_name => mld_d_decmc64_bld

  implicit none

  ! Arguments
  type(psb_dspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)    :: desc_a
  real(psb_dpk_), intent(in)         :: theta
  integer, allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
  integer, intent(out)               :: info

  ! Local variables
  integer, allocatable        :: ils(:), neigh(:), irow(:), icol(:)
  real(psb_dpk_), allocatable :: val(:), diag(:)
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
  integer :: nrow, ncol, n_ne, nnz, job, nr, nc, num, ir, ic, ip, nr2, nc2, nz2, nr3, nc4, nz3
  integer           :: irmax, icmax, idx, mxk2
  real(psb_dpk_)    :: rmax, cmax
  character(len=20) :: name, ch_err
  integer :: icntl(10), infov(10)
  real(psb_dpk_) :: dcntl(10) 
  integer :: liw, ldw
  integer, allocatable        :: mark(:), mark2(:), mark3(:)

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
  nr  = a%get_nrows()
  nnz = a%get_nzeros()
  call mld_d_match(a,nagg1,mark,info)
  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Fatal error: naggr>ncol')
    goto 9999
  end if

  ! Now call matching a second time

  ! First prepare an aggregated matrix
  call a%csclip(acoo,info,jmax=nr)
  nz2 = acoo%get_nzeros()
  do i=1, nz2
    acoo%ia(i) = mark(acoo%ia(i))
    acoo%ja(i) = mark(acoo%ja(i))
  end do
  call acoo%set_nrows(nagg1)
  call acoo%set_ncols(nagg1)
  call acoo%set_dupl(psb_dupl_add_)
  call acoo%fix(info)
  call atmp%mv_from(acoo)
  call atmp%clip_diag(info)
  
  call mld_d_match(atmp,nagg2,mark2,info)
  
  if (.false.) then 
    call atmp%mv_to(acoo)
    nz3 = acoo%get_nzeros()
    do i=1, nz3
      acoo%ia(i) = mark2(acoo%ia(i))
      acoo%ja(i) = mark2(acoo%ja(i))
    end do
    call acoo%set_nrows(nagg2)
    call acoo%set_ncols(nagg2)
    call acoo%set_dupl(psb_dupl_add_)
    call acoo%fix(info)
    call atmp%mv_from(acoo)
    call atmp%clip_diag(info)
    call mld_d_match(atmp,nagg3,mark3,info)
    
    naggr = nagg3
    do i=1, nr
      ilaggr(i) = mark3(mark2(mark(i)))
    end do
  else 
    naggr = nagg2
    do i=1, nr
      ilaggr(i) = (mark2(mark(i)))
    end do
  end if

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

