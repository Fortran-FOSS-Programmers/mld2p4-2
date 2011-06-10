module mld_mc64_tools_mod

#ifdef HAVE_MC64_
contains
  subroutine marking(nr,perm,mark,info)
    implicit none 
    integer, intent(in)  :: nr, perm(:)
    integer, intent(out) :: mark(:)

    integer, allocatable :: pinv(:)
    integer :: i,j,k, info, nagg

    allocate(pinv(nr),stat=info) 

    if (info /= 0) then 
      write(0,*) 'Allocation error'
      info = -1
      return
    end if
    pinv = 0
    do i=1, nr
      if ((perm(i) > 0).and.(perm(i)<=nr)) then 
        pinv(perm(i)) = i
      end if
    end do
    do i=1,nr
      if (pinv(i) == 0) pinv(i) = -3
    end do
    mark(1:nr) = 0
    nagg = 0
    do i=1, nr
      if (mark(i) == 0) then 
        if ((perm(i) <= 0).or.(perm(i)>nr)) then 
          ! Unmatched node
          mark(i) = -1
        else if (perm(i) == i) then 
          ! Self-matched node
          mark(i) = -2
        else 
          if (mark(perm(i)) <= 0) then 
            ! adjacent node unmarked
            nagg          = nagg + 1 
            mark(i)       = nagg
            mark(perm(i)) = nagg
          else 
            ! Adjacent node marked, check dual of I
            if (pinv(i) < 0) then 
              ! unmatched dual 
              mark(i) = -4
            else 
              if (mark(pinv(i))  == 0 ) then 
                nagg          = nagg + 1 
                mark(i)       = nagg
                mark(perm(i)) = nagg
              else 
                mark(i)    = -3
              end if

            end if
          end if
        end if
      end if
    end do


  end subroutine marking

  subroutine mld_d_match(a,nagg1,mark,info)
    use psb_base_mod
    implicit none 
    type(psb_dspmat_type), intent(in) :: a
    integer, intent(out)              :: nagg1 
    integer, allocatable, intent(out) :: mark(:)
    integer, intent(out)              :: info
    !

    type(psb_dspmat_type)       :: atmp
    type(psb_d_csc_sparse_mat)  :: acsc, acred
    type(psb_d_csr_sparse_mat)  :: acsr
    integer :: icnt,nlp,k,n,ia,isz,naggr,i,j,m, nz, i1, i2, ns, irdp
    integer :: num_agg, num_unagg
    real(psb_dpk_)  :: cpling, tcl
    logical :: recovery
    integer :: debug_level, debug_unit
    integer :: ictxt,np,me,err_act
    integer :: nrow, ncol, n_ne, nnz, job, nr, nc, num, ir, ic, ip, nr2, nc2, nz2
    integer           :: irmax, icmax, idx, mxk2
    real(psb_dpk_)    :: rmax, cmax
    character(len=20) :: name, ch_err
    integer :: icntl(10), infov(10)
    real(psb_dpk_) :: dcntl(10) 
    integer :: liw, ldw
    integer, allocatable        :: iwork(:), perm(:), perm2(:),mark2(:)
    integer, allocatable        :: ind_agg(:), ind_unagg(:)
    real(psb_dpk_), allocatable :: dwork(:)

    if (psb_get_errstatus() /= 0) return 
    info = psb_success_
    name = 'mld_d_match'
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    nr = a%get_nrows()
    !
    ! Purely local aggregation, clip away halo
    ! Also, take out diagonal for sake of mc64
    !
    call a%csclip(atmp,info,jmax=nr)
    call atmp%clip_diag(info)
    call atmp%mv_to(acsc)  
    nr  = acsc%get_nrows()
    nc  = acsc%get_ncols()
    nnz = acsc%get_nzeros()

    call mc64id(icntl,dcntl)
    icntl(1) = -1
    icntl(2) = -1
    !
    ! Job default. Should add an option somewhere. 
    ! 
    job      = 5  
    select case(job)
    case (5) 
      liw = 3*nc+2*nr
      ldw = nc+2*nr+nnz
    case default
      ! Unknown job, set to 5 
      job = 5
      liw = 3*nc+2*nr
      ldw = nc+2*nr+nnz
    end select

    allocate(iwork(liw),dwork(ldw),perm(nr),mark(nr),&
         & perm2(nr),mark2(nr),stat=info)
    if (info /= 0) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/ldw,0,0,0,0/),&
           & a_err='real(psb_dpk_)')
      goto 9999
    end if
    
    perm  = 0
    perm2 = 0 
    ! First aggregation

    call mc64ad(job,nr,nc,nnz,acsc%icp,acsc%ia,acsc%val,&
         & num,perm,liw,iwork,ldw,dwork,icntl,dcntl,infov)

    if (infov(1) < 0) then 
      write(psb_err_unit,*) 'MC64 returned a failure ',infov(1:2)
      info=psb_err_internal_error_
      call psb_errpush(info,name)
      goto 9999
    end if

   
    ! Mark nodes
    call marking(nr,perm,mark,info) 

    if (info < 0) then 
      write(psb_err_unit,*) 'marking returned a failure ',info
      info=psb_err_internal_error_
      call psb_errpush(info,name)
      goto 9999
    end if


    nagg1     = 0
    num_agg   = 0
    num_unagg = 0
    do i=1, nr
      nagg1 = max(nagg1,mark(i))
      if (mark(i) > 0) num_agg   = num_agg + 1
      if (mark(i) < 0) num_unagg = num_unagg + 1
    end do
    allocate(ind_agg(num_agg), ind_unagg(num_unagg), stat=info)
    if (info /= 0) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/num_agg+num_unagg,0,0,0,0/),&
           & a_err='integer')
      goto 9999
    end if
    i1 = 0
    i2 = 0
    nagg1 = 0
    do i=1, nr
      if (mark(i) > 0) then 
        nagg1 = max(nagg1,mark(i))
        i1 = i1 + 1
        ind_agg(i1) = i 
      end if
      if (mark(i) < 0) then 
        i2 = i2 + 1 
        ind_unagg(i2) = i
      end if
    end do

    !
    ! Build reduced matrix and perform second aggregation step 
    !
    if (num_agg < nr) then 
      nr2 = num_unagg
      nc2 = num_unagg
      call acred%allocate(nr2,nc2,nnz)
      acred%icp(1) = 1
      irdp         = 1
      do i=1, nc2
        ic = ind_unagg(i)
        j  = acsc%icp(ic)
        nz = acsc%icp(ic+1)-j
        do k=1, nz 
          ir = acsc%ia(j)
          if (ir /= ic) then 
            ip = psb_ibsrch(ir,num_unagg,ind_unagg) 
            if (ip > 0) then 
              acred%ia(irdp)  = ip
              acred%val(irdp) = acsc%val(j)
              irdp = irdp + 1 
            end if
          end if
          j = j + 1 
        end do
        acred%icp(i+1) = irdp
      end do
      nz2 = irdp - 1
      call mc64ad(job,nr2,nc2,nz2,acred%icp,acred%ia,acred%val,&
           & num,perm2,liw,iwork,ldw,dwork,icntl,dcntl,infov)

      call marking(nc2,perm2,mark2,info) 
      if (info < 0) then 
        write(psb_err_unit,*) 'marking returned a failure ',info
        info=psb_err_internal_error_
        call psb_errpush(info,name)
        goto 9999
      end if

      !
      ! Update marks and rebuild list of unaggregated indices
      !
      num_unagg = 0
      mxk2      = 0
      do i=1, nc2 
        if (mark2(i) > 0) then 
          mark(ind_unagg(i))   = nagg1+mark2(i)
          mxk2 = max(mxk2,mark2(i))
        else
          mark(ind_unagg(i))   = mark2(i)
          num_unagg            = num_unagg + 1 
          ind_unagg(num_unagg) = ind_unagg(i) 
        end if
      end do
      call acred%free()
      nagg1 = nagg1 + mxk2 
    end if
    !
    ! Fixup leftover unaggregated nodes, if any
    !

    !
    ! A CSR copy comes in handy. Can we do it
    ! with less memory? 
    !
    call acsc%cp_to_fmt(acsr,info)
    do i=1, num_unagg 
      ! Look at row/col neighbours of row unagg(i)
      k = ind_unagg(i)
      ! Max along the row
      rmax  = dzero
      irmax = -1
      do j=acsr%irp(k),acsr%irp(k+1)-1
        if (abs(acsr%val(j)) > rmax) then 
          rmax  = abs(acsr%val(j))
          irmax = acsr%ja(j)
        end if
      end do
      ! Max along the column
      cmax  = dzero
      icmax = -1
      do j=acsc%icp(k),acsc%icp(k+1)-1
        if (abs(acsc%val(j)) > cmax) then 
          cmax  = abs(acsc%val(j))
          icmax = acsc%ia(j)
        end if
      end do
      if (rmax > cmax) then 
        idx = irmax
      else
        idx = icmax
      end if
      if (mark(idx) < 0) then 
        nagg1 = nagg1 + 1 
        mark(idx) = nagg1
        mark(k)   = nagg1
      else
        mark(k)   = mark(idx)
      end if
    end do

    call acsr%free()
    call acsc%free()

    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return


  end subroutine mld_d_match

#endif  
end module mld_mc64_tools_mod
