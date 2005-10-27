! File:  psb_dipcoo2csc.f90 
! Subroutine: 
! Parameters:

subroutine psb_dipcoo2csc(a,info,clshr)
  use psb_spmat_type
  use psb_const_mod
  use psb_serial_mod, only : psb_fixcoo
  use psb_error_mod
  implicit none

  !....Parameters...
  Type(psb_dspmat_type), intent(inout) :: A
  Integer, intent(out)                 :: info
  logical, optional                    :: clshr

  integer, pointer    :: iaux(:), itemp(:)
  !locals
  logical             :: clshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act,nc,icl
  Integer, Parameter  :: maxtry=8
  logical, parameter  :: debug=.false.
  character(len=20)                 :: name, ch_err

  name='psb_ipcoo2csc'
  info  = 0
  call psb_erractionsave(err_act)

  if(debug) write(0,*)'Inside ipcoo2csc',a%fida,a%m
  if (a%fida /= 'COO') then 
    write(0,*) 'ipcoo2csc Invalid input ',a%fida
    info = -1
    call psb_errpush(info,name)
    goto 9999
  end if
  if (present(clshr)) then 
    clshr_ = clshr
  else
    clshr_ = .false.
  end if

  call psb_fixcoo(a,info,idir=1)
  nc  = a%k
  nza = a%infoa(psb_nnz_)
  allocate(iaux(nc+1))
  if(debug) write(0,*)'DIPCOO2CSC: out of fixcoo',nza,nc,size(a%ia2),size(iaux)

  itemp => a%ia1
  a%ia1 => a%ia2
  a%ia2 => iaux

  !
  ! This routine can be used in two modes:
  ! 1. Normal: just look at the col indices and trust them. This
  !    implies putting in empty cols where needed. In this case you
  !    can get in trouble if A%M < A%IA1(NZA)
  ! 2. Shrink mode: disregard the actual value of the col indices,
  !    just treat them as ident markers. In this case you can get in
  !    trouble when the number of distinct col indices is greater 
  !    than A%M
  !
  !

  a%ia2(1) = 1

  if (nza <= 0) then 
    do i=1,nc
      a%ia2(i+1) = a%ia2(i)
    end do
  else

    if (clshr_) then 
      
      
      j = 1
      i = 1
      icl = itemp(j) 

      do j=1, nza
        if (itemp(j) /= icl) then 
          a%ia2(i+1) = j
          icl = itemp(j) 
          i = i + 1
          if (i>nc) then 
            write(0,*) 'IPCOO2CSC: CLSHR=.true. : ',&
             & i, nc,' Expect trouble!'
            exit
          end if
        endif
      enddo 
!      write(0,*) 'Exit from loop',j,nza,i
      do 
        if (i>=nc+1) exit
        a%ia2(i+1) = j
        i = i + 1
      end do

    else
      
      if (nc < itemp(nza)) then 
        write(0,*) 'IPCOO2CSC: CLSHR=.false. : ',&
             &nc,itemp(nza),' Expect trouble!'
      end if
             

      j = 1 
      i = 1
      icl = itemp(j) 

      outer: do 
        inner: do 
          if (i >= icl) exit inner
          if (i>nc) then 
            write(0,*) 'strange situation: i>nc ',i,nc,j,nza,icl,idl
            exit outer
          end if
          a%ia2(i+1) = a%ia2(i) 
          i = i + 1
        end do inner
        j = j + 1
        if (j > nza) exit
        if (itemp(j) /= icl) then 
          a%ia2(i+1) = j
          icl = itemp(j) 
          i = i + 1
        endif
      enddo outer
      !
      ! Cleanup empty cols at the end
      !
      if (j /= (nza+1)) then 
        write(0,*) 'IPCOO2CSC : Problem from loop :',j,nza
      endif
      do 
        if (i>=nc+1) exit
        a%ia2(i+1) = j
        i = i + 1
      end do

    endif

  end if

!!$  write(0,*) 'IPcoo2csc end loop ',i,nc,a%ia2(nc+1),nza
  a%fida='CSC'

  deallocate(itemp)
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return

end Subroutine psb_dipcoo2csc