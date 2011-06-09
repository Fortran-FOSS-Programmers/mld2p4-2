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
      if (perm(i) > 0) then 
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
        if (perm(i) < 0) then 
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
#endif  
end module mld_mc64_tools_mod
