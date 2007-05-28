!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
!*****************************************************************************
!*                                                                           *
!* This is where the action takes place.                                     *
!* ASMATBLD does the setup: building the prec descriptor plus retrieving     *
!*                           matrix rows if needed                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!* some open code does the renumbering                                       *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine psb_zilu_bld(a,desc_a,p,upd,info,blck)
  use psb_base_mod
  use psb_prec_mod, mld_protect_name => psb_zilu_bld

  implicit none
  !                                                                               
  !     .. Scalar Arguments ..                                                    
  integer, intent(out)                      :: info
  !     .. array Arguments ..                                                     
  type(psb_zspmat_type), intent(in), target :: a
  type(psb_zbaseprc_type), intent(inout)    :: p
  type(psb_desc_type), intent(in)           :: desc_a
  character, intent(in)                     :: upd
  type(psb_zspmat_type), intent(in), optional :: blck
  !     .. Local Scalars ..                                                       
  integer  ::    i, j, jj, k, kk, m
  integer  ::    int_err(5)
  character ::        trans, unitd
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6, t7, t8
  logical, parameter :: debugprt=.false., debug=.false., aggr_dump=.false.
  integer   nztota, nztotb, nztmp, nzl, nnr, ir, err_act,&
       & n_row, nrow_a, ind, iind, i1,i2,ia
  integer :: ictxt,np,me
  character(len=20)      :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_zilu_bld'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  trans = 'N'
  unitd = 'U'

  if (allocated(p%av)) then 
    if (size(p%av) < bp_ilu_avsz) then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)
    endif
  end if
  if (.not.allocated(p%av)) then 
    allocate(p%av(max_avsz),stat=info)
    if (info /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    end if
  endif
!!$  call psb_csprt(50+me,a,head='% (A)')    

  nrow_a = psb_cd_get_local_rows(desc_a)
  nztota = psb_sp_get_nnzeros(a)
  if (present(blck)) then 
    nztota = nztota + psb_sp_get_nnzeros(blck)
  end if
  if (debug) write(0,*)me,': out get_nnzeros',nztota,a%m,a%k
  if (debug) call psb_barrier(ictxt)

  n_row  = p%desc_data%matrix_data(psb_n_row_)
  p%av(l_pr_)%m  = n_row
  p%av(l_pr_)%k  = n_row
  p%av(u_pr_)%m  = n_row
  p%av(u_pr_)%k  = n_row
  call psb_sp_all(n_row,n_row,p%av(l_pr_),nztota,info)
  if (info == 0) call psb_sp_all(n_row,n_row,p%av(u_pr_),nztota,info)
  if(info/=0) then
    info=4010
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (allocated(p%d)) then 
    if (size(p%d) < n_row) then 
      deallocate(p%d)
    endif
  endif
  if (.not.allocated(p%d)) then 
    allocate(p%d(n_row),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

  endif


  !
  ! Ok, factor the matrix.  
  !
  t5 = psb_wtime()
  call psb_ilu_fct(a,p%av(l_pr_),p%av(u_pr_),p%d,info,blck=blck)
  if(info/=0) then
    info=4010
    ch_err='psb_ilu_fct'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  if (debugprt) then 
    !
    ! Print out the factors on file.
    !
    open(80+me)

    call psb_csprt(80+me,p%av(l_pr_),head='% Local L factor')
    write(80+me,*) '% Diagonal: ',p%av(l_pr_)%m
    do i=1,p%av(l_pr_)%m
      write(80+me,*) i,i,p%d(i)
    enddo
    call psb_csprt(80+me,p%av(u_pr_),head='% Local U factor')

    close(80+me)
  endif

!!$  call psb_csprt(60+me,a,head='% (A)')    


  !    ierr = MPE_Log_event( ifcte, 0, "st SIMPLE" )
  t6 = psb_wtime()
  !
  !    write(0,'(i3,1x,a,3(1x,g18.9))') me,'renum/factor time',t3-t2,t6-t5
  !    if (me==0) write(0,'(a,3(1x,g18.9))') 'renum/factor time',t3-t2,t6-t5

  if (psb_sp_getifld(psb_upd_,p%av(u_pr_),info) /= psb_upd_perm_) then
    call psb_sp_trim(p%av(u_pr_),info)
  endif

  if (psb_sp_getifld(psb_upd_,p%av(l_pr_),info) /= psb_upd_perm_) then
    call psb_sp_trim(p%av(l_pr_),info)
  endif


  if (debug) write(0,*) me,'End of ilu_bld'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_zilu_bld

