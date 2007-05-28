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
subroutine mld_zsludist_bld(a,desc_a,p,info)
  use psb_base_mod
  use psb_prec_mod, mld_protect_name => mld_zsludist_bld

  implicit none 

  type(psb_zspmat_type), intent(inout)      :: a
  type(psb_desc_type), intent(in)        :: desc_a
  type(psb_zbaseprc_type), intent(inout) :: p
  integer, intent(out)                   :: info

  integer            :: i,j,nza,nzb,nzt,ictxt,me,np,err_act,&
       &                mglob,ifrst,ibcheck,nrow,ncol,npr,npc, ip
  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='mld_zslu_bld'
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)

  if (toupper(a%fida) /= 'CSR') then 
    write(0,*) 'Unimplemented input to SLU_BLD'
    goto 9999
  endif


  !
  ! WARN we need to check for a BLOCK distribution. 
  !
  nrow = psb_cd_get_local_rows(desc_a)
  ncol = psb_cd_get_local_cols(desc_a)
  ifrst   = desc_a%loc_to_glob(1) 
  ibcheck = desc_a%loc_to_glob(nrow) - ifrst + 1 
  ibcheck = ibcheck - nrow
  call psb_amx(ictxt,ibcheck)
  if (ibcheck > 0) then 
    write(0,*) 'Warning: does not look like a BLOCK distribution'
  endif

  mglob = psb_cd_get_global_rows(desc_a)
  nzt   = psb_sp_get_nnzeros(a)

  call psb_loc_to_glob(a%ia1(1:nzt),desc_a,info,iact='I')
  
  npr = np
  npc = 1
  ip = floor(sqrt(dble(np)))
  do 
    if (ip <= 1) exit
    if (mod(np,ip)==0) then 
      npr = np/ip
      npc = ip
      exit
    end if
    ip = ip - 1
  end do
!!$  write(0,*) 'Process grid : ',npr,npc
  call mld_zsludist_factor(mglob,nrow,nzt,ifrst,&
       & a%aspk,a%ia2,a%ia1,p%iprcparm(slud_ptr_),&
       & npr, npc, info)
  if (info /= 0) then
    ch_err='psb_slud_fact'
    call psb_errpush(4110,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
    goto 9999
  end if
  
  call psb_glob_to_loc(a%ia1(1:nzt),desc_a,info,iact='I')

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zsludist_bld

