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
subroutine mld_zprecinit(p,ptype,info,nlev)

  use psb_base_mod
  use mld_prec_mod, psb_protect_name => mld_zprecinit

  implicit none
  type(mld_zprec_type), intent(inout)    :: p
  character(len=*), intent(in)           :: ptype
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: nlev

  integer                                :: nlev_, ilev_

  info = 0
  
  if (allocated(p%baseprecv)) then 
    call mld_precfree(p,info) 
    if (info /=0) then 
      ! Do we want to do something? 
    endif
  endif

  select case(toupper(ptype(1:len_trim(ptype))))
  case ('NONE','NOPREC') 
    nlev_ = 1
    ilev_ = 1
    allocate(p%baseprecv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(ifpsz,p%baseprecv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(dfpsz,p%baseprecv(ilev_)%dprcparm,info)
    if (info /= 0) return
    p%baseprecv(ilev_)%iprcparm(:) = 0
    p%baseprecv(ilev_)%iprcparm(prec_type_)     = noprec_
    p%baseprecv(ilev_)%iprcparm(sub_solve_)     = f_none_
    p%baseprecv(ilev_)%iprcparm(sub_restr_)      = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_prol_)       = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_ren_)       = 0
    p%baseprecv(ilev_)%iprcparm(n_ovr_)      = 0
    p%baseprecv(ilev_)%iprcparm(smooth_sweeps_) = 1

  case ('DIAG')
    nlev_ = 1
    ilev_ = 1
    allocate(p%baseprecv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(ifpsz,p%baseprecv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(dfpsz,p%baseprecv(ilev_)%dprcparm,info)
    if (info /= 0) return
    p%baseprecv(ilev_)%iprcparm(prec_type_)     = diag_
    p%baseprecv(ilev_)%iprcparm(sub_solve_)     = f_none_
    p%baseprecv(ilev_)%iprcparm(sub_restr_)      = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_prol_)       = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_ren_)       = 0 
    p%baseprecv(ilev_)%iprcparm(n_ovr_)      = 0
    p%baseprecv(ilev_)%iprcparm(smooth_sweeps_) = 1

  case ('BJAC') 
    nlev_ = 1
    ilev_ = 1
    allocate(p%baseprecv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(ifpsz,p%baseprecv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(dfpsz,p%baseprecv(ilev_)%dprcparm,info)
    if (info /= 0) return
    p%baseprecv(ilev_)%iprcparm(:)            = 0
    p%baseprecv(ilev_)%iprcparm(prec_type_)      = bjac_
    p%baseprecv(ilev_)%iprcparm(sub_solve_)      = ilu_n_
    p%baseprecv(ilev_)%iprcparm(sub_restr_)       = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_prol_)        = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_ren_)        = 0
    p%baseprecv(ilev_)%iprcparm(n_ovr_)       = 0
    p%baseprecv(ilev_)%iprcparm(sub_fill_in_) = 0
    p%baseprecv(ilev_)%iprcparm(smooth_sweeps_)  = 1

  case ('ASM','AS')
    nlev_ = 1
    ilev_ = 1
    allocate(p%baseprecv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(ifpsz,p%baseprecv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(dfpsz,p%baseprecv(ilev_)%dprcparm,info)
    if (info /= 0) return
    p%baseprecv(ilev_)%iprcparm(prec_type_)      = as_ 
    p%baseprecv(ilev_)%iprcparm(sub_solve_)      = ilu_n_
    p%baseprecv(ilev_)%iprcparm(sub_restr_)       = psb_halo_
    p%baseprecv(ilev_)%iprcparm(sub_prol_)        = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_ren_)        = 0
    p%baseprecv(ilev_)%iprcparm(n_ovr_)       = 1
    p%baseprecv(ilev_)%iprcparm(sub_fill_in_) = 0
    p%baseprecv(ilev_)%iprcparm(smooth_sweeps_)  = 1


  case ('MLD', 'ML')
    
    if (present(nlev)) then 
      nlev_ = max(1,nlev)
    else
      nlev_ = 2
    end if
    if (nlev_ == 1) then 
      write(0,*) 'Warning: requested ML preconditioner with NLEV=1'
    endif
    ilev_ = 1
    allocate(p%baseprecv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(ifpsz,p%baseprecv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(dfpsz,p%baseprecv(ilev_)%dprcparm,info)
    if (info /= 0) return
    p%baseprecv(ilev_)%iprcparm(prec_type_)      = as_ 
    p%baseprecv(ilev_)%iprcparm(sub_solve_)      = ilu_n_
    p%baseprecv(ilev_)%iprcparm(sub_restr_)       = psb_halo_
    p%baseprecv(ilev_)%iprcparm(sub_prol_)        = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_ren_)        = 0
    p%baseprecv(ilev_)%iprcparm(n_ovr_)       = 1
    p%baseprecv(ilev_)%iprcparm(sub_fill_in_) = 0
    p%baseprecv(ilev_)%iprcparm(smooth_sweeps_)  = 1
    if (nlev_ == 1) return 

    do ilev_ = 2, nlev_ -1 
      if (info == 0) call psb_realloc(ifpsz,p%baseprecv(ilev_)%iprcparm,info)
      if (info == 0) call psb_realloc(dfpsz,p%baseprecv(ilev_)%dprcparm,info)
      if (info /= 0) return
      p%baseprecv(ilev_)%iprcparm(prec_type_)       = bjac_
      p%baseprecv(ilev_)%iprcparm(sub_restr_)        = psb_none_
      p%baseprecv(ilev_)%iprcparm(sub_prol_)         = psb_none_
      p%baseprecv(ilev_)%iprcparm(sub_ren_)         = 0
      p%baseprecv(ilev_)%iprcparm(n_ovr_)        = 0
      p%baseprecv(ilev_)%iprcparm(ml_type_)      = mult_ml
      p%baseprecv(ilev_)%iprcparm(aggr_alg_)     = dec_aggr_
      p%baseprecv(ilev_)%iprcparm(aggr_kind_)    = tent_prol_
      p%baseprecv(ilev_)%iprcparm(coarse_mat_)   = distr_mat_
      p%baseprecv(ilev_)%iprcparm(smooth_pos_)     = post_smooth_
      p%baseprecv(ilev_)%iprcparm(aggr_eig_)    = max_norm_
      p%baseprecv(ilev_)%iprcparm(sub_solve_)       = ilu_n_
      p%baseprecv(ilev_)%iprcparm(sub_fill_in_)  = 0
      p%baseprecv(ilev_)%iprcparm(smooth_sweeps_)   = 1
      p%baseprecv(ilev_)%dprcparm(aggr_damp_) = 4.d0/3.d0         
    end do
    ilev_ = nlev_
    if (info == 0) call psb_realloc(ifpsz,p%baseprecv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(dfpsz,p%baseprecv(ilev_)%dprcparm,info)
    if (info /= 0) return
    p%baseprecv(ilev_)%iprcparm(prec_type_)       = bjac_
    p%baseprecv(ilev_)%iprcparm(sub_restr_)        = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_prol_)         = psb_none_
    p%baseprecv(ilev_)%iprcparm(sub_ren_)         = 0
    p%baseprecv(ilev_)%iprcparm(n_ovr_)        = 0
    p%baseprecv(ilev_)%iprcparm(ml_type_)      = mult_ml
    p%baseprecv(ilev_)%iprcparm(aggr_alg_)     = dec_aggr_
    p%baseprecv(ilev_)%iprcparm(aggr_kind_)    = tent_prol_
    p%baseprecv(ilev_)%iprcparm(coarse_mat_)   = distr_mat_
    p%baseprecv(ilev_)%iprcparm(smooth_pos_)     = post_smooth_
    p%baseprecv(ilev_)%iprcparm(aggr_eig_)    = max_norm_
    p%baseprecv(ilev_)%iprcparm(sub_solve_)       = umf_
    p%baseprecv(ilev_)%iprcparm(sub_fill_in_)  = 0
    p%baseprecv(ilev_)%iprcparm(smooth_sweeps_)   = 4
    p%baseprecv(ilev_)%dprcparm(aggr_damp_) = 4.d0/3.d0         

  case default
    write(0,*) 'Unknown preconditioner type request "',ptype,'"'
    info = 2

  end select


end subroutine mld_zprecinit
