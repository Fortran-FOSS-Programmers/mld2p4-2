!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 
!  
!        Salvatore Filippone    Cranfield University, UK
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta, IT
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
subroutine mld_s_base_onelev_check(lv,info)
  
  use psb_base_mod
  use mld_s_onelev_mod, mld_protect_name => mld_s_base_onelev_check

  Implicit None

  ! Arguments
  class(mld_s_onelev_type), intent(inout) :: lv 
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_)           :: err_act
  character(len=20) :: name='s_base_onelev_check'

  call psb_erractionsave(err_act)
  info = psb_success_

  call mld_check_def(lv%parms%sweeps_pre,&
       & 'Jacobi sweeps',ione,is_int_non_negative)
  call mld_check_def(lv%parms%sweeps_post,&
       & 'Jacobi sweeps',ione,is_int_non_negative)

  if (allocated(lv%sm)) then 
    call lv%sm%check(info)
  else 
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  if (allocated(lv%sm2a)) then 
    call lv%sm2a%check(info)
  else if (.not.inner_check(lv%sm2,lv%sm)) then 
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains
  function inner_check(smp,sm) result(res)
    implicit none
    logical :: res
    class(mld_s_base_smoother_type), intent(in), pointer :: smp 
    class(mld_s_base_smoother_type), intent(in), target  :: sm

    res = associated(smp, sm)
  end function inner_check
  
end subroutine mld_s_base_onelev_check
