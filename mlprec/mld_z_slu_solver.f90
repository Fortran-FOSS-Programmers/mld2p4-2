!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010, 2010
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
!
!
!
!
!
!

module mld_z_slu_solver

  use iso_c_binding
  use mld_z_prec_type

  type, extends(mld_z_base_solver_type) :: mld_z_slu_solver_type
    type(c_ptr)                 :: lufactors=c_null_ptr
    integer(c_long_long)        :: symbsize=0, numsize=0
  contains
    procedure, pass(sv) :: build   => z_slu_solver_bld
    procedure, pass(sv) :: apply_a => z_slu_solver_apply
    procedure, pass(sv) :: free    => z_slu_solver_free
    procedure, pass(sv) :: seti    => z_slu_solver_seti
    procedure, pass(sv) :: setc    => z_slu_solver_setc
    procedure, pass(sv) :: setr    => z_slu_solver_setr
    procedure, pass(sv) :: descr   => z_slu_solver_descr
    procedure, pass(sv) :: sizeof  => z_slu_solver_sizeof
  end type mld_z_slu_solver_type


  private :: z_slu_solver_bld, z_slu_solver_apply, &
       &  z_slu_solver_free,   z_slu_solver_seti, &
       &  z_slu_solver_setc,   z_slu_solver_setr,&
       &  z_slu_solver_descr,  z_slu_solver_sizeof


  interface 
    function mld_zslu_fact(n,nnz,values,rowptr,colind,&
         & lufactors)&
         & bind(c,name='mld_zslu_fact') result(info)
      use iso_c_binding
      integer(c_int), value :: n,nnz
      integer(c_int)        :: info
      !integer(c_long_long)  :: ssize, nsize
      integer(c_int)        :: rowptr(*),colind(*)
      complex(c_double_complex)        :: values(*)
      type(c_ptr)           :: lufactors
    end function mld_zslu_fact
  end interface

  interface 
    function mld_zslu_solve(itrans,n,x, b, ldb, lufactors)&
         & bind(c,name='mld_zslu_solve') result(info)
      use iso_c_binding
      integer(c_int)        :: info
      integer(c_int), value :: itrans,n,ldb
      complex(c_double_complex)        :: x(*), b(ldb,*)
      type(c_ptr), value    :: lufactors
    end function mld_zslu_solve
  end interface

  interface 
    function mld_zslu_free(lufactors)&
         & bind(c,name='mld_zslu_free') result(info)
      use iso_c_binding
      integer(c_int)        :: info
      type(c_ptr), value    :: lufactors
    end function mld_zslu_free
  end interface

contains

  subroutine z_slu_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_z_slu_solver_type), intent(in) :: sv
    complex(psb_dpk_),intent(inout)         :: x(:)
    complex(psb_dpk_),intent(inout)         :: y(:)
    complex(psb_dpk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)          :: trans
    complex(psb_dpk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info

    integer    :: n_row,n_col
    complex(psb_dpk_), pointer :: ww(:)
    integer    :: ictxt,np,me,i, err_act
    character          :: trans_
    character(len=20)  :: name='z_slu_solver_apply'

    call psb_erractionsave(err_act)

    info = psb_success_

    trans_ = psb_toupper(trans)
    select case(trans_)
    case('N')
    case('T','C')
    case default
      call psb_errpush(psb_err_iarg_invalid_i_,name)
      goto 9999
    end select

    n_row = desc_data%get_local_rows()
    n_col = desc_data%get_local_cols()

    if (n_col <= size(work)) then 
      ww => work(1:n_col)
    else
      allocate(ww(n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/n_col,0,0,0,0/),&
             & a_err='complex(psb_dpk_)')
        goto 9999      
      end if
    endif

    select case(trans_)
    case('N')
      info = mld_zslu_solve(0,n_row,ww,x,n_row,sv%lufactors)
    case('T')
      info = mld_zslu_solve(1,n_row,ww,x,n_row,sv%lufactors)
    case('C')
      info = mld_zslu_solve(2,n_row,ww,x,n_row,sv%lufactors)
    case default
      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in ILU subsolve')
      goto 9999
    end select

    if (info == psb_success_) call psb_geaxpby(alpha,ww,beta,y,desc_data,info)


    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,a_err='Error in subsolve')
      goto 9999
    endif

    if (n_col > size(work)) then 
      deallocate(ww)
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine z_slu_solver_apply

  subroutine z_slu_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)

    use psb_base_mod

    Implicit None

    ! Arguments
    type(psb_zspmat_type), intent(in), target           :: a
    Type(psb_desc_type), Intent(in)                     :: desc_a 
    class(mld_z_slu_solver_type), intent(inout)         :: sv
    character, intent(in)                               :: upd
    integer, intent(out)                                :: info
    type(psb_zspmat_type), intent(in), target, optional :: b
    class(psb_z_base_sparse_mat), intent(in), optional  :: amold
    class(psb_z_base_vect_type), intent(in), optional   :: vmold
    ! Local variables
    type(psb_zspmat_type) :: atmp
    type(psb_z_csr_sparse_mat) :: acsr
    integer :: n_row,n_col, nrow_a, nztota
    integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
    character(len=20)  :: name='z_slu_solver_bld', ch_err
    
    info=psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = desc_a%get_context()
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'


    n_row  = desc_a%get_local_rows()
    n_col  = desc_a%get_local_cols()

    if (psb_toupper(upd) == 'F') then 

      call a%cscnv(atmp,info,type='coo')
      call psb_rwextd(n_row,atmp,info,b=b) 
      call atmp%cscnv(info,type='csr',dupl=psb_dupl_add_)
      call atmp%mv_to(acsr)
      nrow_a = acsr%get_nrows()
      nztota = acsr%get_nzeros()
      ! Fix the entres to call C-base SuperLU
      acsr%ja(:)  = acsr%ja(:)  - 1
      acsr%irp(:) = acsr%irp(:) - 1
      info = mld_zslu_fact(nrow_a,nztota,acsr%val,&
           & acsr%irp,acsr%ja,sv%lufactors)

      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='mld_zslu_fact'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      call acsr%free()
      call atmp%free()
    else
      ! ? 
        info=psb_err_internal_error_
        call psb_errpush(info,name)
        goto 9999
      
    end if

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' end'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_slu_solver_bld


  subroutine z_slu_solver_seti(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_z_slu_solver_type), intent(inout) :: sv 
    integer, intent(in)                    :: what 
    integer, intent(in)                    :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='z_slu_solver_seti'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(what) 
    case default
!!$      write(0,*) name,': Error: invalid WHAT'
!!$      info = -2
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_slu_solver_seti

  subroutine z_slu_solver_setc(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_z_slu_solver_type), intent(inout) :: sv
    integer, intent(in)                    :: what 
    character(len=*), intent(in)           :: val
    integer, intent(out)                   :: info
    Integer :: err_act, ival
    character(len=20)  :: name='z_slu_solver_setc'

    info = psb_success_
    call psb_erractionsave(err_act)


    call mld_stringval(val,ival,info)
    if (info == psb_success_) call sv%set(what,ival,info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info, name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_slu_solver_setc
  
  subroutine z_slu_solver_setr(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_z_slu_solver_type), intent(inout) :: sv 
    integer, intent(in)                    :: what 
    real(psb_dpk_), intent(in)             :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='z_slu_solver_setr'

    call psb_erractionsave(err_act)
    info = psb_success_

    select case(what)
    case default
!!$      write(0,*) name,': Error: invalid WHAT'
!!$      info = -2
!!$      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_slu_solver_setr

  subroutine z_slu_solver_free(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_z_slu_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='z_slu_solver_free'

    call psb_erractionsave(err_act)

    
    info = mld_zslu_free(sv%lufactors)
    
    if (info /= psb_success_) goto 9999
    sv%lufactors = c_null_ptr


    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_slu_solver_free

  subroutine z_slu_solver_descr(sv,info,iout,coarse)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_z_slu_solver_type), intent(in) :: sv
    integer, intent(out)                     :: info
    integer, intent(in), optional            :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_z_slu_solver_descr'
    integer :: iout_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif
    
    write(iout_,*) '  SuperLU Sparse Factorization Solver. '

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_slu_solver_descr

  function z_slu_solver_sizeof(sv) result(val)
    use psb_base_mod
    implicit none 
    ! Arguments
    class(mld_z_slu_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i

    val = 2*psb_sizeof_int + psb_sizeof_dp
    val = val + sv%symbsize
    val = val + sv%numsize
    return
  end function z_slu_solver_sizeof

end module mld_z_slu_solver
