!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
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
! File: mld_prec_type.f90
!
! Module: mld_prec_type
!
!  This module defines: 
!  - the mld_prec_type data structure containing the preconditioner and related
!    data structures;
!  - integer constants defining the preconditioner;
!  - character constants describing the preconditioner (used by the routines
!    printing out a preconditioner description);
!  - the interfaces to the routines for the management of the preconditioner
!    data structure (see below).
!
!  It contains routines for
!  - converting character constants defining the preconditioner into integer
!    constants; 
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_d_prec_type

  use mld_base_prec_type
  !
  ! Type: mld_Tprec_type.
  !
  !  It is the data type containing all the information about the multilevel
  !  preconditioner (here and in the following 'T' denotes 'd', 's', 'c' and
  !  'z', according to the real/complex, single/double precision version of
  !  MLD2P4). It consists of an array of 'one-level' intermediate data structures
  !  of type mld_Tonelev_type, each containing the information needed to apply
  !  the smoothing and the coarse-space correction at a generic level.
  !
  !  type mld_Tprec_type
  !    type(mld_Tonelev_type), allocatable :: precv(:) 
  !  end type mld_Tprec_type
  ! 
  !  Note that the levels are numbered in increasing order starting from
  !  the finest one and the number of levels is given by size(precv(:)).
  !
  !
  ! Type: mld_Tonelev_type.
  !
  !  It is the data type containing the necessary items for the	current
  !  level (essentially, the base preconditioner, the current-level	matrix
  !  and the restriction and prolongation operators).
  !
  !  type mld_Tonelev_type
  !    type(mld_Tbaseprec_type)       :: prec
  !    integer, allocatable           :: iprcparm(:)
  !    real(psb_Tpk_), allocatable    :: rprcparm(:)
  !    type(psb_Tspmat_type)          :: ac
  !    type(psb_desc_type)            :: desc_ac
  !    type(psb_Tspmat_type), pointer :: base_a    => null()
  !    type(psb_desc_type), pointer   :: base_desc => null()
  !    type(psb_Tlinmap_type)         :: map
  !  end type mld_Tonelev_type
  !
  !  Note that psb_Tpk denotes the kind of the real data type to be chosen
  !  according to single/double precision version of MLD2P4.
  !
  !   prec         -  type(mld_Tbaseprec_type). 
  !                   The current level preconditioner (aka smoother).
  !   iprcparm     -  integer, dimension(:), allocatable.
  !                   The integer parameters defining the multilevel strategy.
  !   rprcparm     -  real(psb_Ypk_), dimension(:), allocatable.
  !                   The real parameters defining the multilevel strategy.
  !   ac           -  The local part of the current-level matrix, built by
  !                   coarsening the previous-level matrix.
  !   desc_ac      -  type(psb_desc_type).
  !                   The communication descriptor associated to the matrix
  !                   stored in ac.
  !   base_a       -  type(psb_zspmat_type), pointer.
  !                   Pointer (really a pointer!) to the local part of the current 
  !                   matrix (so we have a unified treatment of residuals).
  !                   We need this to avoid passing explicitly the current matrix
  !                   to the routine which applies the preconditioner.
  !   base_desc    -  type(psb_desc_type), pointer.
  !                   Pointer to the communication descriptor associated to the
  !                   matrix pointed by base_a.
  !   map          -  Stores the maps (restriction and prolongation) between the
  !                   vector spaces associated to the index spaces of the previous
  !                   and current levels.
  !
  ! 
  ! Type: mld_Tbaseprec_type.
  ! 
  !  It holds the smoother (base preconditioner) at a single level.
  !
  !  type mld_Tbaseprec_type
  !    type(psb_Tspmat_type), allocatable :: av(:)
  !    IntrType(psb_Tpk_), allocatable    :: d(:)
  !    type(psb_desc_type)                :: desc_data
  !    integer, allocatable               :: iprcparm(:)
  !    real(psb_Tpk_), allocatable        :: rprcparm(:)
  !    integer, allocatable               :: perm(:),  invperm(:)
  !  end type mld_sbaseprec_type
  !
  !  Note that IntrType denotes the real or complex data type, and psb_Tpk denotes
  !  the kind of the real or complex type, according to the real/complex, single/double
  !  precision version of MLD2P4.
  !
  !    av         -  type(psb_Tspmat_type), dimension(:), allocatable(:).
  !                  The sparse matrices needed to apply the preconditioner at
  !                  the current level ilev. 
  !      av(mld_l_pr_)     -  The L factor of the ILU factorization of the local
  !                           diagonal block of the current-level matrix A(ilev).
  !      av(mld_u_pr_)     -  The U factor of the ILU factorization of the local
  !                           diagonal block of A(ilev), except its diagonal entries
  !                           (stored in d).
  !      av(mld_ap_nd_)    -  The entries of the local part of A(ilev) outside
  !                           the diagonal block, for block-Jacobi sweeps.
  !   d            -  real/complex(psb_Tpk_), dimension(:), allocatable.
  !                   The diagonal entries of the U factor in the ILU factorization
  !                   of A(ilev).
  !   desc_data    -  type(psb_desc_type).
  !                   The communication descriptor associated to the base preconditioner,
  !                   i.e. to the sparse matrices needed to apply the base preconditioner
  !                   at the current level.
  !   iprcparm     -  integer, dimension(:), allocatable.
  !                   The integer parameters defining the base preconditioner K(ilev)
  !                   (the iprcparm entries and values are specified below).
  !   rprcparm     -  real(psb_Tpk_), dimension(:), allocatable.
  !                   The real parameters defining the base preconditioner K(ilev)
  !                   (the rprcparm entries and values are specified below).
  !   perm         -  integer, dimension(:), allocatable.
  !                   The row and column permutations applied to the local part of
  !                   A(ilev) (defined only if iprcparm(mld_sub_ren_)>0). 
  !   invperm      -  integer, dimension(:), allocatable.
  !                   The inverse of the permutation stored in perm.
  !
  !   Note that when the LU factorization of the (local part of the) matrix A(ilev) is
  !   computed instead of the ILU one, by using UMFPACK, SuperLU or SuperLU_dist, the
  !   corresponding L and U factors are stored in data structures provided by those
  !   packages and pointed by prec%iprcparm(mld_umf_ptr), prec%iprcparm(mld_slu_ptr)
  !   or prec%iprcparm(mld_slud_ptr).
  !

  type mld_d_base_solver_type
  contains
    procedure, pass(sv) :: check => d_base_solver_check
    procedure, pass(sv) :: dump  => d_base_solver_dmp
    procedure, pass(sv) :: build => d_base_solver_bld
    procedure, pass(sv) :: apply => d_base_solver_apply
    procedure, pass(sv) :: free  => d_base_solver_free
    procedure, pass(sv) :: seti  => d_base_solver_seti
    procedure, pass(sv) :: setc  => d_base_solver_setc
    procedure, pass(sv) :: setr  => d_base_solver_setr
    generic, public     :: set   => seti, setc, setr
    procedure, pass(sv) :: default => d_base_solver_default
    procedure, pass(sv) :: descr =>   d_base_solver_descr
    procedure, pass(sv) :: sizeof =>  d_base_solver_sizeof
  end type mld_d_base_solver_type

  type  mld_d_base_smoother_type
    class(mld_d_base_solver_type), allocatable :: sv
  contains
    procedure, pass(sm) :: check => d_base_smoother_check
    procedure, pass(sm) :: dump  => d_base_smoother_dmp
    procedure, pass(sm) :: build => d_base_smoother_bld
    procedure, pass(sm) :: apply => d_base_smoother_apply
    procedure, pass(sm) :: free  => d_base_smoother_free
    procedure, pass(sm) :: seti  => d_base_smoother_seti
    procedure, pass(sm) :: setc  => d_base_smoother_setc
    procedure, pass(sm) :: setr  => d_base_smoother_setr
    generic, public     :: set   => seti, setc, setr
    procedure, pass(sm) :: default => d_base_smoother_default
    procedure, pass(sm) :: descr =>   d_base_smoother_descr
    procedure, pass(sm) :: sizeof =>  d_base_smoother_sizeof
  end type mld_d_base_smoother_type

  type, extends(psb_d_base_prec_type)   :: mld_dbaseprec_type
    integer, allocatable                :: iprcparm(:) 
    real(psb_dpk_), allocatable         :: rprcparm(:) 
  end type mld_dbaseprec_type
         
  type mld_donelev_type
    class(mld_d_base_smoother_type), allocatable :: sm
    type(mld_dml_parms)             :: parms 
!!$    type(mld_dbaseprec_type)        :: prec
!!$    integer, allocatable            :: iprcparm(:) 
!!$    real(psb_dpk_), allocatable     :: rprcparm(:) 
    type(psb_dspmat_type)           :: ac
    type(psb_desc_type)             :: desc_ac
    type(psb_dspmat_type), pointer  :: base_a    => null() 
    type(psb_desc_type), pointer    :: base_desc => null() 
    type(psb_dlinmap_type)          :: map
  contains
    procedure, pass(lv) :: descr   => d_base_onelev_descr
    procedure, pass(lv) :: default => d_base_onelev_default
    procedure, pass(lv) :: check => d_base_onelev_check
    procedure, pass(lv) :: dump  => d_base_onelev_dump
    procedure, pass(lv) :: seti  => d_base_onelev_seti
    procedure, pass(lv) :: setr  => d_base_onelev_setr
    procedure, pass(lv) :: setc  => d_base_onelev_setc
    generic, public     :: set   => seti, setr, setc
  end type mld_donelev_type

  type, extends(psb_dprec_type)         :: mld_dprec_type
    integer                             :: ictxt
    type(mld_donelev_type), allocatable :: precv(:) 
  contains
    procedure, pass(prec)               :: d_apply2v => mld_d_apply2v
    procedure, pass(prec)               :: d_apply1v => mld_d_apply1v
    procedure, pass(prec)               :: dump      => mld_d_dump
  end type mld_dprec_type

  private :: d_base_solver_bld,  d_base_solver_apply, &
       &  d_base_solver_free,    d_base_solver_seti, &
       &  d_base_solver_setc,    d_base_solver_setr, &
       &  d_base_solver_descr,   d_base_solver_sizeof, &
       &  d_base_solver_default, d_base_solver_check,&
       &  d_base_solver_dmp, &
       &  d_base_smoother_bld,   d_base_smoother_apply, &
       &  d_base_smoother_free,  d_base_smoother_seti, &
       &  d_base_smoother_setc,  d_base_smoother_setr,&
       &  d_base_smoother_descr, d_base_smoother_sizeof, &
       &  d_base_smoother_default, d_base_smoother_check, &
       &  d_base_smoother_dmp, &
       &  d_base_onelev_seti, d_base_onelev_setc, &
       &  d_base_onelev_setr, d_base_onelev_check, &
       &  d_base_onelev_default, d_base_onelev_dump, &
       &  d_base_onelev_descr


  !
  ! Interfaces to routines for checking the definition of the preconditioner,
  ! for printing its description and for deallocating its data structure
  !

  interface mld_precfree
    module procedure mld_dbase_precfree, mld_d_onelev_precfree, mld_dprec_free
  end interface

  interface mld_nullify_baseprec
    module procedure mld_nullify_dbaseprec
  end interface

  interface mld_nullify_onelevprec
    module procedure  mld_nullify_d_onelevprec
  end interface

  interface mld_precdescr
    module procedure mld_dfile_prec_descr
  end interface

  interface mld_sizeof
    module procedure mld_dprec_sizeof, mld_dbaseprec_sizeof, mld_d_onelev_prec_sizeof
  end interface

  interface mld_precaply
    subroutine mld_dprecaply(prec,x,y,desc_data,info,trans,work)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      import mld_dprec_type
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_dprec_type), intent(in) :: prec
      real(psb_dpk_),intent(in)        :: x(:)
      real(psb_dpk_),intent(inout)     :: y(:)
      integer, intent(out)             :: info
      character(len=1), optional       :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_dprecaply
    subroutine mld_dprecaply1(prec,x,desc_data,info,trans)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      import mld_dprec_type
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_dprec_type), intent(in) :: prec
      real(psb_dpk_),intent(inout)     :: x(:)
      integer, intent(out)             :: info
      character(len=1), optional       :: trans
    end subroutine mld_dprecaply1
  end interface

contains
  !
  ! Function returning the size of the mld_prec_type data structure
  !

  function mld_dprec_sizeof(prec) result(val)
    use psb_sparse_mod
    implicit none 
    type(mld_dprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    val = val + psb_sizeof_int
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
!!$        write(0,*) 'At level ',i, mld_sizeof(prec%precv(i))
        val = val + mld_sizeof(prec%precv(i))
      end do
    end if
  end function mld_dprec_sizeof

  function mld_dbaseprec_sizeof(prec) result(val)
    implicit none 
    type(mld_dbaseprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
      if (prec%iprcparm(mld_prec_status_) == mld_prec_built_) then 
        select case(prec%iprcparm(mld_sub_solve_)) 
        case(mld_ilu_n_,mld_ilu_t_)
          ! do nothing
        case(mld_slu_)
        case(mld_umf_)
        case(mld_sludist_)
        case default
        end select
        
      end if
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_dp * size(prec%rprcparm)
!!$    if (allocated(prec%d))        val = val + psb_sizeof_dp * size(prec%d)
!!$    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
!!$    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
!!$                                  val = val + psb_sizeof(prec%desc_data)
!!$    if (allocated(prec%av))  then 
!!$      do i=1,size(prec%av)
!!$        val = val + psb_sizeof(prec%av(i))
!!$      end do
!!$    end if


  end function mld_dbaseprec_sizeof

  function mld_d_onelev_prec_sizeof(prec) result(val)
    implicit none 
    type(mld_donelev_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
!!$    val = val + mld_sizeof(prec%prec)
!!$    if (allocated(prec%iprcparm)) &
!!$         &  val = val + psb_sizeof_int * size(prec%iprcparm)
!!$    if (allocated(prec%ilaggr)) &
!!$         &  val = val + psb_sizeof_int * size(prec%ilaggr)
!!$    if (allocated(prec%nlaggr)) &
!!$         &  val = val + psb_sizeof_int * size(prec%nlaggr)
!!$    if (allocated(prec%rprcparm)) val = val + psb_sizeof_dp * size(prec%rprcparm)
    val = val + psb_sizeof(prec%desc_ac)
    val = val + psb_sizeof(prec%ac)
    val = val + psb_sizeof(prec%map) 
    if (allocated(prec%sm))  val = val + prec%sm%sizeof()
!!$    write(0,*) 'ONelev sizes: ',psb_sizeof(prec%desc_ac),&
!!$         & psb_sizeof(prec%ac),  psb_sizeof(prec%map), prec%sm%sizeof() 
  end function mld_d_onelev_prec_sizeof

  !
  ! Subroutine: mld_file_prec_descr
  ! Version: real
  !
  !  This routine prints a description of the preconditioner to the standard 
  !  output or to a file. It must be called after the preconditioner has been
  !  built by mld_precbld.
  !
  ! Arguments:
  !  p       -  type(mld_Tprec_type), input.
  !             The preconditioner data structure to be printed out.
  !  info    -  integer, output.
  !             error code.
  !  iout    -  integer, input, optional.
  !             The id of the file where the preconditioner description
  !             will be printed. If iout is not present, then the standard
  !             output is condidered.
  !
  subroutine mld_dfile_prec_descr(p,info,iout)
    implicit none 
    ! Arguments
    type(mld_dprec_type), intent(in) :: p
    integer, intent(out)             :: info
    integer, intent(in), optional    :: iout

    ! Local variables
    integer      :: ilev, nlev
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_file_prec_descr'
    integer :: iout_

    info = psb_success_
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if
    if (iout_ < 0) iout_ = 6 

    ictxt = p%ictxt

    if (allocated(p%precv)) then
!!$      ictxt = psb_cd_get_context(p%precv(1)%prec%desc_data)

      call psb_info(ictxt,me,np)

      !
      ! The preconditioner description is printed by processor psb_root_.
      ! This agrees with the fact that all the parameters defining the
      ! preconditioner have the same values on all the procs (this is
      ! ensured by mld_precbld).
      !
      if (me == psb_root_) then
        nlev = size(p%precv)
        do ilev = 1, nlev 
          if (.not.allocated(p%precv(ilev)%sm)) then 
            info = 3111
            write(iout_,*) ' ',name,&
                 & ': error: inconsistent MLPREC part, should call MLD_PRECINIT'
            return
          endif
        end do

        write(iout_,*) 
        write(iout_,'(a)') 'Preconditioner description'
        if (nlev >= 1) then
          !
          ! Print description of base preconditioner
          !
          if (nlev > 1) then
            write(iout_,*) 'Multilevel Schwarz'
            write(iout_,*) 
            write(iout_,*) 'Base preconditioner (smoother) details'
          endif
          call p%precv(1)%sm%descr(info,iout=iout_)
        end if

        if (nlev > 1) then

          !
          ! Print multilevel details
          !
          write(iout_,*) 
          write(iout_,*) 'Multilevel details'
          write(iout_,*) ' Number of levels: ',nlev

          !
          ! Currently, all the preconditioner parameters must have
          ! the same value at levels
          ! 2,...,nlev-1, hence only the values at level 2 are printed
          !

          ilev=2
          call p%precv(ilev)%parms%descr(iout_,info)

          !
          ! Coarse matrices are different at levels 2,...,nlev-1, hence related
          ! info is printed separately
          !
          write(iout_,*) 
          do ilev = 2, nlev-1
            write(iout_,*) ' Level ',ilev
            call p%precv(ilev)%descr(info,iout=iout_)
          end do

          !
          ! Print coarsest level details
          !

          ilev = nlev
          write(iout_,*) 
          write(iout_,*) ' Level ',ilev,' (coarsest)'

          call p%precv(ilev)%parms%descr(iout_,info,coarse=.true.)
          call p%precv(ilev)%descr(info,iout=iout_,coarse=.true.)
!!$          call mld_ml_new_coarse_descr(iout_,ilev,&
!!$               & p%precv(ilev)%iprcparm,&
!!$               & p%precv(ilev)%map%naggr,info,&
!!$               & dprcparm=p%precv(ilev)%rprcparm)
!!$          call p%precv(ilev)%sm%descr(info,iout=iout_)
        end if

      endif
      write(iout_,*) 
    else
      write(iout_,*) trim(name), &
           & ': Error: no base preconditioner available, something is wrong!'
      info = -2
      return
    endif


  end subroutine mld_dfile_prec_descr

  !
  ! Subroutines: mld_Tbase_precfree, mld_T_onelev_precfree, mld_Tprec_free
  ! Version: real/complex
  !
  !  These routines deallocate the mld_Tbaseprec_type, mld_Tonelev_type and
  !  mld_Tprec_type data structures.
  !
  ! Arguments:
  !  p       -  type(mld_Tbaseprec_type/mld_Tonelev_type/mld_Tprec_type), input.
  !             The data structure to be deallocated.
  !  info    -  integer, output.
  !             error code.
  !
  subroutine mld_dbase_precfree(p,info)
    implicit none 

    type(mld_dbaseprec_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = psb_success_

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff

!!$    if (allocated(p%d)) then 
!!$      deallocate(p%d,stat=info)
!!$    end if
!!$
!!$    if (allocated(p%av))  then 
!!$      do i=1,size(p%av) 
!!$        call p%av(i)%free()
!!$        if (info /= psb_success_) then 
!!$          ! Actually, we don't care here about this.
!!$          ! Just let it go.
!!$          ! return
!!$        end if
!!$      enddo
!!$      deallocate(p%av,stat=info)
!!$    end if
!!$
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if

!!$    if (allocated(p%perm)) then 
!!$      deallocate(p%perm,stat=info)
!!$    endif
!!$
!!$    if (allocated(p%invperm)) then 
!!$      deallocate(p%invperm,stat=info)
!!$    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(mld_prec_status_) == mld_prec_built_) then       
        if (p%iprcparm(mld_sub_solve_) == mld_slu_) then 
          call mld_dslu_free(p%iprcparm(mld_slu_ptr_),info)
        end if
        if (p%iprcparm(mld_sub_solve_) == mld_sludist_) then 
          call mld_dsludist_free(p%iprcparm(mld_slud_ptr_),info)
        end if
!!$        if (p%iprcparm(mld_sub_solve_) == mld_umf_) then 
!!$          call mld_dumf_free(p%iprcparm(mld_umf_symptr_),&
!!$               & p%iprcparm(mld_umf_numptr_),info)
!!$        end if
      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call mld_nullify_baseprec(p)

  end subroutine mld_dbase_precfree

  subroutine d_base_onelev_descr(lv,info,iout,coarse)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(in) :: lv
    integer, intent(out)                :: info
    integer, intent(in), optional       :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_base_onelev_descr'
    integer      :: iout_
    logical      :: coarse_


    call psb_erractionsave(err_act)

    if (present(coarse)) then 
      coarse_ = coarse
    else
      coarse_ = .false.
    end if
    if (present(iout)) then 
      iout_ = iout
    else 
      iout_ = 6
    end if

    if (lv%parms%ml_type > mld_no_ml_) then
      if (allocated(lv%map%naggr)) then
        write(iout_,*) '  Size of coarse matrix: ', &
             &  sum(lv%map%naggr(:))
        write(iout_,*) '  Sizes of aggregates: ', &
             &  lv%map%naggr(:)
      end if
      if (lv%parms%aggr_kind /= mld_no_smooth_) then
        write(iout_,*) '  Damping omega: ', &
             & lv%parms%aggr_omega_val
      end if
    end if
    if (allocated(lv%sm)) &
         & call lv%sm%descr(info,iout=iout_)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_onelev_descr

  subroutine mld_d_onelev_precfree(p,info)
    use psb_sparse_mod
    implicit none 

    type(mld_donelev_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = psb_success_

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff
    call p%sm%free(info)
!!$    call mld_precfree(p%prec,info)
    
    call p%ac%free()
    if (psb_is_ok_desc(p%desc_ac)) &
         & call psb_cdfree(p%desc_ac,info)
    
!!$    if (allocated(p%rprcparm)) then 
!!$      deallocate(p%rprcparm,stat=info)
!!$    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    !
    ! free explicitly map???
    ! For now thanks to allocatable semantics
    ! works anyway. 
    !

    call mld_nullify_onelevprec(p)
  end subroutine mld_d_onelev_precfree

  subroutine mld_nullify_dbaseprec(p)
    implicit none 

    type(mld_dbaseprec_type), intent(inout) :: p


  end subroutine mld_nullify_dbaseprec

  subroutine mld_nullify_d_onelevprec(p)
    implicit none 

    type(mld_donelev_type), intent(inout) :: p

    nullify(p%base_a) 
    nullify(p%base_desc) 

  end subroutine mld_nullify_d_onelevprec

  subroutine mld_dprec_free(p,info)
  
    use psb_sparse_mod
    
    implicit none
    
    ! Arguments
    type(mld_dprec_type), intent(inout) :: p
    integer, intent(out)                :: info
    
    ! Local variables
    integer             :: me,err_act,i
    character(len=20)   :: name
    
    if(psb_get_errstatus().ne.0) return 
    info=psb_success_
    name = 'mld_dprecfree'
    call psb_erractionsave(err_act)
    
    me=-1
    
    if (allocated(p%precv)) then 
      do i=1,size(p%precv) 
        call mld_precfree(p%precv(i),info)
      end do
      deallocate(p%precv)
    end if
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine mld_dprec_free


  subroutine d_base_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
    use psb_sparse_mod
    type(psb_desc_type), intent(in)             :: desc_data
    class(mld_d_base_smoother_type), intent(in) :: sm
    real(psb_dpk_),intent(in)                   :: x(:)
    real(psb_dpk_),intent(inout)                :: y(:)
    real(psb_dpk_),intent(in)                   :: alpha,beta
    character(len=1),intent(in)                 :: trans
    integer, intent(in)                         :: sweeps
    real(psb_dpk_),target, intent(inout)        :: work(:)
    integer, intent(out)                        :: info
    
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_apply'

    call psb_erractionsave(err_act)
    info = psb_success_
    if (allocated(sm%sv)) then 
      call sm%sv%apply(alpha,x,beta,y,desc_data,trans,work,info)
    else
      info = 1121
    endif
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
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
    
  end subroutine d_base_smoother_apply

  subroutine d_base_smoother_check(sm,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    if (allocated(sm%sv)) then 
      call sm%sv%check(info)
    else 
      info=3111
      call psb_errpush(info,name)
      goto 9999
    end if

    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_check


  subroutine d_base_smoother_seti(sm,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    integer, intent(in)                            :: what 
    integer, intent(in)                            :: val
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_seti'

    call psb_erractionsave(err_act)
    info = psb_success_

    if (allocated(sm%sv)) then 
      call sm%sv%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_seti

  subroutine d_base_smoother_setc(sm,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    integer, intent(in)                            :: what 
    character(len=*), intent(in)                   :: val
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_setc'

    call psb_erractionsave(err_act)

    info = psb_success_

    if (allocated(sm%sv)) then 
      call sm%sv%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_setc
  
  subroutine d_base_smoother_setr(sm,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    integer, intent(in)                            :: what 
    real(psb_dpk_), intent(in)                     :: val
    integer, intent(out)                           :: info
    Integer :: err_act
    character(len=20)  :: name='d_base_smoother_setr'

    call psb_erractionsave(err_act)


    info = psb_success_

    if (allocated(sm%sv)) then 
      call sm%sv%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_setr

  subroutine d_base_smoother_bld(a,desc_a,sm,upd,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    type(psb_dspmat_type), intent(in), target     :: a
    Type(psb_desc_type), Intent(in)                :: desc_a 
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    character, intent(in)                          :: upd
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_bld'

    call psb_erractionsave(err_act)

    info = psb_success_
    if (allocated(sm%sv)) then 
      call sm%sv%build(a,desc_a,upd,info)
    else
      info = 1121
      call psb_errpush(info,name)
    endif
    if (info /= psb_success_) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_bld


  subroutine d_base_smoother_free(sm,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_free'

    call psb_erractionsave(err_act)
    info = psb_success_
    
    if (allocated(sm%sv)) then 
      call sm%sv%free(info)
    end if
    if (info == psb_success_) deallocate(sm%sv,stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
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
  end subroutine d_base_smoother_free

  subroutine d_base_smoother_descr(sm,info,iout)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(in) :: sm
    integer, intent(out)                        :: info
    integer, intent(in), optional               :: iout

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_base_smoother_descr'
    integer :: iout_


    call psb_erractionsave(err_act)
    info = psb_success_

    if (present(iout)) then 
      iout_ = iout
    else 
      iout_ = 6
    end if

    write(iout_,*) 'Base smoother with local solver'
    if (allocated(sm%sv)) then 
      call sm%sv%descr(info,iout)
      if (info /= psb_success_) then 
        info = psb_err_from_subroutine_ 
        call psb_errpush(info,name,a_err='Local solver')
        goto 9999
      end if
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
  end subroutine d_base_smoother_descr

  function d_base_smoother_sizeof(sm) result(val)
    implicit none 
    ! Arguments
    class(mld_d_base_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_)                    :: val
    integer             :: i
    
    val = 0
    if (allocated(sm%sv)) then 
      val = sm%sv%sizeof()
    end if

    return
  end function d_base_smoother_sizeof

  subroutine d_base_smoother_default(sm) 
    implicit none 
    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm
    ! Do nothing for base version

    if (allocated(sm%sv)) call sm%sv%default()

    return
  end subroutine d_base_smoother_default



  subroutine d_base_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_sparse_mod
    type(psb_desc_type), intent(in)           :: desc_data
    class(mld_d_base_solver_type), intent(in) :: sv
    real(psb_dpk_),intent(in)                 :: x(:)
    real(psb_dpk_),intent(inout)              :: y(:)
    real(psb_dpk_),intent(in)                 :: alpha,beta
    character(len=1),intent(in)               :: trans
    real(psb_dpk_),target, intent(inout)      :: work(:)
    integer, intent(out)                      :: info
    
    Integer :: err_act
    character(len=20)  :: name='d_base_solver_apply'

    call psb_erractionsave(err_act)
    
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_base_solver_apply

  subroutine d_base_solver_bld(a,desc_a,sv,upd,info,b)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    type(psb_dspmat_type), intent(in), target   :: a
    Type(psb_desc_type), Intent(in)              :: desc_a 
    class(mld_d_base_solver_type), intent(inout) :: sv
    character, intent(in)                        :: upd
    integer, intent(out)                         :: info
    type(psb_dspmat_type), intent(in), target, optional  :: b
    Integer :: err_act
    character(len=20)  :: name='d_base_solver_bld'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_bld

  subroutine d_base_solver_check(sv,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_check'

    call psb_erractionsave(err_act)
    info = psb_success_


    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_check

  subroutine d_base_solver_seti(sv,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv 
    integer, intent(in)                          :: what 
    integer, intent(in)                          :: val
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_seti'
    
    ! Correct action here is doing nothing. 
    info = 0
    
    return
  end subroutine d_base_solver_seti

  subroutine d_base_solver_setc(sv,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv
    integer, intent(in)                          :: what 
    character(len=*), intent(in)                 :: val
    integer, intent(out)                         :: info
    Integer           :: err_act, ival 
    character(len=20) :: name='d_base_solver_setc'

    call psb_erractionsave(err_act)

    info = psb_success_

    call mld_stringval(val,ival,info)
    if (info == psb_success_) call sv%set(what,ival,info)

    if (info /= psb_success_) goto 9999


    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_setc
  
  subroutine d_base_solver_setr(sv,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv 
    integer, intent(in)                          :: what 
    real(psb_dpk_), intent(in)                   :: val
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_setr'

    
    ! Correct action here is doing nothing. 
    info = 0
    
    return
  end subroutine d_base_solver_setr

  subroutine d_base_solver_free(sv,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_free'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_free

  subroutine d_base_solver_descr(sv,info,iout)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(in) :: sv
    integer, intent(out)                      :: info
    integer, intent(in), optional             :: iout

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_base_solver_descr'
    integer      :: iout_


    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_descr

  function d_base_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_d_base_solver_type), intent(in) :: sv
    integer(psb_long_int_k_)                  :: val
    integer             :: i
    val = 0

    return
  end function d_base_solver_sizeof

  subroutine d_base_solver_default(sv) 
    implicit none 
    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv
    ! Do nothing for base version

    return
  end subroutine d_base_solver_default


  subroutine mld_d_apply2v(prec,x,y,desc_data,info,trans,work)
    use psb_sparse_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(mld_dprec_type), intent(in) :: prec
    real(psb_dpk_),intent(in)         :: x(:)
    real(psb_dpk_),intent(inout)      :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    real(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer           :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_dprec_type)
      call mld_precaply(prec,x,y,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
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

  end subroutine mld_d_apply2v

  subroutine mld_d_apply1v(prec,x,desc_data,info,trans)
    use psb_sparse_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(mld_dprec_type), intent(in) :: prec
    real(psb_dpk_),intent(inout)      :: x(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    Integer           :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_dprec_type)
      call mld_precaply(prec,x,desc_data,info,trans)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
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

  end subroutine mld_d_apply1v

  subroutine d_base_onelev_check(lv,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_onelev_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    call mld_check_def(lv%parms%sweeps,&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)
    call mld_check_def(lv%parms%sweeps_pre,&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)
    call mld_check_def(lv%parms%sweeps_post,&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)

    
    if (allocated(lv%sm)) then 
      call lv%sm%check(info)
    else 
      info=3111
      call psb_errpush(info,name)
      goto 9999
    end if

    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_onelev_check


  subroutine d_base_onelev_default(lv)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 

    lv%parms%sweeps          = 1
    lv%parms%sweeps_pre      = 1
    lv%parms%sweeps_post     = 1
    lv%parms%ml_type         = mld_mult_ml_
    lv%parms%aggr_alg        = mld_dec_aggr_
    lv%parms%aggr_kind       = mld_smooth_prol_
    lv%parms%coarse_mat      = mld_distr_mat_
    lv%parms%smoother_pos    = mld_twoside_smooth_
    lv%parms%aggr_omega_alg  = mld_eig_est_
    lv%parms%aggr_eig        = mld_max_norm_
    lv%parms%aggr_filter     = mld_no_filter_mat_
    lv%parms%aggr_omega_val  = dzero
    lv%parms%aggr_thresh     = dzero
    
    if (allocated(lv%sm)) call lv%sm%default()

    return

  end subroutine d_base_onelev_default


  subroutine d_base_onelev_seti(lv,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 
    integer, intent(in)                          :: what 
    integer, intent(in)                          :: val
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_onelev_seti'

    call psb_erractionsave(err_act)
    info = psb_success_

    select case (what) 

    case (mld_smoother_sweeps_)
      lv%parms%sweeps      = val
      lv%parms%sweeps_pre  = val
      lv%parms%sweeps_post = val

    case (mld_smoother_sweeps_pre_)
      lv%parms%sweeps_pre  = val

    case (mld_smoother_sweeps_post_)
      lv%parms%sweeps_post = val

    case (mld_ml_type_)
      lv%parms%ml_type       = val

    case (mld_aggr_alg_)
      lv%parms%aggr_alg      = val

    case (mld_aggr_kind_)
      lv%parms%aggr_kind     = val

    case (mld_coarse_mat_)
      lv%parms%coarse_mat    = val

    case (mld_smoother_pos_)
      lv%parms%smoother_pos  = val

    case (mld_aggr_omega_alg_)
      lv%parms%aggr_omega_alg= val

    case (mld_aggr_eig_)
      lv%parms%aggr_eig      = val

    case (mld_aggr_filter_)
      lv%parms%aggr_filter   = val

    case (mld_coarse_solve_)
      lv%parms%coarse_solve    = val

    case default
      if (allocated(lv%sm)) then 
        call lv%sm%set(what,val,info)
      end if
      if (info /= psb_success_) goto 9999
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
  end subroutine d_base_onelev_seti

  subroutine d_base_onelev_setc(lv,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 
    integer, intent(in)                            :: what 
    character(len=*), intent(in)                   :: val
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_onelev_setc'
    integer :: ival 

    call psb_erractionsave(err_act)

    info = psb_success_

    call mld_stringval(val,ival,info)
    if (info == psb_success_) call lv%set(what,ival,info)

    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_onelev_setc
  
  subroutine d_base_onelev_setr(lv,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 
    integer, intent(in)                            :: what 
    real(psb_dpk_), intent(in)                     :: val
    integer, intent(out)                           :: info
    Integer :: err_act
    character(len=20)  :: name='d_base_onelev_setr'

    call psb_erractionsave(err_act)


    info = psb_success_
    
    select case (what) 

    case (mld_aggr_omega_val_)
      lv%parms%aggr_omega_val= val
   
    case (mld_aggr_thresh_)
      lv%parms%aggr_thresh   = val

    case default
      if (allocated(lv%sm)) then 
        call lv%sm%set(what,val,info)
      end if
      if (info /= psb_success_) goto 9999
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
  end subroutine d_base_onelev_setr

  subroutine mld_d_dump(prec,info,istart,iend,prefix,head,ac,smoother,solver)
    use psb_sparse_mod
    implicit none 
    class(mld_dprec_type), intent(in) :: prec
    integer, intent(out)             :: info
    integer, intent(in), optional    :: istart, iend
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: smoother, solver,ac
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    !  len of prefix_ 

    info = 0

    iln = size(prec%precv)
    if (present(istart)) then 
      il1 = max(1,istart)
    else
      il1 = 2
    end if
    if (present(iend)) then 
      iln = min(iln, iend)
    end if

    do lev=il1, iln
      call prec%precv(lev)%dump(lev,info,prefix=prefix,head=head,&
           & ac=ac,smoother=smoother,solver=solver)
    end do

  end subroutine mld_d_dump
  

  subroutine d_base_onelev_dump(lv,level,info,prefix,head,ac,smoother,solver)
    use psb_sparse_mod
    implicit none 
    class(mld_donelev_type), intent(in) :: lv
    integer, intent(in)              :: level
    integer, intent(out)             :: info
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: ac, smoother, solver
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    logical :: ac_
    !  len of prefix_ 

    info = 0

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_lev_d"
    end if

    if (associated(lv%base_desc)) then 
      icontxt = psb_cd_get_context(lv%base_desc)
      call psb_info(icontxt,iam,np)
    else 
      icontxt = -1 
      iam     = -1
    end if
    if (present(ac)) then 
      ac_ = ac
    else
      ac_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    if (level >= 2) then 
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_ac.mtx'
      write(0,*) 'Filename ',fname
      if (ac_) call lv%ac%print(fname,head=head)
    end if
    if (allocated(lv%sm)) &
         & call lv%sm%dump(icontxt,level,info,smoother=smoother,solver=solver)

  end subroutine d_base_onelev_dump

  subroutine d_base_smoother_dmp(sm,ictxt,level,info,prefix,head,smoother,solver)
    use psb_sparse_mod
    implicit none 
    class(mld_d_base_smoother_type), intent(in) :: sm
    integer, intent(in)              :: ictxt,level
    integer, intent(out)             :: info
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: smoother, solver
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    logical :: smoother_
    !  len of prefix_ 

    info = 0

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_smth_d"
    end if

    call psb_info(ictxt,iam,np)

    if (present(smoother)) then 
      smoother_ = smoother
    else
      smoother_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    ! At base level do nothing for the smoother
    if (allocated(sm%sv)) &
         & call sm%sv%dump(ictxt,level,info,solver=solver)

  end subroutine d_base_smoother_dmp

  subroutine d_base_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
    use psb_sparse_mod
    implicit none 
    class(mld_d_base_solver_type), intent(in) :: sv
    integer, intent(in)              :: ictxt,level
    integer, intent(out)             :: info
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: solver
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    logical :: solver_
    !  len of prefix_ 

    info = 0

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_slv_d"
    end if

    call psb_info(ictxt,iam,np)

    if (present(solver)) then 
      solver_ = solver
    else
      solver_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    ! At base level do nothing for the solver

  end subroutine d_base_solver_dmp

!!$
!!$
!!$  subroutine mld_d_precdump_fact(prec,info,istart,iend,prefix,head)
!!$    use psb_base_mod
!!$    implicit none 
!!$    type(mld_dprec_type), intent(in) :: prec
!!$    integer, intent(out)             :: info
!!$    integer, intent(in), optional    :: istart, iend
!!$    character(len=*), intent(in), optional :: prefix,head
!!$    integer :: i, j, il1, iln, lname, lev
!!$    integer :: icontxt,iam, np
!!$    character(len=80)  :: prefix_
!!$    character(len=120) :: fname ! len should be at least 20 more than
!!$    !  len of prefix_ 
!!$
!!$    info = 0
!!$
!!$    if (.not.mld_is_asb(prec)) then 
!!$      info = -1
!!$      write(psb_err_unit,*) 'Trying to dump a non-built preconditioner'
!!$      return
!!$    end if
!!$
!!$    il1 = 1
!!$    iln = size(prec%precv)
!!$    if (present(istart)) then 
!!$      il1 = max(1,istart)
!!$    end if
!!$    if (present(iend)) then 
!!$      iln = min(iln, iend)
!!$    end if
!!$    if (present(prefix)) then 
!!$      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
!!$    else
!!$      prefix_ = "dump_fact_d"
!!$    end if
!!$
!!$    icontxt = psb_cd_get_context(prec%precv(1)%prec%desc_data)
!!$    call psb_info(icontxt,iam,np)
!!$    lname = len_trim(prefix_)
!!$    fname = trim(prefix_)
!!$    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
!!$    lname = lname + 5
!!$    do lev=il1, iln
!!$      write(fname(lname+1:),'(a,i3.3,a)')'_l',lev,'_lower.mtx'
!!$      if (psb_is_asb(prec%precv(lev)%prec%av(mld_l_pr_))) &
!!$           & call psb_csprt(fname,prec%precv(lev)%prec%av(mld_l_pr_),head=head)
!!$      write(fname(lname+1:),'(a,i3.3,a)')'_l',lev,'_diag.mtx'
!!$      if (allocated(prec%precv(lev)%prec%d)) &
!!$           & call psb_geprt(fname,prec%precv(lev)%prec%d,head=head)
!!$      write(fname(lname+1:),'(a,i3.3,a)')'_l',lev,'_upper.mtx'
!!$      if (psb_is_asb(prec%precv(lev)%prec%av(mld_u_pr_))) &
!!$           & call psb_csprt(fname,prec%precv(lev)%prec%av(mld_u_pr_),head=head)
!!$    end do
!!$
!!$  end subroutine mld_d_precdump_fact

end module mld_d_prec_type
