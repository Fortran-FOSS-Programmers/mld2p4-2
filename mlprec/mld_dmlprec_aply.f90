!!$
!!$ 
!!$                           MLD2P4  version 1.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.2)
!!$  
!!$  (C) Copyright 2008
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
!!$                      Alfredo Buttari      University of Rome Tor Vergata
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
! File: mld_dmlprec_aply.f90
!
! Subroutine: mld_dmlprec_aply
! Version:    real
!
!  This routine computes
!  
!                        Y = beta*Y + alpha*op(M^(-1))*X,
!  where 
!  - M is a multilevel domain decomposition (Schwarz) preconditioner associated
!    to a certain matrix A and stored in the array baseprecv,
!  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  For each level we have as many submatrices as processes (except for the coarsest
!  level where we might have a replicated index space) and each process takes care
!  of one submatrix.
!
!  The multilevel preconditioner M is regarded as an array of 'base preconditioners',
!  each representing the part of the preconditioner associated to a certain level.
!  For each level ilev, the base preconditioner K(ilev) is stored in baseprecv(ilev)
!  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
!  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
!  aggregation.
!
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level and A(1) is the matrix A.
!
!  For a general description of (parallel) multilevel preconditioners see
!    -  B.F. Smith, P.E. Bjorstad & W.D. Gropp,
!       Domain decomposition: parallel multilevel methods for elliptic partial
!       differential equations,
!       Cambridge University Press, 1996.
!    -  K. Stuben,
!       Algebraic Multigrid (AMG): An Introduction with Applications,
!       GMD Report N. 70, 1999.
!
!
! Arguments:
!   alpha      -   real(psb_dpk_), input.
!                  The scalar alpha.
!   baseprecv  -   type(mld_dbaseprc_type), dimension(:), input.
!                  The array of base preconditioner data structures containing the
!                  local parts of the preconditioners to be applied at each level.
!      Note that nlev = size(baseprecv) = number of levels.
!      baseprecv(ilev)%av  -  type(psb_dspmat_type), dimension(:), allocatable(:).
!                             The sparse matrices needed to apply the preconditioner 
!                             at level ilev. 
!         baseprecv(ilev)%av(mld_l_pr_)    -  The L factor of the ILU factorization of
!                                             the local diagonal block of A(ilev).
!         baseprecv(ilev)%av(mld_u_pr_)    -  The U factor of the ILU factorization of the
!                                             local diagonal block of A(ilev), except its
!                                             diagonal entries (stored in baseprecv(ilev)%d).
!         baseprecv(ilev)%av(mld_ap_nd_)   -  The entries of the local part of A(ilev)
!                                             outside the diagonal block, for block-Jacobi
!                                             sweeps.
!         baseprecv(ilev)%av(mld_ac_)      -  The local part of the matrix A(ilev).
!         baseprecv(ilev)%av(mld_sm_pr_)   -  The smoothed prolongator.   
!                                             It maps vectors (ilev) ---> (ilev-1).
!         baseprecv(ilev)%av(mld_sm_pr_t_) -  The smoothed prolongator transpose.   
!                                             It maps vectors (ilev-1) ---> (ilev).
!      baseprecv(ilev)%d         -  real(psb_dpk_), dimension(:), allocatable.
!                                   The diagonal entries of the U factor in the ILU
!                                   factorization of A(ilev).
!      baseprecv(ilev)%desc_data -  type(psb_desc_type).
!                                   The communication descriptor associated to the base
!                                   preconditioner, i.e. to the sparse matrices needed
!                                   to apply the base preconditioner at the current level.
!      baseprecv(ilev)%desc_ac   -  type(psb_desc_type).
!                                   The communication descriptor associated to the sparse
!                                   matrix A(ilev), stored in baseprecv(ilev)%av(mld_ac_).
!      baseprecv(ilev)%iprcparm  -  integer, dimension(:), allocatable.
!                                   The integer parameters defining the base
!                                   preconditioner K(ilev).
!      baseprecv(ilev)%rprcparm  -  real(psb_dpk_), dimension(:), allocatable.
!                                   The real parameters defining the base preconditioner
!                                   K(ilev).
!      baseprecv(ilev)%perm      -  integer, dimension(:), allocatable.
!                                   The row and column permutations applied to the local
!                                   part of A(ilev) (defined only if baseprecv(ilev)%
!                                   iprcparm(mld_sub_ren_)>0). 
!      baseprecv(ilev)%invperm   -  integer, dimension(:), allocatable.
!                                   The inverse of the permutation stored in
!                                   baseprecv(ilev)%perm.
!      baseprecv(ilev)%mlia      -  integer, dimension(:), allocatable.
!                                   The aggregation map (ilev-1) --> (ilev).
!                                   In case of non-smoothed aggregation, it is used
!                                   instead of mld_sm_pr_.
!      baseprecv(ilev)%nlaggr    -  integer, dimension(:), allocatable.
!                                   The number of aggregates (rows of A(ilev)) on the
!                                   various processes. 
!      baseprecv(ilev)%base_a    -  type(psb_dspmat_type), pointer.
!                                   Pointer (really a pointer!) to the base matrix of
!                                   the current level, i.e. the local part of A(ilev);
!                                   so we have a unified treatment of residuals. We
!                                   need this to avoid passing explicitly the matrix
!                                   A(ilev) to the routine which applies the
!                                   preconditioner.
!      baseprecv(ilev)%base_desc -  type(psb_desc_type), pointer.
!                                   Pointer to the communication descriptor associated
!                                   to the sparse matrix pointed by base_a.  
!                  
!   x          -  real(psb_dpk_), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  real(psb_dpk_), input.
!                 The scalar beta.
!   y          -  real(psb_dpk_), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned.
!   trans      -  character, optional.
!                 If trans='N','n' then op(M^(-1)) = M^(-1);
!                 if trans='T','t' then op(M^(-1)) = M^(-T) (transpose of M^(-1)).
!   work       -  real(psb_dpk_), dimension (:), optional, target.
!                 Workspace. Its size must be at least 4*psb_cd_get_local_cols(desc_data).
!   info       -  integer, output.
!                 Error code.
!
!   Note that when the LU factorization of the matrix A(ilev) is computed instead of
!   the ILU one, by using UMFPACK or SuperLU, the corresponding L and U factors
!   are stored in data structures provided by UMFPACK or SuperLU and pointed by
!   baseprecv(ilev)%iprcparm(mld_umf_ptr) or baseprecv(ilev)%iprcparm(mld_slu_ptr),
!   respectively.
!  
subroutine mld_dmlprec_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_inner_mod, mld_protect_name => mld_dmlprec_aply

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)      :: desc_data
  type(mld_dbaseprc_type), intent(in) :: baseprecv(:)
  real(psb_dpk_),intent(in)         :: alpha,beta
  real(psb_dpk_),intent(in)         :: x(:)
  real(psb_dpk_),intent(inout)      :: y(:)
  character, intent(in)               :: trans
  real(psb_dpk_),target             :: work(:)
  integer, intent(out)                :: info

  ! Local variables
  integer           :: ictxt, np, me, err_act
  integer           :: debug_level, debug_unit
  character(len=20) :: name
  character         :: trans_

  name='mld_dmlprec_aply'
  info = 0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Entry  ', size(baseprecv)

  trans_ = psb_toupper(trans)

  select case(baseprecv(2)%iprcparm(mld_ml_type_)) 

  case(mld_no_ml_)
    !
    ! No preconditioning, should not really get here
    ! 
    call psb_errpush(4001,name,a_err='mld_no_ml_ in mlprc_aply?')
    goto 9999      

  case(mld_add_ml_)
    !
    ! Additive multilevel
    !

    call add_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans_,work,info)

  case(mld_mult_ml_)
    ! 
    !  Multiplicative multilevel (multiplicative among the levels, additive inside
    !  each level)
    !
    !  Pre/post-smoothing versions.
    !  Note that the transpose switches pre <-> post.
    !

    select case(baseprecv(2)%iprcparm(mld_smooth_pos_))

    case(mld_post_smooth_)

      select case (trans_) 
      case('N')
        call mlt_post_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans_,work,info)
      case('T','C')
        call mlt_pre_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans_,work,info)
      case default
        info = 4001
        call psb_errpush(info,name,a_err='invalid trans')
        goto 9999      
      end select

    case(mld_pre_smooth_)

      select case (trans_) 
      case('N')
        call mlt_pre_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans_,work,info)
      case('T','C')
        call mlt_post_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans_,work,info)
      case default
        info = 4001
        call psb_errpush(info,name,a_err='invalid trans')
        goto 9999      
      end select

    case(mld_twoside_smooth_)

      call mlt_twoside_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans_,work,info)

    case default
      info = 4013
      call psb_errpush(info,name,a_err='invalid smooth_pos',&
           &  i_Err=(/baseprecv(2)%iprcparm(mld_smooth_pos_),0,0,0,0/))
      goto 9999      

    end select

  case default
    info = 4013
    call psb_errpush(info,name,a_err='invalid mltype',&
         &  i_Err=(/baseprecv(2)%iprcparm(mld_ml_type_),0,0,0,0/))
    goto 9999      

  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains
  !
  ! Subroutine: add_ml_aply
  ! Version:    real
  ! Note:       internal subroutine of mld_dmlprec_aply.
  !
  !  This routine computes
  !  
  !                        Y = beta*Y + alpha*op(M^(-1))*X,
  !  where 
  !  - M is an additive multilevel domain decomposition (Schwarz) preconditioner
  !    associated to a certain matrix A and stored in the array baseprecv,
  !  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
  !  - X and Y are vectors,
  !  - alpha and beta are scalars.
  !
  !  The preconditioner M is additive both through the levels and inside each
  !  level.
  !
  !  For each level we have as many submatrices as processes (except for the coarsest
  !  level where we might have a replicated index space) and each process takes care
  !  of one submatrix. 
  !
  !  The multilevel preconditioner M is regarded as an array of 'base preconditioners',
  !  each representing the part of the preconditioner associated to a certain level.
  !  For each level ilev, the base preconditioner K(ilev) is stored in baseprecv(ilev)
  !  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
  !  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
  !  aggregation.
  !
  !  The levels are numbered in increasing order starting from the finest one, i.e.
  !  level 1 is the finest level and A(1) is the matrix A. 
  !
  !  For details on the additive multilevel Schwarz preconditioner see the
  !  Algorithm 3.1.1 in the book:
  !    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
  !    Domain decomposition: parallel multilevel methods for elliptic partial
  !    differential equations, Cambridge University Press, 1996.
  !
  !  For a description of the arguments see mld_dmlprec_aply.
  !
  !  A sketch of the algorithm implemented in this routine is provided below
  !  (AV(ilev; sm_pr_) denotes the smoothed prolongator from level ilev to
  !  level ilev-1, while AV(ilev; sm_pr_t_) denotes its transpose, i.e. the
  !  corresponding restriction operator from level ilev-1 to level ilev).
  !
  !   1. ! Apply the base preconditioner at level 1.
  !      ! The sum over the subdomains is carried out in the
  !      ! application of K(1).
  !        X(1) = Xest
  !        Y(1) = (K(1)^(-1))*X(1)
  !
  !    2.  DO ilev=2,nlev
  !
  !         ! Transfer X(ilev-1) to the next coarser level.
  !           X(ilev) = AV(ilev; sm_pr_t_)*X(ilev-1)
  !
  !         ! Apply the base preconditioner at the current level.
  !         ! The sum over the subdomains is carried out in the
  !         ! application of K(ilev).
  !           Y(ilev) = (K(ilev)^(-1))*X(ilev)
  !
  !        ENDDO
  !
  !    3.  DO ilev=nlev-1,1,-1
  !
  !         ! Transfer Y(ilev+1) to the next finer level.
  !           Y(ilev) = AV(ilev+1; sm_pr_)*Y(ilev+1)
  !
  !        ENDDO
  !     
  !    4.  Yext = beta*Yext + alpha*Y(1)
  !
  subroutine add_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)

    implicit none 

    ! Arguments
    type(psb_desc_type),intent(in)      :: desc_data
    type(mld_dbaseprc_type), intent(in) :: baseprecv(:)
    real(psb_dpk_),intent(in)         :: alpha,beta
    real(psb_dpk_),intent(in)         :: x(:)
    real(psb_dpk_),intent(inout)      :: y(:)
    character, intent(in)               :: trans
    real(psb_dpk_),target             :: work(:)
    integer, intent(out)                :: info

    ! Local variables
    integer            :: n_row,n_col
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: ismth, nlev, ilev, icm
    character(len=20)  :: name

    type psb_mlprec_wrk_type
      real(psb_dpk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
    end type psb_mlprec_wrk_type
    type(psb_mlprec_wrk_type), allocatable  :: mlprec_wrk(:)

    name='add_ml_aply'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Entry  ', size(baseprecv)

    nlev = size(baseprecv)
    allocate(mlprec_wrk(nlev),stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    !
    ! STEP 1
    !
    ! Apply the base preconditioner at the finest level
    !
    allocate(mlprec_wrk(1)%x2l(size(x)),mlprec_wrk(1)%y2l(size(y)), stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/size(x)+size(y),0,0,0,0/),&
           & a_err='real(psb_dpk_)')
      goto 9999      
    end if

    mlprec_wrk(1)%x2l(:) = x(:) 
    mlprec_wrk(1)%y2l(:) = dzero 

    call mld_baseprec_aply(alpha,baseprecv(1),x,beta,y,&
         & baseprecv(1)%base_desc,trans,work,info)
    if (info /=0) then 
      call psb_errpush(4010,name,a_err='baseprec_aply')
      goto 9999
    end if
    !
    ! STEP 2
    !
    ! For each level except the finest one ...
    !
    do ilev = 2, nlev
      n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
      n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%base_desc)
      nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)
      allocate(mlprec_wrk(ilev)%x2l(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
           & stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/2*(nc2l+max(n_row,n_col)),0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if

      ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
      icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)
      
      ! Apply prolongator transpose, i.e. restriction
      call psb_forward_map(done,mlprec_wrk(ilev-1)%x2l,&
           & dzero,mlprec_wrk(ilev)%x2l,&
           & baseprecv(ilev)%map_desc,info,work=work)

      
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if


      !
      ! Apply the base preconditioner
      !
      call mld_baseprec_aply(done,baseprecv(ilev),&
           & mlprec_wrk(ilev)%x2l,dzero,mlprec_wrk(ilev)%y2l,&
           & baseprecv(ilev)%base_desc, trans,work,info)

    enddo

    !
    ! STEP 3
    !
    ! For each level except the finest one ...
    !
    do ilev =nlev,2,-1

      n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
      n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%base_desc)
      nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)
      ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
      icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)

      !
      ! Apply prolongator
      !  
      call psb_backward_map(done,mlprec_wrk(ilev)%y2l,&
           & done,mlprec_wrk(ilev-1)%y2l,&
           & baseprecv(ilev)%map_desc,info,work=work)

      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during prolongation')
        goto 9999
      end if
    end do

    !
    ! STEP 4
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,done,y,baseprecv(1)%base_desc,info)
    if (info /= 0) then
      call psb_errpush(4001,name,a_err='Error on final update')
      goto 9999
    end if

    deallocate(mlprec_wrk,stat=info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
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

  end subroutine add_ml_aply
  !
  ! Subroutine: mlt_pre_ml_aply
  ! Version:    real
  ! Note:       internal subroutine of mld_dmlprec_aply.
  !
  !  This routine computes
  !  
  !                        Y = beta*Y + alpha*op(M^(-1))*X,
  !  where 
  !  - M is a hybrid multilevel domain decomposition (Schwarz) preconditioner
  !    associated to a certain matrix A and stored in the array baseprecv,
  !  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
  !  - X and Y are vectors,
  !  - alpha and beta are scalars.
  !
  !  The preconditioner M is hybrid in the sense that it is multiplicative through the
  !  levels and additive inside a level; pre-smoothing only is applied at each level.
  !
  !  For each level we have as many submatrices as processes (except for the coarsest
  !  level where we might have a replicated index space) and each process takes care
  !  of one submatrix. 
  !
  !  The multilevel preconditioner M is regarded as an array of 'base preconditioners',
  !  each representing the part of the preconditioner associated to a certain level.
  !  For each level ilev, the base preconditioner K(ilev) is stored in baseprecv(ilev)
  !  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
  !  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
  !  aggregation.
  !
  !  The levels are numbered in increasing order starting from the finest one, i.e.
  !  level 1 is the finest level and A(1) is the matrix A. 
  !
  !  For details on the pre-smoothed hybrid multiplicative multilevel Schwarz
  !  preconditioner, see the Algorithm 3.2.1 in the book:
  !    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
  !    Domain decomposition: parallel multilevel methods for elliptic partial
  !    differential equations, Cambridge University Press, 1996.
  !
  !  For a description of the arguments see mld_dmlprec_aply.
  ! 
  !  A sketch of the algorithm implemented in this routine is provided below
  !  (AV(ilev; sm_pr_) denotes the smoothed prolongator from level ilev to
  !  level ilev-1, while AV(ilev; sm_pr_t_) denotes its transpose, i.e. the
  !  corresponding restriction operator from level ilev-1 to level ilev).
  !
  !    1.   X(1) = Xext
  !
  !    2. ! Apply the base preconditioner at the finest level.
  !         Y(1) = (K(1)^(-1))*X(1)
  !
  !    3. ! Compute the residual at the finest level.
  !         TX(1) = X(1) - A(1)*Y(1)
  !
  !    4.   DO ilev=2, nlev
  !
  !          ! Transfer the residual to the current (coarser) level.
  !            X(ilev) = AV(ilev; sm_pr_t_)*TX(ilev-1)
  !
  !          ! Apply the base preconditioner at the current level.
  !          ! The sum over the subdomains is carried out in the
  !          ! application of K(ilev).     
  !            Y(ilev) = (K(ilev)^(-1))*X(ilev)
  !
  !          ! Compute the residual at the current level (except at
  !          ! the coarsest level).
  !            IF (ilev < nlev)
  !               TX(ilev) = (X(ilev)-A(ilev)*Y(ilev))
  !
  !         ENDDO
  !
  !    5.   DO ilev=nlev-1,1,-1
  !
  !          ! Transfer Y(ilev+1) to the next finer level
  !            Y(ilev) = Y(ilev) + AV(ilev+1; sm_pr_)*Y(ilev+1)
  !
  !         ENDDO
  !     
  !    6.  Yext = beta*Yext + alpha*Y(1)
  ! 
  !
  subroutine mlt_pre_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)

    implicit none 

    ! Arguments
    type(psb_desc_type),intent(in)      :: desc_data
    type(mld_dbaseprc_type), intent(in) :: baseprecv(:)
    real(psb_dpk_),intent(in)         :: alpha,beta
    real(psb_dpk_),intent(in)         :: x(:)
    real(psb_dpk_),intent(inout)      :: y(:)
    character, intent(in)               :: trans
    real(psb_dpk_),target             :: work(:)
    integer, intent(out)                :: info

    ! Local variables
    integer            :: n_row,n_col
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: ismth, nlev, ilev, icm
    character(len=20)  :: name

    type psb_mlprec_wrk_type
      real(psb_dpk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
    end type psb_mlprec_wrk_type
    type(psb_mlprec_wrk_type), allocatable  :: mlprec_wrk(:)

    name='mlt_pre_ml_aply'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Entry  ', size(baseprecv)

    nlev = size(baseprecv)
    allocate(mlprec_wrk(nlev),stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    !
    ! STEP 1
    !
    ! Copy the input vector X
    !
    n_col = psb_cd_get_local_cols(desc_data)
    nc2l  = psb_cd_get_local_cols(baseprecv(1)%base_desc)

    allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
         & mlprec_wrk(1)%tx(nc2l), stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
           & a_err='real(psb_dpk_)')
      goto 9999
    end if

    mlprec_wrk(1)%x2l(:) = x
    !
    ! STEP 2
    !
    ! Apply the base preconditioner at the finest level
    !
    call mld_baseprec_aply(done,baseprecv(1),mlprec_wrk(1)%x2l,&
         &  dzero,mlprec_wrk(1)%y2l,baseprecv(1)%base_desc,&
         &  trans,work,info)
    if (info /=0) then
      call psb_errpush(4010,name,a_err=' baseprec_aply')
      goto 9999
    end if

    !
    ! STEP 3
    !
    ! Compute the residual at the finest level
    !
    mlprec_wrk(1)%tx = mlprec_wrk(1)%x2l

    call psb_spmm(-done,baseprecv(1)%base_a,mlprec_wrk(1)%y2l,&
         & done,mlprec_wrk(1)%tx,baseprecv(1)%base_desc,info,&
         & work=work,trans=trans)
    if (info /=0) then
      call psb_errpush(4001,name,a_err=' fine level residual')
      goto 9999
    end if

    !
    ! STEP 4
    !
    ! For each level but the finest one ...
    !
    do ilev = 2, nlev

      n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
      n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%base_desc)
      nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)
      ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
      icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)

      allocate(mlprec_wrk(ilev)%tx(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
           &   mlprec_wrk(ilev)%x2l(nc2l), stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if

      ! Apply prolongator transpose, i.e. restriction      
      call psb_forward_map(done,mlprec_wrk(ilev-1)%tx,&
           & dzero,mlprec_wrk(ilev)%x2l,&
           & baseprecv(ilev)%map_desc,info,work=work)
      
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      !
      ! Apply the base preconditioner
      !
      call mld_baseprec_aply(done,baseprecv(ilev),mlprec_wrk(ilev)%x2l,&
           & dzero,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%base_desc,trans,work,info)

      !
      ! Compute the residual (at all levels but the coarsest one)
      !
      if (ilev < nlev) then
        mlprec_wrk(ilev)%tx = mlprec_wrk(ilev)%x2l
        if (info == 0) call psb_spmm(-done,baseprecv(ilev)%base_a,&
             & mlprec_wrk(ilev)%y2l,done,mlprec_wrk(ilev)%tx,&
             & baseprecv(ilev)%base_desc,info,work=work,trans=trans)
      endif
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error on up sweep residual')
        goto 9999
      end if
    enddo

    !
    ! STEP 5
    !
    ! For each level but the coarsest one ...
    !
    do ilev = nlev-1, 1, -1

      ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
      n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

      !
      ! Apply prolongator
      !  
      call psb_backward_map(done,mlprec_wrk(ilev+1)%y2l,&
           & done,mlprec_wrk(ilev)%y2l,&
           & baseprecv(ilev+1)%map_desc,info,work=work)

      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during prolongation')
        goto 9999
      end if
    enddo

    !
    ! STEP 6
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,&
         &  baseprecv(1)%base_desc,info)
    if (info /=0) then
      call psb_errpush(4001,name,a_err='Error on final update')
      goto 9999
    end if

    deallocate(mlprec_wrk,stat=info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
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
  end subroutine mlt_pre_ml_aply
  !
  ! Subroutine: mlt_post_ml_aply
  ! Version:    real
  ! Note:       internal subroutine of mld_dmlprec_aply.
  !
  !  This routine computes
  !  
  !                        Y = beta*Y + alpha*op(M^(-1))*X,
  !  where 
  !  - M is a hybrid multilevel domain decomposition (Schwarz) preconditioner
  !    associated to a certain matrix A and stored in the array baseprecv,
  !  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
  !  - X and Y are vectors,
  !  - alpha and beta are scalars.
  !
  !  The preconditioner M is hybrid in the sense that it is multiplicative through the
  !  levels and additive inside a level; post-smoothing only is applied at each level.
  !
  !  For each level we have as many submatrices as processes (except for the coarsest
  !  level where we might have a replicated index space) and each process takes care
  !  of one submatrix. 
  !
  !  The multilevel preconditioner M is regarded as an array of 'base preconditioners',
  !  each representing the part of the preconditioner associated to a certain level.
  !  For each level ilev, the base preconditioner K(ilev) is stored in baseprecv(ilev)
  !  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
  !  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
  !  aggregation.
  !
  !  The levels are numbered in increasing order starting from the finest one, i.e.
  !  level 1 is the finest level and A(1) is the matrix A. 
  !
  !  For details on hybrid multiplicative multilevel Schwarz preconditioners, see
  !    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
  !    Domain decomposition: parallel multilevel methods for elliptic partial
  !    differential equations, Cambridge University Press, 1996.
  !
  !  For a description of the arguments see mld_dmlprec_aply.
  !
  !  A sketch of the algorithm implemented in this routine is provided below.
  !  (AV(ilev; sm_pr_) denotes the smoothed prolongator from level ilev to
  !  level ilev-1, while AV(ilev; sm_pr_t_) denotes its transpose, i.e. the
  !  corresponding restriction operator from level ilev-1 to level ilev).
  !
  !    1.  X(1) = Xext
  !
  !    2.  DO ilev=2, nlev
  !
  !         ! Transfer X(ilev-1) to the next coarser level.
  !           X(ilev) = AV(ilev; sm_pr_t_)*X(ilev-1) 
  !
  !        ENDDO
  !   
  !    3.! Apply the preconditioner at the coarsest level.
  !        Y(nlev) = (K(nlev)^(-1))*X(nlev)
  !
  !    4.  DO ilev=nlev-1,1,-1
  !
  !         ! Transfer Y(ilev+1) to the next finer level.
  !           Y(ilev) = AV(ilev+1; sm_pr_)*Y(ilev+1)
  !
  !         ! Compute the residual at the current level and apply to it the
  !         ! base preconditioner. The sum over the subdomains is carried out
  !         ! in the application of K(ilev).
  !           Y(ilev) = Y(ilev) + (K(ilev)^(-1))*(X(ilev)-A(ilev)*Y(ilev))
  !
  !        ENDDO
  !
  !    5.  Yext = beta*Yext + alpha*Y(1)
  ! 
  !
  subroutine mlt_post_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)

    implicit none 

    ! Arguments
    type(psb_desc_type),intent(in)      :: desc_data
    type(mld_dbaseprc_type), intent(in) :: baseprecv(:)
    real(psb_dpk_),intent(in)         :: alpha,beta
    real(psb_dpk_),intent(in)         :: x(:)
    real(psb_dpk_),intent(inout)      :: y(:)
    character, intent(in)               :: trans
    real(psb_dpk_),target             :: work(:)
    integer, intent(out)                :: info

    ! Local variables
    integer            :: n_row,n_col
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: ismth, nlev, ilev, icm
    character(len=20)  :: name

    type psb_mlprec_wrk_type
      real(psb_dpk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
    end type psb_mlprec_wrk_type
    type(psb_mlprec_wrk_type), allocatable  :: mlprec_wrk(:)

    name='mlt_post_ml_aply'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Entry  ', size(baseprecv)

    nlev = size(baseprecv)
    allocate(mlprec_wrk(nlev),stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    !
    ! STEP 1
    !
    ! Copy the input vector X
    !
    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' desc_data status',allocated(desc_data%matrix_data)    

    n_col = psb_cd_get_local_cols(desc_data)
    nc2l  = psb_cd_get_local_cols(baseprecv(1)%base_desc)

    allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
         & mlprec_wrk(1)%tx(nc2l), stat=info)

    call psb_geaxpby(done,x,dzero,mlprec_wrk(1)%tx,&
         & baseprecv(1)%base_desc,info)
    call psb_geaxpby(done,x,dzero,mlprec_wrk(1)%x2l,&
         & baseprecv(1)%base_desc,info)

    !
    ! STEP 2
    !
    ! For each level but the finest one ...
    !
    do ilev=2, nlev

      n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
      n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%base_desc)
      nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)
      ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
      icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)

      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name), &
           & ' starting up sweep ',&
           & ilev,allocated(baseprecv(ilev)%iprcparm),n_row,n_col,&
           & nc2l, nr2l,ismth

      allocate(mlprec_wrk(ilev)%tx(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
           &   mlprec_wrk(ilev)%x2l(nc2l), stat=info)

      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if

      ! Apply prolongator transpose, i.e. restriction
      call psb_forward_map(done,mlprec_wrk(ilev-1)%x2l,&
           & dzero,mlprec_wrk(ilev)%x2l,&
           & baseprecv(ilev)%map_desc,info,work=work)
      
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      !
      ! update x2l
      !
      call psb_geaxpby(done,mlprec_wrk(ilev)%x2l,dzero,mlprec_wrk(ilev)%tx,&
           & baseprecv(ilev)%base_desc,info)
      if (info /= 0) then
        call psb_errpush(4001,name,a_err='Error in update')
        goto 9999
      end if

      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' done up sweep ', ilev

    enddo

    !
    ! STEP 3
    !
    ! Apply the base preconditioner at the coarsest level
    !
    call mld_baseprec_aply(done,baseprecv(nlev),mlprec_wrk(nlev)%x2l, &
         & dzero, mlprec_wrk(nlev)%y2l,baseprecv(nlev)%base_desc,trans,work,info)

    if (info /=0) then
      call psb_errpush(4010,name,a_err='baseprec_aply')
      goto 9999
    end if

    if (debug_level >= psb_debug_inner_) write(debug_unit,*) &
         & me,' ',trim(name), ' done baseprec_aply ', nlev

    !
    ! STEP 4
    !
    ! For each level but the coarsest one ...
    !
    do ilev=nlev-1, 1, -1

      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' starting down sweep',ilev

      ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
      n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

      !
      ! Apply prolongator
      !  
      call psb_backward_map(done,mlprec_wrk(ilev+1)%y2l,&
           & dzero,mlprec_wrk(ilev)%y2l,&
           & baseprecv(ilev+1)%map_desc,info,work=work)

      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during prolongation')
        goto 9999
      end if

      !
      ! Compute the residual
      !
      call psb_spmm(-done,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
           & done,mlprec_wrk(ilev)%tx,baseprecv(ilev)%base_desc,info,&
           & work=work,trans=trans)

      !
      ! Apply the base preconditioner
      !
      if (info == 0) call mld_baseprec_aply(done,baseprecv(ilev),mlprec_wrk(ilev)%tx,&
           & done,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%base_desc,trans,work,info)
      if (info /=0) then
        call psb_errpush(4001,name,a_err=' spmm/baseprec_aply')
        goto 9999
      end if

      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' done down sweep',ilev
    enddo

    !
    ! STEP 5
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,baseprecv(1)%base_desc,info)

    if (info /=0) then
      call psb_errpush(4001,name,a_err=' Final update')
      goto 9999
    end if



    deallocate(mlprec_wrk,stat=info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
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
  end subroutine mlt_post_ml_aply
  !
  ! Subroutine: mlt_twoside_ml_aply
  ! Version:    real
  ! Note:       internal subroutine of mld_dmlprec_aply.
  !
  !  This routine computes
  !  
  !                        Y = beta*Y + alpha*op(M^(-1))*X,
  !  where 
  !  - M is a symmetrized hybrid multilevel domain decomposition (Schwarz)
  !    preconditioner associated to a certain matrix A and stored in the array
  !    baseprecv,
  !  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
  !  - X and Y are vectors,               
  !  - alpha and beta are scalars.
  !
  !  The preconditioner M is hybrid in the sense that it is multiplicative through
  !  the levels and additive inside a level; it is symmetrized since pre-smoothing
  !  and post-smoothing are applied at each level.
  !
  !  For each level we have as many submatrices as processes (except for the coarsest
  !  level where we might have a replicated index space) and each process takes care
  !  of one submatrix.   
  !
  !  The multilevel preconditioner M is regarded as an array of 'base preconditioners',
  !  each representing the part of the preconditioner associated to a certain level.
  !  For each level ilev, the base preconditioner K(ilev) is stored in baseprecv(ilev)
  !  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
  !  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
  !  aggregation.
  !
  !  The levels are numbered in increasing order starting from the finest one, i.e.
  !  level 1 is the finest level and A(1) is the matrix A. 
  !
  !  For details on the symmetrized hybrid multiplicative multilevel Schwarz
  !  preconditioner, see the Algorithm 3.2.2 of the book:
  !    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
  !    Domain decomposition: parallel multilevel methods for elliptic partial
  !    differential equations, Cambridge University Press, 1996.
  !
  !  For a description of the arguments see mld_dmlprec_aply.
  !
  !  A sketch of the algorithm implemented in this routine is provided below.
  !  (AV(ilev; sm_pr_) denotes the smoothed prolongator from level ilev to
  !  level ilev-1, while AV(ilev; sm_pr_t_) denotes its transpose, i.e. the
  !  corresponding restriction operator from level ilev-1 to level ilev).
  !
  !    1.   X(1)  = Xext
  !
  !    2. ! Apply the base peconditioner at the finest level
  !         Y(1)  = (K(1)^(-1))*X(1)
  !
  !    3. ! Compute the residual at the finest level
  !         TX(1) = X(1) - A(1)*Y(1)
  !
  !    4.   DO ilev=2, nlev
  !
  !          ! Transfer the residual to the current (coarser) level
  !            X(ilev) = AV(ilev; sm_pr_t)*TX(ilev-1)
  !    
  !          ! Apply the base preconditioner at the current level.
  !          ! The sum over the subdomains is carried out in the
  !          ! application of K(ilev)
  !            Y(ilev) = (K(ilev)^(-1))*X(ilev)
  !
  !          ! Compute the residual at the current level
  !            TX(ilev) = (X(ilev)-A(ilev)*Y(ilev))
  !
  !         ENDDO
  !
  !    5.   DO ilev=NLEV-1,1,-1
  !
  !          ! Transfer Y(ilev+1) to the next finer level
  !            Y(ilev) = Y(ilev) + AV(ilev+1; sm_pr_)*Y(ilev+1)
  !
  !          ! Compute the residual at the current level and apply to it the
  !          ! base preconditioner. The sum over the subdomains is carried out
  !          ! in the application of K(ilev)     
  !            Y(ilev) = Y(ilev) + (K(ilev)**(-1))*(X(ilev)-A(ilev)*Y(ilev))
  !
  !         ENDDO
  !
  !    6.  Yext = beta*Yext + alpha*Y(1)
  !
  subroutine mlt_twoside_ml_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)

    implicit none 

    ! Arguments
    type(psb_desc_type),intent(in)      :: desc_data
    type(mld_dbaseprc_type), intent(in) :: baseprecv(:)
    real(psb_dpk_),intent(in)         :: alpha,beta
    real(psb_dpk_),intent(in)         :: x(:)
    real(psb_dpk_),intent(inout)      :: y(:)
    character, intent(in)               :: trans
    real(psb_dpk_),target             :: work(:)
    integer, intent(out)                :: info

    ! Local variables
    integer            :: n_row,n_col
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: ismth, nlev, ilev, icm
    character(len=20)  :: name

    type psb_mlprec_wrk_type
      real(psb_dpk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
    end type psb_mlprec_wrk_type
    type(psb_mlprec_wrk_type), allocatable  :: mlprec_wrk(:)

    name='mlt_twoside_ml_aply'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Entry  ', size(baseprecv)

    nlev = size(baseprecv)
    allocate(mlprec_wrk(nlev),stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    ! STEP 1
    !
    ! Copy the input vector X
    !
    n_col = psb_cd_get_local_cols(desc_data)
    nc2l  = psb_cd_get_local_cols(baseprecv(1)%base_desc)

    allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
         & mlprec_wrk(1)%ty(nc2l), mlprec_wrk(1)%tx(nc2l), stat=info)

    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
           & a_err='real(psb_dpk_)')
      goto 9999
    end if

    call psb_geaxpby(done,x,dzero,mlprec_wrk(1)%x2l,&
         & baseprecv(1)%base_desc,info)
    call psb_geaxpby(done,x,dzero,mlprec_wrk(1)%tx,&
         & baseprecv(1)%base_desc,info)

    !
    ! STEP 2
    !
    ! Apply the base preconditioner at the finest level
    !
    call mld_baseprec_aply(done,baseprecv(1),mlprec_wrk(1)%x2l,&
         &  dzero,mlprec_wrk(1)%y2l,baseprecv(1)%base_desc,&
         &  trans,work,info)
    !
    ! STEP 3
    !
    ! Compute the residual at the finest level
    !
    mlprec_wrk(1)%ty = mlprec_wrk(1)%x2l
    if (info == 0) call psb_spmm(-done,baseprecv(1)%base_a,mlprec_wrk(1)%y2l,&
         & done,mlprec_wrk(1)%ty,baseprecv(1)%base_desc,info,&
         & work=work,trans=trans)
    if (info /=0) then
      call psb_errpush(4010,name,a_err='Fine level baseprec/residual')
      goto 9999
    end if

    !
    ! STEP 4
    !
    ! For each level but the finest one ...
    !
    do ilev = 2, nlev

      n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
      n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%base_desc)
      nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)
      ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
      icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)
      allocate(mlprec_wrk(ilev)%tx(nc2l),mlprec_wrk(ilev)%ty(nc2l),&
           &  mlprec_wrk(ilev)%y2l(nc2l),mlprec_wrk(ilev)%x2l(nc2l), stat=info)

      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if

      ! Apply prolongator transpose, i.e. restriction
      call psb_forward_map(done,mlprec_wrk(ilev-1)%ty,&
           & dzero,mlprec_wrk(ilev)%x2l,&
           & baseprecv(ilev)%map_desc,info,work=work)
      
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      call psb_geaxpby(done,mlprec_wrk(ilev)%x2l,dzero,mlprec_wrk(ilev)%tx,&
           & baseprecv(ilev)%base_desc,info)
      !
      ! Apply the base preconditioner
      !
      if (info == 0) call mld_baseprec_aply(done,baseprecv(ilev),&
           & mlprec_wrk(ilev)%x2l,dzero,mlprec_wrk(ilev)%y2l,&
           & baseprecv(ilev)%base_desc,trans,work,info)
      !
      ! Compute the residual (at all levels but the coarsest one)
      !
      if(ilev < nlev) then
        mlprec_wrk(ilev)%ty = mlprec_wrk(ilev)%x2l
        if (info == 0) call psb_spmm(-done,baseprecv(ilev)%base_a,&
             & mlprec_wrk(ilev)%y2l,done,mlprec_wrk(ilev)%ty,&
             & baseprecv(ilev)%base_desc,info,work=work,trans=trans)
      endif
      if (info /=0) then
        call psb_errpush(4001,name,a_err='baseprec_aply/residual')
        goto 9999
      end if

    enddo

    !
    ! STEP 5
    !
    ! For each level but the coarsest one ...
    !
    do ilev=nlev-1, 1, -1

      ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
      n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

      !
      ! Apply prolongator
      !  
      call psb_backward_map(done,mlprec_wrk(ilev+1)%y2l,&
           & done,mlprec_wrk(ilev)%y2l,&
           & baseprecv(ilev+1)%map_desc,info,work=work)

      if (info /=0 ) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      !
      ! Compute the residual
      !
      call psb_spmm(-done,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
           & done,mlprec_wrk(ilev)%tx,baseprecv(ilev)%base_desc,info,&
           & work=work,trans=trans)
      !
      ! Apply the base preconditioner
      !
      if (info == 0) call mld_baseprec_aply(done,baseprecv(ilev),mlprec_wrk(ilev)%tx,&
           & done,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%base_desc, trans, work,info)
      if (info /= 0) then
        call psb_errpush(4001,name,a_err='Error: residual/baseprec_aply')
        goto 9999
      end if
    enddo

    !
    ! STEP 6
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,&
         &   baseprecv(1)%base_desc,info)

    if (info /= 0) then
      call psb_errpush(4001,name,a_err='Error final update')
      goto 9999
    end if



    deallocate(mlprec_wrk,stat=info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
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
  end subroutine mlt_twoside_ml_aply

end subroutine mld_dmlprec_aply

