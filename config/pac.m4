dnl
dnl $Id$
dnl
dnl 20080206
dnl M4 macros for the PSBLAS library and useful for packages using PSBLAS.
dnl

dnl @synopsis PAC_CHECK_LIBS
dnl
dnl Tries to detect the presence of a specific function among various libraries, using AC_CHECK_LIB
dnl repeatedly on the specified libraries.
dnl 
dnl Example use:
dnl
dnl PAC_CHECK_LIBS([atlas blas],
dnl		[dgemm],
dnl		[have_dgemm=yes],
dnl		[have_dgemm=no])
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
dnl 20080211 modified slighty from original.
AC_DEFUN([PAC_CHECK_LIBS],
[
 pac_check_libs_ok=no
 [for pac_check_libs_f in $2 
 do ]
 [for pac_check_libs_l in $1 
 do ]
    if test x"$pac_check_libs_ok" == xno ; then
     AC_CHECK_LIB([$pac_check_libs_l],[$pac_check_libs_f], [pac_check_libs_ok=yes; pac_check_libs_LIBS="-l$pac_check_libs_l"],[],[$5])
    fi
  done
  done
 # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
 [ if test x"$pac_check_libs_ok" = xyes ; then
	$3
 else
        pac_check_libs_ok=no
        $4
 fi
 ]
])dnl 

dnl @synopsis PAC_FORTRAN_FUNC_MOVE_ALLOC( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program with move_alloc (a Fortran 2003 function).
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_FORTRAN_HAVE_MOVE_ALLOC],
ac_exeext=''
ac_ext='f'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([for MOVE_ALLOC intrinsic])
cat > conftest.$ac_ext <<EOF
           program test_move_alloc
               integer, allocatable :: a(:), b(:)
               allocate(a(3))
               call move_alloc(a, b)
               print *, allocated(a), allocated(b)
               print *, b
           end program test_move_alloc
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  AC_MSG_RESULT([yes])
  ifelse([$1], , :, [rm -rf conftest*
  $1])
else
  AC_MSG_RESULT([no])	
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$2], , , [  rm -rf conftest*
  $2
])dnl
fi
rm -f conftest*])



dnl @synopsis PAC_CHECK_HAVE_GFORTRAN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if MPIFC is $FC.
dnl The check will proceed by compiling a small Fortran program
dnl containing the __GNUC__ macro, which should be defined in the
dnl gfortran compiled programs.
dnl
dnl On pass, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_HAVE_GFORTRAN,
ac_exeext=''
ac_ext='F'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[
cat > conftest.$ac_ext <<EOF
           program main
#ifdef __GNUC__ 
              print *, "GCC!"
#else
        this program will fail
#endif
           end

EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  ifelse([$1], , :, [rm -rf conftest*
  $1])
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$2], , , [  rm -rf conftest*
  $2
])dnl
fi
rm -f conftest*])



dnl @synopsis PAC_HAVE_MODERN_GCC( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if the GNU fortran version is suitable for PSBLAS.
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl Note : Will use MPIFC; if unset, will use '$FC'.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_HAVE_MODERN_GCC,
ac_exeext=''
ac_ext='F'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[
cat > conftest.$ac_ext <<EOF
           program main
#if ( __GNUC__ >= 4 && __GNUC_MINOR__ > 2 ) || ( __GNUC__ > 4 )
              print *, "ok"
#else
        this program will fail
#endif
           end

EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  ifelse([$1], , :, [rm -rf conftest*
  $1])
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$2], , , [  rm -rf conftest*
  $2
])dnl
fi
rm -f conftest*])


dnl @synopsis PAC_FORTRAN_CHECK_HAVE_MPI_MOD( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will determine if the fortran compiler MPIFC needs to include mpi.h or needs
dnl to use the mpi module.
dnl
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl Modified Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_CHECK_HAVE_MPI_MOD,
ac_exeext=''
ac_ext='f90'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([MPI Fortran interface])
cat > conftest.$ac_ext <<EOF
           program test
             use mpi
           end program test
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  AC_MSG_RESULT([ use mpi ])
  ifelse([$1], , :, [rm -rf conftest*
  $1])
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
  AC_MSG_RESULT([ include mpif.h ])
ifelse([$2], , , [  rm -rf conftest*
  $2
])dnl
fi
rm -f conftest*])



dnl @synopsis PAC_ARG_WITH_FLAGS(lcase_name, UCASE_NAME)
dnl
dnl Test for --with-lcase_name="compiler/loader flags".  if defined, prepend 
dnl flags to standard UCASE_NAME definition.
dnl
dnl Use this macro to facilitate additional special flags that should be
dnl passed on to the preprocessor/compilers/loader.
dnl
dnl NOTE : Renamed after TAC_ARG_WITH_FLAGS as in the Trilinos-8.0.4 package.
dnl 
dnl NOTE : This macro works in a way the user should invoke
dnl         --with-flags=...
dnl	   only once, otherwise the first one will take effect.
dnl
dnl Example use:
dnl 
dnl PAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
dnl 
dnl tests for --with-cxxflags and pre-pends to CXXFLAGS
dnl 
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl @notes  Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_FLAGS],
[
AC_MSG_CHECKING([whether additional [$2] flags should be added (should be invoked only once)])
dnl AC_MSG_CHECKING([whether additional [$2] flags should be added])
AC_ARG_WITH($1,
AC_HELP_STRING([--with-$1], 
[additional [$2] flags to be added: will prepend to [$2]]),
[
$2="${withval} ${$2}"
AC_MSG_RESULT([$2 = ${$2}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis PAC_ARG_WITH_LIBS
dnl
dnl Test for --with-libs="name(s)".
dnl 
dnl Prepends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl note: Renamed after PAC_ARG_WITH_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WITH_LIBS
dnl 
dnl tests for --with-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([PAC_ARG_WITH_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(libs,
AC_HELP_STRING([--with-libs], 
[List additional link flags  here.  For example, --with-libs=-lspecial_system_lib
or --with-libs=-L/path/to/libs]),
[
LIBS="${withval} ${LIBS}"
AC_MSG_RESULT([LIBS = ${LIBS}])
],
AC_MSG_RESULT(no)
)
]
)


dnl @synopsis PAC_ARG_WITH_EXTRA_LIBS
dnl
dnl Test for --with-extra-libs="name(s)".
dnl 
dnl Appends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl note: Renamed after PAC_ARG_WITH_EXTRA_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WITH_EXTRA_LIBS
dnl 
dnl tests for --with-extra-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([PAC_ARG_WITH_EXTRA_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(extra-libs,
AC_HELP_STRING([--with-extra-libs], 
[List additional link flags  here.  For example, --with-extra-libs=-lspecial_system_lib
or --with-extra-libs=-L/path/to/libs]),
[
EXTRA_LIBS="${withval}"
AC_MSG_RESULT([EXTRA_LIBS = ${EXTRA_LIBS}])
],
AC_MSG_RESULT(no)
)
]
)

dnl @synopsis PAC_ARG_WITH_PSBLAS
dnl
dnl Test for --with-psblas="pathname".
dnl 
dnl Defines the path to PSBLAS build dir.
dnl
dnl note: Renamed after PAC_ARG_WITH_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WIT_PSBLAS
dnl 
dnl tests for --with-psblas and pre-pends to PSBLAS_PATH
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_PSBLAS],
[
AC_ARG_WITH(psblas,
AC_HELP_STRING([--with-psblas=DIR], [The install directory for PSBLAS, for example,
 --with-psblas=/opt/packages/psblas-3.3]),
[pac_cv_psblas_dir=$withval],
[pac_cv_psblas_dir=''])
AC_ARG_WITH(psblas-incdir, AC_HELP_STRING([--with-psblas-incdir=DIR], [Specify the directory for PSBLAS includes.]),
        [pac_cv_psblas_incdir=$withval],
        [pac_cv_psblas_incdir=''])
AC_ARG_WITH(psblas-libdir, AC_HELP_STRING([--with-psblas-libdir=DIR], [Specify the directory for PSBLAS library.]),
        [pac_cv_psblas_libdir=$withval],
        [pac_cv_psblas_libdir=''])
if test x"$pac_cv_psblas_incdir" == "x" ; then
   pac_cv_psblas_incdir="$pac_cv_psblas_dir/include";
fi
if test x"$pac_cv_psblas_libdir" == "x" ; then
   pac_cv_psblas_libdir="$pac_cv_psblas_dir/lib";
fi

])

dnl @synopsis PAC_FORTRAN_HAVE_PSBLAS( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program using the PSBLAS library
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_HAVE_PSBLAS,
ac_objext='.o'
ac_ext='f90'
ac_compile='${MPIFC-$FC} -c -o conftest${ac_objext} $FMFLAG$PSBLAS_DIR/include $FMFLAG$PSBLAS_DIR/lib conftest.$ac_ext  1>&5'
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([for working installation of PSBLAS])
cat > conftest.$ac_ext <<EOF
           program test
	       use psb_base_mod
           end program test
EOF
if AC_TRY_EVAL(ac_compile) && test -s conftest${ac_objext}; then
  ifelse([$1], , :, [rm -rf conftest*
  $1])
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$2], , , [  rm -rf conftest*
  $2
])dnl
fi
rm -f conftest*])


dnl @synopsis PAC_FORTRAN_PSBLAS_VERSION( )
dnl
dnl Will try to compile, link and run  a program using the PSBLAS library. \
dnl  Checks for version major,  minor and patchlevel
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_PSBLAS_VERSION,
ac_exeext=''
ac_objext='.o'
ac_ext='f90'
ac_compile='${MPIFC-$FC} -c -o conftest${ac_objext} $FMFLAG$PSBLAS_DIR/include $FMFLAG$PSBLAS_DIR/lib conftest.$ac_ext  1>&5'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $FMFLAG$PSBLAS_DIR/include -L$PSBLAS_DIR/lib -lpsb_base $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([for version of PSBLAS])
cat > conftest.$ac_ext <<EOF
           program test
	       use psb_base_mod
               print *,psb_version_major_
           end program test
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
 pac_cv_psblas_major=`./conftest${ac_exeext} | sed 's/^ *//'`
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
  pac_cv_psblas_major="unknown";
fi
cat > conftest.$ac_ext <<EOF
           program test
	       use psb_base_mod
               print *,psb_version_minor_
           end program test
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
 pac_cv_psblas_minor=`./conftest${ac_exeext} | sed 's/^ *//'`
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
  pac_cv_psblas_minor="unknown";
fi
cat > conftest.$ac_ext <<EOF
           program test
	       use psb_base_mod
               print *,psb_patchlevel_
           end program test
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
 pac_cv_psblas_patchlevel=`./conftest${ac_exeext} | sed 's/^ *//'`
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
  pac_cv_psblas_patchlevel="unknown";
fi
rm -f conftest*
AC_MSG_RESULT([Done])
]
)


dnl @synopsis PAC_FORTRAN_TEST_TR15581( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the TR15581 Fortran extension support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_TR15581,
ac_exeext=''
ac_ext='f90'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran allocatables TR15581])
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
cat > conftest.$ac_ext <<EOF
module conftest
  type outer
    integer,  allocatable :: v(:)
  end type outer

  interface foo
    module procedure foov, food
  end interface
contains

  subroutine foov(a,b)

    implicit none
    integer, allocatable, intent(inout) :: a(:)
    integer, allocatable, intent(out) :: b(:)


    allocate(b(size(a)))

  end subroutine foov
  subroutine food(a,b)

    implicit none
    type(outer), intent(inout) :: a
    type(outer), intent(out) :: b


    allocate(b%v(size(a%v)))

  end subroutine food

end module conftest



program testtr15581
  use conftest
  type(outer) :: da, db
  integer, allocatable :: a(:), b(:)

  allocate(a(10),da%v(10))
  a = (/ (i,i=1,10) /)
  da%v = (/ (i,i=1,10) /)
  call foo(a,b)
  call foo(da,db)
  write(*,*) b
  write(*,*) db%v

end program testtr15581
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  AC_MSG_RESULT([yes])
  ifelse([$1], , :, [
  $1])
else
  AC_MSG_RESULT([no])
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$2], , , [  
  $2
])dnl
fi
cd ..
rm -fr tmpdir_$i])
dnl @synopsis PAC_FORTRAN_TEST_VOLATILE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the VOLATILE Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_VOLATILE,
ac_exeext=''
ac_ext='f90'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran VOLATILE])
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
cat > conftest.$ac_ext <<EOF
program conftest
  integer, volatile :: i, j
end program conftest
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  AC_MSG_RESULT([yes])
  ifelse([$1], , :, [
  $1])
else
  AC_MSG_RESULT([no])
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$2], , , [  
  $2
])dnl
fi
cd ..
rm -fr tmpdir_$i])

dnl @synopsis PAC_CHECK_BLACS
dnl
dnl Will try to find the BLACS
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_BLACS,
[AC_ARG_WITH(blacs, AC_HELP_STRING([--with-blacs=LIB], [Specify BLACSLIBNAME or -lBLACSLIBNAME or the absolute library filename.]),
        [psblas_cv_blacs=$withval],
        [psblas_cv_blacs=''])

case $psblas_cv_blacs in
	yes | "") ;;
	-* | */* | *.a | *.so | *.so.* | *.o) 
	     BLACS_LIBS="$psblas_cv_blacs" ;;
	*) BLACS_LIBS="-l$psblas_cv_blacs" ;;
esac

#
# Test user-defined BLACS
#
if test x"$psblas_cv_blacs" != "x" ; then
      save_LIBS="$LIBS";
      AC_LANG([Fortran])
      LIBS="$BLACS_LIBS $LIBS"
      AC_MSG_CHECKING([for dgesd2d in $BLACS_LIBS])
      AC_TRY_LINK_FUNC(dgesd2d, [psblas_cv_blacs_ok=yes], [psblas_cv_blacs_ok=no;BLACS_LIBS=""])
      AC_MSG_RESULT($psblas_cv_blacs_ok)

     if test x"$psblas_cv_blacs_ok" == x"yes";  then 
     AC_MSG_CHECKING([for blacs_pinfo in $BLACS_LIBS])
     AC_TRY_LINK_FUNC(blacs_pinfo, [psblas_cv_blacs_ok=yes], [psblas_cv_blacs_ok=no;BLACS_LIBS=""])
     AC_MSG_RESULT($psblas_cv_blacs_ok)
     fi 
     LIBS="$save_LIBS";
fi
AC_LANG([C])	

######################################
# System BLACS with PESSL default names. 
######################################
if test x"$BLACS_LIBS" == "x" ; then
   AC_LANG([Fortran])
   PAC_CHECK_LIBS([blacssmp blacsp2 blacs], 
	[dgesd2d],
	[psblas_cv_blacs_ok=yes; LIBS="$LIBS $pac_check_libs_LIBS "  ]
	[BLACS_LIBS="$pac_check_libs_LIBS" ]
	AC_MSG_NOTICE([BLACS libraries detected.]),[]
    )
    if test x"$BLACS_LIBS" != "x"; then 
          save_LIBS="$LIBS";
          LIBS="$BLACS_LIBS $LIBS"
          AC_MSG_CHECKING([for blacs_pinfo in $BLACS_LIBS])
          AC_LANG([Fortran])
	  AC_TRY_LINK_FUNC(blacs_pinfo, [psblas_cv_blacs_ok=yes], [psblas_cv_blacs_ok=no;BLACS_LIBS=""])
          AC_MSG_RESULT($psblas_cv_blacs_ok)
          LIBS="$save_LIBS";	
    fi 
fi
######################################
# Maybe we're looking at PESSL BLACS?#
######################################
if  test x"$BLACS_LIBS" != "x" ; then
    save_LIBS="$LIBS";
    LIBS="$BLACS_LIBS $LIBS"
    AC_MSG_CHECKING([for PESSL BLACS])
    AC_LANG([Fortran])
    AC_TRY_LINK_FUNC(esvemonp, [psblas_cv_pessl_blacs=yes], [psblas_cv_pessl_blacs=no])
    AC_MSG_RESULT($psblas_cv_pessl_blacs)
    LIBS="$save_LIBS";
fi    
if test "x$psblas_cv_pessl_blacs" == "xyes";  then
   FDEFINES="$psblas_cv_define_prepend-DHAVE_ESSL_BLACS $FDEFINES"
fi 
    

##############################################################################
#	Netlib BLACS library with default names
##############################################################################

if test x"$BLACS_LIBS" == "x" ; then
   save_LIBS="$LIBS";
   AC_LANG([Fortran])
   PAC_CHECK_LIBS([ blacs_MPI-LINUX-0 blacs_MPI-SP5-0 blacs_MPI-SP4-0 blacs_MPI-SP3-0 blacs_MPI-SP2-0 blacsCinit_MPI-ALPHA-0 blacsCinit_MPI-IRIX64-0 blacsCinit_MPI-RS6K-0 blacsCinit_MPI-SPP-0 blacsCinit_MPI-SUN4-0 blacsCinit_MPI-SUN4SOL2-0 blacsCinit_MPI-T3D-0 blacsCinit_MPI-T3E-0 
	], 
	[dgesd2d],
	[psblas_cv_blacs_ok=yes; LIBS="$LIBS $pac_check_libs_LIBS " 
	psblas_have_netlib_blacs=yes;  ]
	[BLACS_LIBS="$pac_check_libs_LIBS" ]
	AC_MSG_NOTICE([BLACS libraries detected.]),[]
    )
    
    if test x"$BLACS_LIBS" != "x" ; then	
      AC_LANG([Fortran])	   
      PAC_CHECK_LIBS([ blacsF77init_MPI-LINUX-0 blacsF77init_MPI-SP5-0 blacsF77init_MPI-SP4-0 blacsF77init_MPI-SP3-0 blacsF77init_MPI-SP2-0 blacsF77init_MPI-ALPHA-0 blacsF77init_MPI-IRIX64-0 blacsF77init_MPI-RS6K-0 blacsF77init_MPI-SPP-0 blacsF77init_MPI-SUN4-0 blacsF77init_MPI-SUN4SOL2-0 blacsF77init_MPI-T3D-0 blacsF77init_MPI-T3E-0 
 	], 
	[blacs_pinfo],
	[psblas_cv_blacs_ok=yes; LIBS="$pac_check_libs_LIBS $LIBS" ]
	[BLACS_LIBS="$pac_check_libs_LIBS $BLACS_LIBS" ]
	AC_MSG_NOTICE([Netlib BLACS Fortran initialization libraries detected.]),[]
       )
    fi

    if test x"$BLACS_LIBS" != "x" ; then	
    
      AC_LANG([C])
      PAC_CHECK_LIBS([ blacsCinit_MPI-LINUX-0 blacsCinit_MPI-SP5-0 blacsCinit_MPI-SP4-0 blacsCinit_MPI-SP3-0 blacsCinit_MPI-SP2-0 blacsCinit_MPI-ALPHA-0 blacsCinit_MPI-IRIX64-0 blacsCinit_MPI-RS6K-0 blacsCinit_MPI-SPP-0 blacsCinit_MPI-SUN4-0 blacsCinit_MPI-SUN4SOL2-0 blacsCinit_MPI-T3D-0 blacsCinit_MPI-T3E-0 
	], 
	[Cblacs_pinfo],
	[psblas_cv_blacs_ok=yes; LIBS="$pac_check_libs_LIBS $LIBS" ]
	[BLACS_LIBS="$BLACS_LIBS $pac_check_libs_LIBS" ]
	AC_MSG_NOTICE([Netlib BLACS C initialization libraries detected.]),[]
       )
    fi
    LIBS="$save_LIBS";	
fi

if test x"$BLACS_LIBS" == "x" ; then
	AC_MSG_ERROR([
	No BLACS library detected! $PACKAGE_NAME will be unusable.
	Please make sure a BLACS implementation is accessible (ex.: --with-blacs="-lblacsname -L/blacs/dir" )
	])
else 
      save_LIBS="$LIBS";
      LIBS="$BLACS_LIBS $LIBS"
      AC_MSG_CHECKING([for ksendid in $BLACS_LIBS])
      AC_LANG([Fortran])
      AC_TRY_LINK_FUNC(ksendid, [psblas_cv_have_sendid=yes],[psblas_cv_have_sendid=no])
      AC_MSG_RESULT($psblas_cv_have_sendid)
      LIBS="$save_LIBS"
      AC_LANG([C])
      if test "x$psblas_cv_have_sendid" == "xyes";  then
        FDEFINES="$psblas_cv_define_prepend-DHAVE_KSENDID $FDEFINES"
      fi 
fi
])dnl


dnl @synopsis PAC_MAKE_IS_GNUMAKE
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
define(PAC_MAKE_IS_GNUMAKE,[
AC_MSG_CHECKING(for gnumake)
MAKE=${MAKE:-make}

if $MAKE --version 2>&1 | grep -e"GNU Make" >/dev/null; then 
    AC_MSG_RESULT(yes)
    psblas_make_gnumake='yes'
else
    AC_MSG_RESULT(no)
    psblas_make_gnumake='no'
fi
])dnl

dnl @synopsis PAC_CHECK_UMFPACK
dnl
dnl Will try to find the UMFPACK library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_UMFPACK,
[AC_ARG_WITH(umfpack, AC_HELP_STRING([--with-umfpack=LIBNAME], [Specify the library name for UMFPACK library. 
Default: "-lumfpack -lamd"]),
        [mld2p4_cv_umfpack=$withval],
        [mld2p4_cv_umfpack='-lumfpack -lamd'])
AC_ARG_WITH(umfpackdir, AC_HELP_STRING([--with-umfpackdir=DIR], [Specify the directory for UMFPACK library and includes.]),
        [mld2p4_cv_umfpackdir=$withval],
        [mld2p4_cv_umfpackdir=''])
AC_ARG_WITH(umfpackincdir, AC_HELP_STRING([--with-umfpackincdir=DIR], [Specify the directory for UMFPACK includes.]),
        [mld2p4_cv_umfpackincdir=$withval],
        [mld2p4_cv_umfpackincdir=''])
AC_ARG_WITH(umfpacklibdir, AC_HELP_STRING([--with-umfpacklibdir=DIR], [Specify the directory for UMFPACK library.]),
        [mld2p4_cv_umfpacklibdir=$withval],
        [mld2p4_cv_umfpacklibdir=''])

AC_LANG([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
if test "x$mld2p4_cv_umfpackincdir" != "x"; then 
 AC_MSG_NOTICE([umfp include dir $mld2p4_cv_umfpackincdir])
 UMF_INCLUDES="-I$mld2p4_cv_umfpackincdir"
elif test "x$mld2p4_cv_umfpackdir" != "x"; then 
 AC_MSG_NOTICE([umfp dir $mld2p4_cv_umfpackdir])
 UMF_INCLUDES="-I$mld2p4_cv_umfpackdir" 
fi
CPPFLAGS="$UMF_INCLUDES $CPPFLAGS"
AC_CHECK_HEADER([umfpack.h],
 [pac_umf_header_ok=yes],
 [pac_umf_header_ok=no; UMF_INCLUDES=""])
if test "x$mld2p4_cv_umfpacklibdir" != "x"; then 
   LIBS="-L$mld2p4_cv_umfpacklibdir $LIBS"
   UMF_LIBDIR="-L$mld2p4_cv_umfpacklibdir"
elif test "x$mld2p4_cv_umfpackdir" != "x"; then 
   LIBS="-L$mld2p4_cv_umfpackdir $LIBS"
   UMF_LIBDIR="-L$mld2p4_cv_umfpackdir"
fi
if test "x$pac_umf_header_ok" == "xno" ; then
  unset ac_cv_header_umfpack_h
  UMF_INCLUDES="-I$mld2p4_cv_umfpackdir" 
  CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"
  AC_CHECK_HEADER([umfpack.h],
  [pac_umf_header_ok=yes],
  [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi
if test "x$pac_umf_header_ok" == "xno" ; then
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_umfpack_h
  UMF_INCLUDES="-I$mld2p4_cv_umfpackdir/include -I$mld2p4_cv_umfpackdir/Include "
  CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"

  AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_INCLUDES])
  AC_CHECK_HEADER([umfpack.h],
    [pac_umf_header_ok=yes],
    [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi
if test "x$pac_umf_header_ok" == "xno" ; then
dnl Maybe new structure with UMFPACK UFconfig AMD? 
   unset ac_cv_header_umfpack_h
   UMF_INCLUDES="-I$mld2p4_cv_umfpackdir/UFconfig -I$mld2p4_cv_umfpackdir/UMFPACK/Include -I$mld2p4_cv_umfpackdir/AMD/Include"
   CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"
   AC_CHECK_HEADER([umfpack.h],
     [pac_umf_header_ok=yes],
     [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi


if test "x$pac_umf_header_ok" == "xyes" ; then 
      UMF_LIBS="$mld2p4_cv_umfpack $UMF_LIBDIR"
      LIBS="$UMF_LIBS -lm $LIBS";
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [mld2p4_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [mld2p4_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     if test "x$pac_umf_lib_ok" == "xno" ; then 
        dnl Maybe Lib or lib? 
        UMF_LIBDIR="-L$mld2p4_cv_umfpackdir/Lib -L$mld2p4_cv_umfpackdir/lib"
        UMF_LIBS="$mld2p4_cv_umfpack $UMF_LIBDIR -lm $SAVE_LIBS"
        LIBS="$UMF_LIBS"
        
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [mld2p4_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [mld2p4_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     fi
     if test "x$pac_umf_lib_ok" == "xno" ; then 
        dnl Maybe UMFPACK/Lib? 
        UMF_LIBDIR="-L$mld2p4_cv_umfpackdir/AMD/Lib -L$mld2p4_cv_umfpackdir/UMFPACK/Lib"
        UMF_LIBS="$mld2p4_cv_umfpack $UMF_LIBDIR -lm $SAVE_LIBS"
             LIBS="$UMF_LIBS"
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [mld2p4_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [mld2p4_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     fi
fi
LIBS="$SAVE_LIBS";
CPPFLAGS="$SAVE_CPPFLAGS";
])dnl 

dnl @synopsis PAC_CHECK_SUPERLU
dnl
dnl Will try to find the SUPERLU library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_SUPERLU,
[AC_ARG_WITH(superlu, AC_HELP_STRING([--with-superlu=LIBNAME], [Specify the library name for SUPERLU library.
Default: "-lsuperlu"]),
        [mld2p4_cv_superlu=$withval],
        [mld2p4_cv_superlu='-lsuperlu'])
AC_ARG_WITH(superludir, AC_HELP_STRING([--with-superludir=DIR], [Specify the directory for SUPERLU library and includes.]),
        [mld2p4_cv_superludir=$withval],
        [mld2p4_cv_superludir=''])
AC_ARG_WITH(superluincdir, AC_HELP_STRING([--with-superluincdir=DIR], [Specify the directory for SUPERLU includes.]),
        [mld2p4_cv_superluincdir=$withval],
        [mld2p4_cv_superluincdir=''])
AC_ARG_WITH(superlulibdir, AC_HELP_STRING([--with-superlulibdir=DIR], [Specify the directory for SUPERLU library.]),
        [mld2p4_cv_superlulibdir=$withval],
        [mld2p4_cv_superlulibdir=''])
AC_LANG([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
if test "x$mld2p4_cv_superluincdir" != "x"; then 
 AC_MSG_NOTICE([slu include dir $mld2p4_cv_superluincdir])
 SLU_INCLUDES="-I$mld2p4_cv_superluincdir"
elif test "x$mld2p4_cv_superludir" != "x"; then 
 AC_MSG_NOTICE([slu dir $mld2p4_cv_superludir])
 SLU_INCLUDES="-I$mld2p4_cv_superludir"
fi
if test "x$mld2p4_cv_superlulibdir" != "x"; then 
   SLU_LIBS="-L$mld2p4_cv_superlulibdir"
elif test "x$mld2p4_cv_superludir" != "x"; then 
   SLU_LIBS="-L$mld2p4_cv_superludir"
fi

LIBS="$SLU_LIBS $LIBS"
CPPFLAGS="$SLU_INCLUDES $CPPFLAGS"
AC_CHECK_HEADER([slu_ddefs.h],
		[pac_slu_header_ok=yes],
		[pac_slu_header_ok=no; SLU_INCLUDES=""])
if test "x$pac_slu_header_ok" == "xno" ; then 
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_slu_ddefs_h
  SLU_INCLUDES="-I$mld2p4_cv_superludir/include -I$mld2p4_cv_superludir/Include "
  CPPFLAGS="$SLU_INCLUDES $SAVE_CPPFLAGS"

 AC_CHECK_HEADER([slu_ddefs.h],
		 [pac_slu_header_ok=yes],
		 [pac_slu_header_ok=no; SLU_INCLUDES=""])
fi

if test "x$pac_slu_header_ok" == "xyes" ; then 
 SLU_LIBS="$mld2p4_cv_superlu $SLU_LIBS"
 LIBS="$SLU_LIBS -lm $LIBS";
 AC_MSG_CHECKING([for superlu_malloc in $SLU_LIBS])
 AC_TRY_LINK_FUNC(superlu_malloc, 
		  [mld2p4_cv_have_superlu=yes;pac_slu_lib_ok=yes;],
		  [mld2p4_cv_have_superlu=no;pac_slu_lib_ok=no; SLU_LIBS=""; ])
 if test "x$pac_slu_lib_ok" == "xno" ; then 
    dnl Maybe lib?
    SLU_LIBS="$mld2p4_cv_superlu -L$mld2p4_cv_superludir/lib";
    LIBS="$SLU_LIBS -lm $SAVE_LIBS";
    AC_TRY_LINK_FUNC(superlu_malloc, 
		     [mld2p4_cv_have_superlu=yes;pac_slu_lib_ok=yes;],
		     [mld2p4_cv_have_superlu=no;pac_slu_lib_ok=no; SLU_LIBS=""; SLU_INCLUDES=""])
 fi
 AC_MSG_RESULT($pac_slu_lib_ok)
fi
LIBS="$SAVE_LIBS";
CPPFLAGS="$SAVE_CPPFLAGS";
])dnl 

dnl @synopsis PAC_CHECK_SUPERLU_Dist
dnl
dnl Will try to find the SUPERLU_Dist library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_SUPERLUDIST,
[AC_ARG_WITH(superludist, AC_HELP_STRING([--with-superludist=LIBNAME], [Specify the libname for SUPERLUDIST library. Requires you also specify SuperLU. Default: "-lsuperlu_dist"]),
        [mld2p4_cv_superludist=$withval],
        [mld2p4_cv_superludist='-lsuperlu_dist'])
AC_ARG_WITH(superludistdir, AC_HELP_STRING([--with-superludistdir=DIR], [Specify the directory for SUPERLUDIST library and includes.]),
        [mld2p4_cv_superludistdir=$withval],
        [mld2p4_cv_superludistdir=''])

AC_ARG_WITH(superludistincdir, AC_HELP_STRING([--with-superludistincdir=DIR], [Specify the directory for SUPERLUDIST includes.]),
        [mld2p4_cv_superludistincdir=$withval],
        [mld2p4_cv_superludistincdir=''])

AC_ARG_WITH(superludistlibdir, AC_HELP_STRING([--with-superludistlibdir=DIR], [Specify the directory for SUPERLUDIST library.]),
        [mld2p4_cv_superludistlibdir=$withval],
        [mld2p4_cv_superludistlibdir=''])

AC_LANG([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
save_CC="$CC"
CC=${MPICC}
if test "x$mld2p4_cv_superludistincdir" != "x"; then 
 AC_MSG_NOTICE([sludist dir $mld2p4_cv_superludistincdir]) 
 SLUDIST_INCLUDES="-I$mld2p4_cv_superludistincdir"
elif test "x$mld2p4_cv_superludistdir" != "x"; then 
 AC_MSG_NOTICE([sludist dir $mld2p4_cv_superludistdir]) 
 SLUDIST_INCLUDES="-I$mld2p4_cv_superludistdir"
fi
if test "x$mld2p4_cv_superludistlibdir" != "x"; then 
   SLUDIST_LIBS="-L$mld2p4_cv_superludistlibdir"
elif test "x$mld2p4_cv_superludistdir" != "x"; then 
   SLUDIST_LIBS="-L$mld2p4_cv_superludir"
fi

LIBS="$SLUDIST_LIBS $LIBS"
CPPFLAGS="$SLUDIST_INCLUDES $CPPFLAGS"

AC_CHECK_HEADER([superlu_ddefs.h],
 [pac_sludist_header_ok=yes],
 [pac_sludist_header_ok=no; SLUDIST_INCLUDES=""])
if test "x$pac_sludist_header_ok" == "xno" ; then 
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_superlu_ddefs_h
  SLUDIST_INCLUDES="-I$mld2p4_cv_superludistdir/include -I$mld2p4_cv_superludistdir/Include "
  CPPFLAGS="$SLUDIST_INCLUDES $SAVE_CPPFLAGS"

 AC_CHECK_HEADER([superlu_ddefs.h],
		 [pac_sludist_header_ok=yes],
		 [pac_sludist_header_ok=no; SLUDIST_INCLUDES=""])
fi

if test "x$pac_sludist_header_ok" == "xyes" ; then 
      SLUDIST_LIBS="$mld2p4_cv_superludist $SLUDIST_LIBS"
      LIBS="$SLUDIST_LIBS -lm $LIBS";
      AC_MSG_CHECKING([for superlu_malloc_dist in $SLUDIST_LIBS])
      AC_TRY_LINK_FUNC(superlu_malloc_dist, 
       [mld2p4_cv_have_superludist=yes;pac_sludist_lib_ok=yes;],
       [mld2p4_cv_have_superludist=no;pac_sludist_lib_ok=no; 
          SLUDIST_LIBS=""; ])
  if test "x$pac_sludist_lib_ok" == "xno" ; then 
     dnl Maybe lib?
     SLUDIST_LIBS="$mld2p4_cv_superludist  -L$mld2p4_cv_superludistdir/lib";
     LIBS="$SLUDIST_LIBS -lm $SAVE_LIBS";
     AC_TRY_LINK_FUNC(superlu_malloc_dist, 
		     [mld2p4_cv_have_superludist=yes;pac_sludist_lib_ok=yes;],
		     [mld2p4_cv_have_superludist=no;pac_sludist_lib_ok=no; 
		      SLUDIST_LIBS="";SLUDIST_INCLUDES=""])
 fi

 AC_MSG_RESULT($pac_sludist_lib_ok)
fi
 LIBS="$save_LIBS";
 CPPFLAGS="$save_CPPFLAGS";
 CC="$save_CC";
])dnl 

dnl @synopsis PAC_ARG_SERIAL_MPI
dnl
dnl Test for --with-serial-mpi={yes|no}
dnl 
dnl 
dnl
dnl Example use:
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_SERIAL_MPI],
[
AC_MSG_CHECKING([whether we want serial (fake) mpi])
AC_ARG_ENABLE(serial,
AC_HELP_STRING([--enable-serial], 
[Specify whether to enable a fake mpi library to run in serial mode. ]),
[
pac_cv_serial_mpi="yes";
]
dnl ,
dnl [pac_cv_serial_mpi="no";]
)
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   AC_MSG_RESULT([yes.])
else
 pac_cv_serial_mpi="no";
 AC_MSG_RESULT([no.])
fi
]
)
