#ifndef _INCLUDE_BASE_LIBMESH_CONFIG_H
#define _INCLUDE_BASE_LIBMESH_CONFIG_H 1
 
/* include/base/libmesh_config.h. Generated automatically at end of configure. */
/* include/base/libmesh_config.h.tmp.  Generated from libmesh_config.h.in by configure.  */
/* include/base/libmesh_config.h.in.  Generated from configure.in by autoheader.  */

/* definition of the final detected unordered_map type */
#ifndef LIBMESH_BEST_UNORDERED_MAP
#define LIBMESH_BEST_UNORDERED_MAP std::tr1::unordered_map
#endif

/* definition of the final detected unordered_multimap type */
#ifndef LIBMESH_BEST_UNORDERED_MULTIMAP
#define LIBMESH_BEST_UNORDERED_MULTIMAP std::tr1::unordered_multimap
#endif

/* definition of the final detected unordered_set type */
#ifndef LIBMESH_BEST_UNORDERED_SET
#define LIBMESH_BEST_UNORDERED_SET std::tr1::unordered_set
#endif

/* This compiler is known not to support some iostream functionality */
/* #undef BROKEN_IOSTREAM */

/* Configuration information. */
#ifndef LIBMESH_CONFIGURE_INFO
#define LIBMESH_CONFIGURE_INFO "./configure run on Fri Apr 15 15:11:44 EDT 2016"
#endif

/* Flag indicating if double-precision (double) should be used for most
   floating-point calculations */
#ifndef LIBMESH_DEFAULT_DOUBLE_PRECISION
#define LIBMESH_DEFAULT_DOUBLE_PRECISION 1
#endif

/* Flag indicating if quadruple-precision (__float128) should be used for most
   floating-point calculations */
/* #undef DEFAULT_QUADRUPLE_PRECISION */

/* Data type to be used for most floating-point calculations */
#ifndef LIBMESH_DEFAULT_SCALAR_TYPE
#define LIBMESH_DEFAULT_SCALAR_TYPE double
#endif

/* Flag indicating if single-precision (float) should be used for most
   floating-point calculations */
/* #undef DEFAULT_SINGLE_PRECISION */

/* Flag indicating if triple-precision (long double) should be used for most
   floating-point calculations */
/* #undef DEFAULT_TRIPLE_PRECISION */

/* workaround for potentially missing hash<T*> */
#ifndef LIBMESH_DEFINE_HASH_POINTERS
#define LIBMESH_DEFINE_HASH_POINTERS /**/
#endif

/* workaround for potentially missing hash<string> */
#ifndef LIBMESH_DEFINE_HASH_STRING
#define LIBMESH_DEFINE_HASH_STRING /**/
#endif

/* PETSc's major version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_MAJOR
#define LIBMESH_DETECTED_PETSC_VERSION_MAJOR 3
#endif

/* PETSc's minor version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_MINOR
#define LIBMESH_DETECTED_PETSC_VERSION_MINOR 2
#endif

/* PETSc's subminor version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR
#define LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR 0
#endif

/* Integer indicating the highest spatial dimensionality supported by libMesh
   */
/* #undef DIM */

/* Flag indicating if the library should be built with AMR support */
#ifndef LIBMESH_ENABLE_AMR
#define LIBMESH_ENABLE_AMR 1
#endif

/* Flag indicating if the library should be built with Dirichlet boundary
   constraint support */
#ifndef LIBMESH_ENABLE_DIRICHLET
#define LIBMESH_ENABLE_DIRICHLET 1
#endif

/* Flag indicating if the library should be built to throw C++ exceptions on
   unexpected errors */
#ifndef LIBMESH_ENABLE_EXCEPTIONS
#define LIBMESH_ENABLE_EXCEPTIONS 1
#endif

/* Flag indicating if the library should use ghosted local vectors */
#ifndef LIBMESH_ENABLE_GHOSTED
#define LIBMESH_ENABLE_GHOSTED 1
#endif

/* Flag indicating if the library should offer higher order p-FEM shapes */
#ifndef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
#define LIBMESH_ENABLE_HIGHER_ORDER_SHAPES 1
#endif

/* Flag indicating if the library should be built with infinite elements */
/* #undef ENABLE_INFINITE_ELEMENTS */

/* Flag indicating if the library should be built with node constraints
   support */
/* #undef ENABLE_NODE_CONSTRAINTS */

/* Flag indicating if the library should use the experimental ParallelMesh as
   its default Mesh type */
/* #undef ENABLE_PARMESH */

/* Flag indicating if the library should be built with performance logging
   support */
/* #undef ENABLE_PERFORMANCE_LOGGING */

/* Flag indicating if the library should be built with periodic boundary
   condition support */
#ifndef LIBMESH_ENABLE_PERIODIC
#define LIBMESH_ENABLE_PERIODIC 1
#endif

/* Flag indicating if the library should be built with reference counting
   support */
#ifndef LIBMESH_ENABLE_REFERENCE_COUNTING
#define LIBMESH_ENABLE_REFERENCE_COUNTING 1
#endif

/* Flag indicating if the library should be built with second derivatives */
#ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
#define LIBMESH_ENABLE_SECOND_DERIVATIVES 1
#endif

/* Flag indicating if the library should be built to write stack trace files
   on unexpected errors */
/* #undef ENABLE_TRACEFILES */

/* Flag indicating if the library should be built with variational smoother
   support */
#ifndef LIBMESH_ENABLE_VSMOOTHER
#define LIBMESH_ENABLE_VSMOOTHER 1
#endif

/* Flag indicating whether the library shall be compiled to use the Trilinos
   solver collection */
/* #undef HAVE_AZTECOO */

/* define if the Boost library is available */
#ifndef LIBMESH_HAVE_BOOST
#define LIBMESH_HAVE_BOOST /**/
#endif

/* Flag indicating bzip2/bunzip2 are available for handling compressed .bz2
   files */
#ifndef LIBMESH_HAVE_BZIP
#define LIBMESH_HAVE_BZIP 1
#endif

/* Define to 1 if you have the <csignal> header file. */
#ifndef LIBMESH_HAVE_CSIGNAL
#define LIBMESH_HAVE_CSIGNAL 1
#endif

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef LIBMESH_HAVE_DLFCN_H
#define LIBMESH_HAVE_DLFCN_H 1
#endif

/* Flag indicating whether the library will be compiled with Eigen support */
/* #undef HAVE_EIGEN */

/* Flag indicating whether the library will be compiled with Exodus support */
#ifndef LIBMESH_HAVE_EXODUS_API
#define LIBMESH_HAVE_EXODUS_API 1
#endif

/* define if the compiler supports __gnu_cxx::hash_map */
/* #undef HAVE_EXT_HASH_MAP */

/* define if the compiler supports __gnu_cxx::hash_multimap */
/* #undef HAVE_EXT_HASH_MULTIMAP */

/* define if the compiler supports __gnu_cxx::hash_set */
/* #undef HAVE_EXT_HASH_SET */

/* Define to 1 if you have the <fenv.h> header file. */
#ifndef LIBMESH_HAVE_FENV_H
#define LIBMESH_HAVE_FENV_H 1
#endif

/* Flag indicating whether the library will be compiled with FPARSER support
   */
/* #undef HAVE_FPARSER */

/* define if the compiler supports GCC C++ ABI name demangling */
#ifndef LIBMESH_HAVE_GCC_ABI_DEMANGLE
#define LIBMESH_HAVE_GCC_ABI_DEMANGLE 1
#endif

/* Define to 1 if you have the <getopt.h> header file. */
#ifndef LIBMESH_HAVE_GETOPT_H
#define LIBMESH_HAVE_GETOPT_H 1
#endif

/* Flag indicating if the library should be built with calls to getpwuid() */
#ifndef LIBMESH_HAVE_GETPWUID
#define LIBMESH_HAVE_GETPWUID 1
#endif

/* define if the compiler supports glibc backtrace */
#ifndef LIBMESH_HAVE_GLIBC_BACKTRACE
#define LIBMESH_HAVE_GLIBC_BACKTRACE 1
#endif

/* Flag indicating whether the library will be compiled with GLPK support */
/* #undef HAVE_GLPK */

/* Flag indicating whether the library will be compiled with GMV support */
#ifndef LIBMESH_HAVE_GMV
#define LIBMESH_HAVE_GMV 1
#endif

/* Flag indicating whether or not gzstreams are available */
#ifndef LIBMESH_HAVE_GZSTREAM
#define LIBMESH_HAVE_GZSTREAM 1
#endif

/* define if the compiler supports std::hash_map */
/* #undef HAVE_HASH_MAP */

/* define if the compiler supports std::hash_multimap */
/* #undef HAVE_HASH_MULTIMAP */

/* define if the compiler supports std::hash_set */
/* #undef HAVE_HASH_SET */

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef LIBMESH_HAVE_INTTYPES_H
#define LIBMESH_HAVE_INTTYPES_H 1
#endif

/* Flag indicating whether the library will be compiled with LASPACK support
   */
#ifndef LIBMESH_HAVE_LASPACK
#define LIBMESH_HAVE_LASPACK 1
#endif

/* Flag indicating whether the library will be compiled with libHilbert
   support */
#ifndef LIBMESH_HAVE_LIBHILBERT
#define LIBMESH_HAVE_LIBHILBERT 1
#endif

/* define if the compiler has locale */
#ifndef LIBMESH_HAVE_LOCALE
#define LIBMESH_HAVE_LOCALE /**/
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef LIBMESH_HAVE_MEMORY_H
#define LIBMESH_HAVE_MEMORY_H 1
#endif

/* Flag indicating whether the library will be compiled with Metis support */
#ifndef LIBMESH_HAVE_METIS
#define LIBMESH_HAVE_METIS 1
#endif

/* Flag indicating whether the library shall be compiled to use the ML package
   */
/* #undef HAVE_ML */

/* Flag indicating whether or not MPI is available */
#ifndef LIBMESH_HAVE_MPI
#define LIBMESH_HAVE_MPI 1
#endif

/* define if the compiler implements namespaces */
#ifndef LIBMESH_HAVE_NAMESPACES
#define LIBMESH_HAVE_NAMESPACES /**/
#endif

/* Flag indicating whether the library will be compiled with Nemesis support
   */
#ifndef LIBMESH_HAVE_NEMESIS_API
#define LIBMESH_HAVE_NEMESIS_API 1
#endif

/* Flag indicating whether the library will be compiled with Netcdf support */
#ifndef LIBMESH_HAVE_NETCDF
#define LIBMESH_HAVE_NETCDF 1
#endif

/* Flag indicating whether the library shall be compiled to use the Nox solver
   collection */
/* #undef HAVE_NOX */

/* Define if OpenMP is enabled */
/* #undef HAVE_OPENMP */

/* Flag indicating whether the library will be compiled with Parmetis support
   */
#ifndef LIBMESH_HAVE_PARMETIS
#define LIBMESH_HAVE_PARMETIS 1
#endif

/* Flag indicating whether or not PETSc is available */
#ifndef LIBMESH_HAVE_PETSC
#define LIBMESH_HAVE_PETSC 1
#endif

/* Flag indicating whether or not PETSc was compiled with Hypre support */
/* #undef HAVE_PETSC_HYPRE */

/* Define to 1 if you have the <rpc/rpc.h> header file. */
#ifndef LIBMESH_HAVE_RPC_RPC_H
#define LIBMESH_HAVE_RPC_RPC_H 1
#endif

/* define if the compiler supports Run-Time Type Identification */
#ifndef LIBMESH_HAVE_RTTI
#define LIBMESH_HAVE_RTTI /**/
#endif

/* Flag indicating whether the library will be compiled with SFC support */
#ifndef LIBMESH_HAVE_SFCURVES
#define LIBMESH_HAVE_SFCURVES 1
#endif

/* Flag indicating whether or not SLEPc is available */
/* #undef HAVE_SLEPC */

/* define if the compiler has the sstream header */
#ifndef LIBMESH_HAVE_SSTREAM
#define LIBMESH_HAVE_SSTREAM /**/
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef LIBMESH_HAVE_STDINT_H
#define LIBMESH_HAVE_STDINT_H 1
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef LIBMESH_HAVE_STDLIB_H
#define LIBMESH_HAVE_STDLIB_H 1
#endif

/* define if the compiler supports std::unordered_map */
/* #undef HAVE_STD_UNORDERED_MAP */

/* define if the compiler supports std::unordered_multimap */
/* #undef HAVE_STD_UNORDERED_MULTIMAP */

/* define if the compiler supports std::unordered_set */
/* #undef HAVE_STD_UNORDERED_SET */

/* define if the compiler has stringstream functionality */
#ifndef LIBMESH_HAVE_STRINGSTREAM
#define LIBMESH_HAVE_STRINGSTREAM /**/
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef LIBMESH_HAVE_STRINGS_H
#define LIBMESH_HAVE_STRINGS_H 1
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef LIBMESH_HAVE_STRING_H
#define LIBMESH_HAVE_STRING_H 1
#endif

/* define if the compiler has the strstream header */
/* #undef HAVE_STRSTREAM */

/* Define to 1 if you have the <sys/resource.h> header file. */
#ifndef LIBMESH_HAVE_SYS_RESOURCE_H
#define LIBMESH_HAVE_SYS_RESOURCE_H 1
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef LIBMESH_HAVE_SYS_STAT_H
#define LIBMESH_HAVE_SYS_STAT_H 1
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef LIBMESH_HAVE_SYS_TYPES_H
#define LIBMESH_HAVE_SYS_TYPES_H 1
#endif

/* Flag indicating whether the library shall be compiled to use the Threading
   Building Blocks */
/* #undef HAVE_TBB_API */

/* Flag indicating whether the library shall be compiled to use the Tecplot
   interface */
#ifndef LIBMESH_HAVE_TECPLOT_API
#define LIBMESH_HAVE_TECPLOT_API 1
#endif

/* Flag indicating tecplot API understands newer features */
#ifndef LIBMESH_HAVE_TECPLOT_API_112
#define LIBMESH_HAVE_TECPLOT_API_112 1
#endif

/* Flag indicating whether the library will be compiled with Tetgen support */
#ifndef LIBMESH_HAVE_TETGEN
#define LIBMESH_HAVE_TETGEN 1
#endif

/* define if the compiler supports std::tr1::unordered_map */
#ifndef LIBMESH_HAVE_TR1_UNORDERED_MAP
#define LIBMESH_HAVE_TR1_UNORDERED_MAP 1
#endif

/* define if the compiler supports std::tr1::unordered_multimap */
#ifndef LIBMESH_HAVE_TR1_UNORDERED_MULTIMAP
#define LIBMESH_HAVE_TR1_UNORDERED_MULTIMAP 1
#endif

/* define if the compiler supports std::tr1::unordered_set */
#ifndef LIBMESH_HAVE_TR1_UNORDERED_SET
#define LIBMESH_HAVE_TR1_UNORDERED_SET 1
#endif

/* Flag indicating whether the library will be compiled with Triangle support
   */
#ifndef LIBMESH_HAVE_TRIANGLE
#define LIBMESH_HAVE_TRIANGLE 1
#endif

/* Flag indicating whether the library shall be compiled to use the Trilinos
   solver collection */
/* #undef HAVE_TRILINOS */

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef LIBMESH_HAVE_UNISTD_H
#define LIBMESH_HAVE_UNISTD_H 1
#endif

/* Flag indicating whether the library will be compiled with VTK support */
#ifndef LIBMESH_HAVE_VTK
#define LIBMESH_HAVE_VTK 1
#endif

/* Flag indicating headers and libraries for XDR IO are available */
#ifndef LIBMESH_HAVE_XDR
#define LIBMESH_HAVE_XDR 1
#endif

/* Define to 1 if you have the <xmmintrin.h> header file. */
#ifndef LIBMESH_HAVE_XMMINTRIN_H
#define LIBMESH_HAVE_XMMINTRIN_H 1
#endif

/* Flag indicating xz is available for handling compressed .xz files */
#ifndef LIBMESH_HAVE_XZ
#define LIBMESH_HAVE_XZ 1
#endif

/* Define to 1 if you have the <zlib.h> header file. */
#ifndef LIBMESH_HAVE_ZLIB_H
#define LIBMESH_HAVE_ZLIB_H 1
#endif

/* header file for the final detected unordered_map type */
#ifndef LIBMESH_INCLUDE_UNORDERED_MAP
#define LIBMESH_INCLUDE_UNORDERED_MAP <tr1/unordered_map>
#endif

/* header file for the final detected unordered_multimap type */
#ifndef LIBMESH_INCLUDE_UNORDERED_MULTIMAP
#define LIBMESH_INCLUDE_UNORDERED_MULTIMAP <tr1/unordered_map>
#endif

/* header file for the final detected unordered_set type */
#ifndef LIBMESH_INCLUDE_UNORDERED_SET
#define LIBMESH_INCLUDE_UNORDERED_SET <tr1/unordered_set>
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef LIBMESH_PACKAGE_BUGREPORT
#define LIBMESH_PACKAGE_BUGREPORT ""
#endif

/* Define to the full name of this package. */
#ifndef LIBMESH_PACKAGE_NAME
#define LIBMESH_PACKAGE_NAME ""
#endif

/* Define to the full name and version of this package. */
#ifndef LIBMESH_PACKAGE_STRING
#define LIBMESH_PACKAGE_STRING ""
#endif

/* Define to the one symbol short name of this package. */
#ifndef LIBMESH_PACKAGE_TARNAME
#define LIBMESH_PACKAGE_TARNAME ""
#endif

/* Define to the home page for this package. */
#ifndef LIBMESH_PACKAGE_URL
#define LIBMESH_PACKAGE_URL ""
#endif

/* Define to the version of this package. */
#ifndef LIBMESH_PACKAGE_VERSION
#define LIBMESH_PACKAGE_VERSION ""
#endif

/* The size of `double', as computed by sizeof. */
#ifndef LIBMESH_SIZEOF_DOUBLE
#define LIBMESH_SIZEOF_DOUBLE 8
#endif

/* The size of `float', as computed by sizeof. */
#ifndef LIBMESH_SIZEOF_FLOAT
#define LIBMESH_SIZEOF_FLOAT 4
#endif

/* The size of `function_pointer', as computed by sizeof. */
#ifndef LIBMESH_SIZEOF_FUNCTION_POINTER
#define LIBMESH_SIZEOF_FUNCTION_POINTER 8
#endif

/* The size of `int', as computed by sizeof. */
#ifndef LIBMESH_SIZEOF_INT
#define LIBMESH_SIZEOF_INT 4
#endif

/* The size of `long int', as computed by sizeof. */
#ifndef LIBMESH_SIZEOF_LONG_INT
#define LIBMESH_SIZEOF_LONG_INT 8
#endif

/* The size of `short int', as computed by sizeof. */
#ifndef LIBMESH_SIZEOF_SHORT_INT
#define LIBMESH_SIZEOF_SHORT_INT 2
#endif

/* The size of `unsigned int', as computed by sizeof. */
#ifndef LIBMESH_SIZEOF_UNSIGNED_INT
#define LIBMESH_SIZEOF_UNSIGNED_INT 4
#endif

/* The size of `void *', as computed by sizeof. */
#ifndef LIBMESH_SIZEOF_VOID_P
#define LIBMESH_SIZEOF_VOID_P 8
#endif

/* Define to 1 if you have the ANSI C header files. */
#ifndef LIBMESH_STDC_HEADERS
#define LIBMESH_STDC_HEADERS 1
#endif

/* size of subdomain_id */
#ifndef LIBMESH_SUBDOMAIN_ID_BYTES
#define LIBMESH_SUBDOMAIN_ID_BYTES 2
#endif

/* If the compiler supports a TLS storage class define it to that here */
#ifndef LIBMESH_TLS
#define LIBMESH_TLS __thread
#endif

/* Flag indicating if the library should be built using complex numbers */
/* #undef USE_COMPLEX_NUMBERS */

/* Flag indicating if the library should be built using real numbers */
#ifndef LIBMESH_USE_REAL_NUMBERS
#define LIBMESH_USE_REAL_NUMBERS 1
#endif

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#ifndef _libmesh_restrict
#define _libmesh_restrict __restrict
#endif
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif
 
/* once: _INCLUDE_BASE_LIBMESH_CONFIG_H */
#endif
