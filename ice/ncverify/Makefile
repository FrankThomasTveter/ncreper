#! gmake
#
# Makefile    - template Makefile to compile and maintain 
#             a fortran object library and link binaries
#             using the C preprocessor (cpp) on fortran files
#             with the '.F' suffix
#
# Calling:
#	make		- make all binaries
#	make lib	- make library
#	make check_src	- checks if source files have been
#			  changed since last library update
#	make check_libs	- checks if external libaries have been
#			  changed since last binary update
#	make check	- performs both above checks
#
# + Get the list of All fortran source files (suffix .F) 
#   in the directory.
# + Build a library which contains all fortran objects.
#   (the objects are build with selected options)
# + Link binaries from the maintained and external libraries
#
# NOTE:	Names in <> are sympolic names to be replaced by the 
#	actual ones.
#	In [] are sample or recommented settings.
#
# Macro Description (to be adapted):
# MAINS		- list of fortran files containing main programs
#		  [<prog1.F>, <prog2.F>, ...]
# LIB_NAM	- full pathname of the library to be maintained
#		  [$(HOME)/lib/lib<proj>.a]
# BIN		- directory to place the binaries in
#		  [$(HOME)/bin]
# LIBS		- library options: 
#		  -L<lib1home> -l<lib1-name>
# INCL		- include options: 
#		  -I<includedir>
# OPTIONS - fortran compile & link options: 
#				   -g - create debugger info
#
#-------------------------------------------------------------------
#
# List of Main Programs
#	(Put a list of FORTRAN source files containing main programs
#	 here) 
#
 MAINS = ncverify.F90
#
# Full library name
#	(All objects will be stored in a library (using the UNIX 
#	 'ar' command) - put the full pathname of this library here)
#
# LIB_NAM = ../../lib/libbstat.a
#
# Target Directory for Tools Binaries
#	(All binaries will be installed in a dedicated 
#	 'bin' directory - put the path here)
#
#BIN = $(HOME)/bin
#
# Tools Libraries
#	(Specify libraries to be linked in here. (e.g. orbit/attitude libs
#	 MIND the ORDER!)
#
LIBS =/elysium/metop/retper/src/libfmgen4/libfmgen4.a \
      /elysium/metop/retper/src/libfmnuke/libfmnuke.a \
      /elysium/metop/retper/src/libfmmi/libfmmi.a

# LIBS= ../libfmscat/libfmscat.a \
#      ../../atovs/libfmgen4/libfmgen4.a \
#      ../../atovs/libfmnuke/libfmnuke.a \
#       ../../atovs/libfmmi/libfmmi.a
#       /vol/software/polsatproc/3party/lib/libemos.a
#
# LIBS= /home/polarsat/scat/qscat-v1p3/lib/libgen.a \
#       /home/polarsat/libs/lib/libfmgen4.a \
#       /home/polarsat/libs/lib/libfmnuke.a \
#       /home/polarsat/scat/qscat-v1p3/lib/libauxc2.a \
#       /home/polarsat/libs/libsatimg/lib/libsatimg.a \
#       /usr/lib/libtiff.a \
#       /usr/lib/libz.a \
#       /usr/lib/libjpeg.a \
#       /metno/local/lib/libmi_pgi.a 
#
#       ../../lib/libmisc.a
#
#       /metno/local/lib/libtiff.a \
#
# /raid/sat/libsat/libsatimg/lib/libsatimg.a
#
# Include Directory
#	(Specify directories to pick up include files (MIND the ORDER!)
#
 INCL = ./ /opt/netcdf-fortran-4.2/include/
#
# Compilation flags 
#	(Put the compilation flags here (-g is recommented))
#
 OPTIONS = -g -fbounds-check  -fconvert=big-endian -DPREC32
# OPTIONS = -g -fbounds-check
 BINOP =  -L/metno/local/lib -L/usr/lib32  -L/opt/netcdf-fortran-4.2/lib -lnetcdff -Wl,-rpath,/opt/netcdf-fortran-4.2/lib
 COMPILER = gfortran -fdefault-real-8 -fdefault-double-8
#------------------------------------------------------------------------
#		the rest not is not likely to be changed
#------------------------------------------------------------------------
#
# BIN = ./
 BIN_SUFFIX = .bin
#
 include MakeInclude_F90
#
