#
# mbuildopts.sh Shell script for configuring mbuild, the C/C++ stand alone 
#		applications creation script.
#
# usage:        Do not call this file directly; it is sourced by the
#               mbuild shell script.  Modify only if you don't like the
#               defaults after running mbuild.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Note: For the version of system compiler supported with this release,
#       refer to:
#       http://www.mathworks.com/support/compilers/current_release/
#
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_ML_OPT: Build and link with MATLAB Compiler generated library via the system ANSI C/C++ compiler
#
# Copyright 1984-2007 The MathWorks, Inc.
# $Revision: 1.51.4.20 $  $Date: 2009/03/09 18:27:37 $
#----------------------------------------------------------------------------
#
echo hello
htehkouhjju
    MLIBS="-L$TMW_ROOT/runtime/$Arch $OTHER_LIBS"
#   Link against the old set of libraries for HPUX only. For all other 
#   Unix platforms, link against mclmcrrt
    MLIBS="$MLIBS -lmwmclmcrrt"
    MFLAGS="-I$TMW_ROOT/extern/include -DUNIX -DX11"
    MCXXFLAGS="-I$TMW_ROOT/extern/include/cpp $MFLAGS"
    MCXXLIBS="$MLIBS"
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the TMW_ROOT
# root directory.  The mbuild script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            TMW_ROOT="$TMW_ROOT"
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            CC='gcc'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS -D_FILE_OFFSET_BITS=64" 
            CFLAGS="$CFLAGS $MFLAGS"
            CFLAGS="$CFLAGS -pthread"
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            SHLCFLAGS="-fPIC $CFLAGS"
#           
            CXX='g++'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -D_FILE_OFFSET_BITS=64" 
            CXXFLAGS="$CXXFLAGS $MCXXFLAGS -DGLNX86 -DGCC"
            CXXFLAGS="$CXXFLAGS -pthread"
            CXXLIBS="$RPATH $MCXXLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
            CXXLDFLAGS='-pthread'
#
            SHLCXXFLAGS="-fPIC $CXXFLAGS"
#
            LD="$COMPILER"
            LDFLAGS='-pthread'
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            SHLLD="$COMPILER"
            SHLLDFLAGS='-pthread -shared'
            SHLCXXLDFLAGS="$SHLLDFLAGS"
            SHLMAKEDEF='cat'
            SHLPOSTLINK_CMDS=':'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnxa64)
#----------------------------------------------------------------------------
	    RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            CC='gcc'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS $MFLAGS"
            CFLAGS="$CFLAGS -pthread"
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O2 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            SHLCFLAGS="-fPIC $CFLAGS"
#
            CXX='g++'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS $MCXXFLAGS -DGLNXA64 -DGCC"
            CXXFLAGS="$CXXFLAGS -pthread"
            CXXLIBS="$RPATH $MCXXLIBS -lm"
            CXXOPTIMFLAGS='-O2 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
            CXXLDFLAGS='-pthread'
#
            SHLCXXFLAGS="-fPIC $CXXFLAGS"
#
            LD="$COMPILER"
            LDFLAGS='-pthread'
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            SHLLD="$COMPILER"
            SHLLDFLAGS='-pthread -shared'
            SHLCXXLDFLAGS="$SHLLDFLAGS"
            SHLMAKEDEF='cat'
            SHLPOSTLINK_CMDS=':'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol64)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh sol64 12
#----------------------------------------------------------------------------
            ;;
        mac)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh mac 12
#----------------------------------------------------------------------------
            ;;
        maci)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh maci 12
#----------------------------------------------------------------------------
            ;;
        maci64)
#----------------------------------------------------------------------------
            CC='gcc-4.2'
            SDKROOT='/Developer/SDKs/MacOSX10.6.sdk'
            MACOSX_DEPLOYMENT_TARGET='10.5'
            ARCHS='x86_64'
            CFLAGS="-fno-common -no-cpp-precomp -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CFLAGS="$CFLAGS $MFLAGS"
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O2 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            SHLCFLAGS="$CFLAGS"
#           
            CXX=g++-4.2
            CXXFLAGS="-fno-common -no-cpp-precomp -fexceptions -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CXXFLAGS="$CXXFLAGS $MCXXFLAGS -DMACI64"
            CXXLIBS="$MCXXLIBS"
            CXXOPTIMFLAGS='-O2 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
            CXXLDFLAGS="-Wl,-twolevel_namespace -undefined error -arch $ARCHS -Wl,-syslibroot,$SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
#
            SHLCXXFLAGS="$CXXFLAGS"
#
            LD="$COMPILER"
            LDFLAGS="-Wl,-twolevel_namespace -undefined error -arch $ARCHS -Wl,-syslibroot,$SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            LDFLAGS="$LDFLAGS -bind_at_load"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            SHLLD='$COMPILER'
            SHLCXXLDFLAGS="-Wl,-twolevel_namespace -undefined error -dynamiclib -arch $ARCHS -Wl,-syslibroot,$SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            SHLLDFLAGS="$SHLCXXLDFLAGS -Wl,-exported_symbols_list,$EXPFILE"
            SHLMAKEDEF="awk '{printf \"_%s\\n\", \$0;}'"
            SHLPOSTLINK_CMDS=':'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
