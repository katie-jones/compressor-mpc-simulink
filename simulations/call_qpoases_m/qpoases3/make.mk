TOP = $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
BINDIR = ${TOP}/bin
#include ${TOP}/make_linux.mk
include ${TOP}/make_cygwin.mk
#include ${TOP}/make_windows.mk
#include ${TOP}/make_osx.mk
