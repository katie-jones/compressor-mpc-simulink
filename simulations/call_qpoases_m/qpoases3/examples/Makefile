##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2014 by Hans Joachim Ferreau, Andreas Potschka,
##	Christian Kirches et al. All rights reserved.
##
##	qpOASES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpOASES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpOASES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



##
##	Filename:  examples/Makefile
##	Author:    Hans Joachim Ferreau
##	Version:   3.0embedded
##	Date:      2007-2014
##



include ../make.mk


##
##	flags
##

IDIR        =  ../include

IFLAGS      =  -I. \
               -I${IDIR} \
			   -I../testing/

QPOASES_EXES = \
	${BINDIR}/example1${EXE} \
	${BINDIR}/example1b${EXE} \
	${BINDIR}/exampleLP${EXE} \
	${BINDIR}/example3${EXE}


##
##	targets
##

all: ${QPOASES_EXES}


${BINDIR}/example1${EXE}: example1.${OBJEXT}
	@${ECHO} "Creating" $@
	${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/example1b${EXE}: example1b.${OBJEXT}
	@${ECHO} "Creating" $@
	${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/exampleLP${EXE}: exampleLP.${OBJEXT}
	@${ECHO} "Creating" $@
	${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/example3${EXE}: example3.${OBJEXT} ../testing/OQPinterface.${OBJEXT}
	@${ECHO} "Creating" $@
	${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ../testing/OQPinterface.${OBJEXT} ${QPOASES_LINK} ${LINK_LIBRARIES}

clean:
	${RM} -f *.${OBJEXT} ${QPOASES_EXES}

clobber: clean


../testing/OQPinterface.${OBJEXT}: \
	../testing/OQPinterface.c \
	../testing/OQPinterface.h \
	${IDIR}/qpOASES/QProblem.h \
	${IDIR}/qpOASES/QProblemB.h \
	${IDIR}/qpOASES/Utils.h \
	${IDIR}/qpOASES/Constants.h \
	${IDIR}/qpOASES/Types.h \
	${IDIR}/qpOASES/Options.h \
	${IDIR}/qpOASES/Matrices.h \
	${IDIR}/qpOASES/MessageHandling.h
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} ../testing/OQPinterface.c

%.${OBJEXT}: %.c
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} $<


##
##	end of file
##
