##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2012 by Hans Joachim Ferreau, Andreas Potschka,
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
##	Filename:  src/Makefile
##	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
##	Version:   3.0embedded
##	Date:      2007-2014
##


include ../make.mk


##
##	flags
##

IDIR        =  ../include

IFLAGS      =  -I. \
               -I${IDIR}


QPOASES_OBJECTS = \
	QProblem.${OBJEXT} \
	QProblemB.${OBJEXT} \
	Bounds.${OBJEXT} \
	Constraints.${OBJEXT} \
	Indexlist.${OBJEXT} \
	Utils.${OBJEXT} \
	Options.${OBJEXT} \
	Matrices.${OBJEXT} \
	MessageHandling.${OBJEXT}

QPOASES_DEPENDS = \
	${IDIR}/qpOASES/QProblem.h \
	${IDIR}/qpOASES/QProblemB.h \
	${IDIR}/qpOASES/Bounds.h \
	${IDIR}/qpOASES/Constraints.h \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/Matrices.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp


##
##	targets
##

all: ${BINDIR}/libqpOASES.${LIBEXT}


${BINDIR}/libqpOASES.${LIBEXT}: ${QPOASES_OBJECTS}
	${AR} r $@ $?


clean:
	@${ECHO} "Cleaning up (src)"
	@${RM} -f *.${OBJEXT} *.${LIBEXT}

clobber: clean


%.${OBJEXT}: %.c
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} $<

shared:	${QPOASES_OBJECTS}
	${CPP} -shared -Wl,-soname,libqpOASES.so -o libqpOASES.so ${QPOASES_OBJECTS} -lm 

##
##	end of file
##
