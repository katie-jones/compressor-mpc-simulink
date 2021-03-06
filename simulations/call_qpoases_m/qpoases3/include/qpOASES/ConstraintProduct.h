/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2014 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qpOASES/ConstraintProduct.h
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0embedded
 *	\date 2009-2014
 *
 *	Declaration of the ConstraintProduct class which allows to specify a
 *	user-defined function for evaluating the constraint product at the 
 *	current iterate to speed-up QP solution in case of a specially structured
 *	constraint matrix.
 */



#ifndef QPOASES_CONSTRAINT_PRODUCT_H
#define QPOASES_CONSTRAINT_PRODUCT_H


BEGIN_NAMESPACE_QPOASES


/** 
 *	\brief Interface for specifying user-defined evaluations of constraint products.
 *
 *	A class which allows to specify a user-defined function for evaluating the 
 *	constraint product at the current iterate to speed-up QP solution in case 
 *	of a specially structured constraint matrix.
 *
 *	\author Hans Joachim Ferreau
 *	\version 3.0embedded
 *	\date 2009-2012
 */
typedef struct 
{
	int empty;
} ConstraintProduct;


END_NAMESPACE_QPOASES

#endif	/* QPOASES_CONSTRAINT_PRODUCT_H */
