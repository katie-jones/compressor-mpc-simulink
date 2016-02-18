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
 *	\file src/Constraints.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0embedded
 *	\date 2007-2014
 *
 *	Implementation of the Constraints class designed to manage working sets of
 *	constraints within a QProblem.
 */


#include <qpOASES/Constraints.h>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	C o n s t r a i n t s
 */
void ConstraintsCON(	Constraints* THIS,
						int _n
						)
{
	Constraints_init( THIS,_n );
}


/*
 *	c o p y
 */
void ConstraintsCPY(	Constraints* FROM,
						Constraints* TO
						)
{
	int i;

	TO->n = FROM->n;
	TO->noLower = FROM->noLower;
	TO->noUpper = FROM->noUpper;

	if ( FROM->n != 0 )
	{
		for( i=0; i<TO->n; ++i )
		{
			TO->type[i]   = FROM->type[i];
			TO->status[i] = FROM->status[i];
		}
	}

	IndexlistCPY( &(FROM->active),&(TO->active) );
	IndexlistCPY( &(FROM->inactive),&(TO->inactive) );
}



/*
 *	i n i t
 */
returnValue Constraints_init(	Constraints* THIS,
								int _n
								)
{
	int i;

	if ( _n < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( _n > 0 )
	{
		Indexlist_init( &(THIS->active),_n );
		Indexlist_init( &(THIS->inactive),_n );
	}

	THIS->n = _n;
	THIS->noLower = BT_TRUE;
	THIS->noUpper = BT_TRUE;

	assert( THIS->n <= NCMAX );

	if ( THIS->n > 0 )
	{
		for( i=0; i<THIS->n; ++i )
		{
			THIS->type[i]   = ST_UNKNOWN;
			THIS->status[i] = ST_UNDEFINED;
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p C o n s t r a i n t
 */
returnValue Constraints_setupConstraint(	Constraints* THIS,
											int number, SubjectToStatus _status
											)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Add constraint index to respective index list. */
	switch ( _status )
	{
		case ST_INACTIVE:
			if ( Constraints_addIndex( THIS,Constraints_getInactive( THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
			break;

		case ST_LOWER:
			if ( Constraints_addIndex( THIS,Constraints_getActive( THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
			break;

		case ST_UPPER:
			if ( Constraints_addIndex( THIS,Constraints_getActive( THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
			break;

		default:
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A l l I n a c t i v e
 */
returnValue Constraints_setupAllInactive(	Constraints* THIS
											)
{
	return Constraints_setupAll( THIS,ST_INACTIVE );
}


/*
 *	s e t u p A l l L o w e r
 */
returnValue Constraints_setupAllLower(	Constraints* THIS
										)
{
	return Constraints_setupAll( THIS,ST_LOWER );
}


/*
 *	s e t u p A l l U p p e r
 */
returnValue Constraints_setupAllUpper(	Constraints* THIS
										)
{
	return Constraints_setupAll( THIS,ST_UPPER );
}


/*
 *	m o v e A c t i v e T o I n a c t i v e
 */
returnValue Constraints_moveActiveToInactive(	Constraints* THIS,
												int number
												)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of active constraints to that of inactive ones. */
	if ( Constraints_removeIndex( THIS,Constraints_getActive( THIS ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( Constraints_addIndex( THIS,Constraints_getInactive( THIS ),number,ST_INACTIVE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	m o v e I n a c t i v e T o A c t i v e
 */
returnValue Constraints_moveInactiveToActive(	Constraints* THIS,
												int number, SubjectToStatus _status
												)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of inactive constraints to that of active ones. */
	if ( Constraints_removeIndex( THIS,Constraints_getInactive( THIS ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( Constraints_addIndex( THIS,Constraints_getActive( THIS ),number,_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	f l i p F i x e d
 */
returnValue Constraints_flipFixed( Constraints* THIS, int number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( THIS->status != 0 )
		switch (THIS->status[number])
		{
			case ST_LOWER: THIS->status[number] = ST_UPPER; break;
			case ST_UPPER: THIS->status[number] = ST_LOWER; break;
			default: return THROWERROR( RET_MOVING_CONSTRAINT_FAILED );
		}

	return SUCCESSFUL_RETURN;
}


/*
 *	s h i f t
 */
returnValue Constraints_shift( Constraints* THIS, int offset )
{
	int i;
	static Indexlist shiftedActive;
	static Indexlist shiftedInactive;

	Indexlist_init( &shiftedActive,THIS->n );
	Indexlist_init( &shiftedInactive,THIS->n );

	/* consistency check */
	if ( ( offset == 0 ) || ( THIS->n <= 1 ) )
		return SUCCESSFUL_RETURN;

	if ( ( offset < 0 ) || ( offset > THIS->n/2 ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( ( THIS->n % offset ) != 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Shift types and status. */
	for( i=0; i<THIS->n-offset; ++i )
	{
		Constraints_setType( THIS,i,Constraints_getType( THIS,i+offset ) );
		Constraints_setStatus( THIS,i,Constraints_getStatus( THIS,i+offset ) );
	}

	/* 2) Construct shifted index lists of free and fixed variables. */
	for( i=0; i<THIS->n; ++i )
	{
		switch ( Constraints_getStatus( THIS,i ) )
		{
			case ST_INACTIVE:
				if ( Indexlist_addNumber( &shiftedInactive,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_LOWER:
				if ( Indexlist_addNumber( &shiftedActive,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_UPPER:
				if ( Indexlist_addNumber( &shiftedActive,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			default:
				return THROWERROR( RET_SHIFTING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	IndexlistCPY( &shiftedActive,&(THIS->active) );
	IndexlistCPY( &shiftedInactive,&(THIS->inactive) );

	return SUCCESSFUL_RETURN;
}


/*
 *	r o t a t e
 */
returnValue Constraints_rotate( Constraints* THIS, int offset )
{
	int i;
	myStatic SubjectToType   typeTmp[NCMAX];
	myStatic SubjectToStatus statusTmp[NCMAX];
	
	myStatic Indexlist rotatedActive;
	myStatic Indexlist rotatedInactive;

	Indexlist_init( &rotatedActive,THIS->n );
	Indexlist_init( &rotatedInactive,THIS->n );
	
	/* consistency check */
	if ( ( offset == 0 ) || ( offset == THIS->n ) || ( THIS->n <= 1 ) )
		return SUCCESSFUL_RETURN;

	if ( ( offset < 0 ) || ( offset > THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );


	/* 1) Rotate types and status. */
	for( i=0; i<offset; ++i )
	{
		typeTmp[i] = Constraints_getType( THIS,i );
		statusTmp[i] = Constraints_getStatus( THIS,i );
	}

	for( i=0; i<THIS->n-offset; ++i )
	{
		Constraints_setType( THIS,i,Constraints_getType( THIS,i+offset ) );
		Constraints_setStatus( THIS,i,Constraints_getStatus( THIS,i+offset ) );
	}

	for( i=THIS->n-offset; i<THIS->n; ++i )
	{
		Constraints_setType( THIS,i,typeTmp[i-THIS->n+offset] );
		Constraints_setStatus( THIS,i,statusTmp[i-THIS->n+offset] );
	}

	/* 2) Construct shifted index lists of free and fixed variables. */
	for( i=0; i<THIS->n; ++i )
	{
		switch ( Constraints_getStatus( THIS,i ) )
		{
			case ST_INACTIVE:
				if ( Indexlist_addNumber( &rotatedInactive,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_LOWER:
				if ( Indexlist_addNumber( &rotatedActive,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_UPPER:
				if ( Indexlist_addNumber( &rotatedActive,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			default:
				return THROWERROR( RET_ROTATING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	IndexlistCPY( &rotatedActive,&(THIS->active) );
	IndexlistCPY( &rotatedInactive,&(THIS->inactive) );

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue Constraints_print( Constraints* THIS )
{
	#ifndef __XPCTARGET__
	#ifndef __DSPACE__
	static char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	int nIAC = Constraints_getNIAC( THIS );
	int nAC  = Constraints_getNAC( THIS );

	int *IAC_idx, *AC_idx;

	if ( THIS->n == 0 )
		return SUCCESSFUL_RETURN;

	Indexlist_getNumberArray( Constraints_getInactive( THIS ),&IAC_idx );
	Indexlist_getNumberArray( Constraints_getActive( THIS ),&AC_idx );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"Constraints object comprising %d constraints (%d inactive, %d active):\n",THIS->n,nIAC,nAC );
	qpOASES_myPrintf( myPrintfString );

	REFER_NAMESPACE_QPOASES qpOASES_printNI( IAC_idx,nIAC,"inactive" );
	REFER_NAMESPACE_QPOASES qpOASES_printNI( AC_idx, nAC, "active  " );

	#endif
	#endif

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	s e t u p A l l
 */
returnValue Constraints_setupAll( Constraints* THIS, SubjectToStatus _status )
{
	int i;

	/* 1) Place unbounded constraints at the beginning of the index list of inactive constraints. */
	for( i=0; i<THIS->n; ++i )
	{
		if ( Constraints_getType( THIS,i ) == ST_UNBOUNDED )
		{
			if ( Constraints_setupConstraint( THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
		}
	}

	/* 2) Add remaining (i.e. "real" inequality) constraints to the index list of inactive constraints. */
	for( i=0; i<THIS->n; ++i )
	{
		if ( Constraints_getType( THIS,i ) == ST_BOUNDED )
		{
			if ( Constraints_setupConstraint( THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
		}
	}

	/* 3) Place implicit equality constraints at the end of the index list of inactive constraints. */
	for( i=0; i<THIS->n; ++i )
	{
		if ( Constraints_getType( THIS,i ) == ST_EQUALITY )
		{
			if ( Constraints_setupConstraint( THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
		}
	}

	/* 4) Moreover, add all constraints of unknown type. */
	for( i=0; i<THIS->n; ++i )
	{
		if ( Constraints_getType( THIS,i ) == ST_UNKNOWN || Constraints_getType( THIS,i ) == ST_DISABLED )
		{
			if ( Constraints_setupConstraint( THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
		}
	}


	return SUCCESSFUL_RETURN;
}


/*
 *	a d d I n d e x
 */
returnValue Constraints_addIndex(	Constraints* THIS,
									Indexlist* const indexlist,
									int newnumber, SubjectToStatus newstatus
									)
{
	if ( THIS->status != 0 )
	{
		/* consistency check */
		if ( THIS->status[newnumber] == newstatus )
			return THROWERROR( RET_INDEX_ALREADY_OF_DESIRED_STATUS );

		THIS->status[newnumber] = newstatus;
	}
	else
		return THROWERROR( RET_ADDINDEX_FAILED );

	if ( indexlist != 0 )
	{
		if ( Indexlist_addNumber( indexlist,newnumber ) == RET_INDEXLIST_EXCEEDS_MAX_LENGTH )
			return THROWERROR( RET_ADDINDEX_FAILED );
	}
	else
		return THROWERROR( RET_INVALID_ARGUMENTS );

	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e I n d e x
 */
returnValue Constraints_removeIndex(	Constraints* THIS,
										Indexlist* const indexlist,
										int removenumber
										)
{
	if ( THIS->status != 0 )
		THIS->status[removenumber] = ST_UNDEFINED;
	else
		return THROWERROR( RET_REMOVEINDEX_FAILED );

	if ( indexlist != 0 )
	{
		if ( Indexlist_removeNumber( indexlist,removenumber ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_REMOVEINDEX_FAILED );
	}
	else
		return THROWERROR( RET_INVALID_ARGUMENTS );

	return SUCCESSFUL_RETURN;
}


/*
 *	s w a p I n d e x
 */
returnValue Constraints_swapIndex(	Constraints* THIS, Indexlist* const indexlist,
									int number1, int number2
									)
{
	/* consistency checks */
	if ( THIS->status != 0 )
	{
		if ( THIS->status[number1] != THIS->status[number2] )
			return THROWERROR( RET_SWAPINDEX_FAILED );
	}
	else
		return THROWERROR( RET_SWAPINDEX_FAILED );

	if ( number1 == number2 )
	{
		THROWWARNING( RET_NOTHING_TO_DO );
		return SUCCESSFUL_RETURN;
	}

	if ( indexlist != 0 )
	{
		if ( Indexlist_swapNumbers( indexlist,number1,number2 ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SWAPINDEX_FAILED );
	}
	else
		return THROWERROR( RET_INVALID_ARGUMENTS );

	return SUCCESSFUL_RETURN;
}



END_NAMESPACE_QPOASES


/*
 *	end of file
 */
