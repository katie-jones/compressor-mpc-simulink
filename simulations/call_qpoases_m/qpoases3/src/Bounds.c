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
 *	\file src/Bounds.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0embedded
 *	\date 2007-2014
 *
 *	Implementation of the Bounds class designed to manage working sets of
 *	bounds within a QProblem.
 */


#include <qpOASES/Bounds.h>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	B o u n d s
 */
void BoundsCON( Bounds* THIS, int _n )
{
	Bounds_init( THIS,_n );
}


/*
 *	c o p y
 */
void BoundsCPY(	Bounds* FROM,
				Bounds* TO
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

	IndexlistCPY( &(FROM->freee),&(TO->freee) );
	IndexlistCPY( &(FROM->fixed),&(TO->fixed) );
}



/*
 *	i n i t
 */
returnValue Bounds_init(	Bounds* THIS,
							int _n
							)
{
	int i;
	
	if ( _n < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( _n > 0 )
	{
		Indexlist_init( &(THIS->freee),_n );
		Indexlist_init( &(THIS->fixed),_n );
	}

	THIS->n = _n;
	THIS->noLower = BT_TRUE;
	THIS->noUpper = BT_TRUE;

	assert( THIS->n <= NVMAX );

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
 *	s e t u p B o u n d
 */
returnValue Bounds_setupBound(	Bounds* THIS, int number, SubjectToStatus _status
								)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Add bound index to respective index list. */
	switch ( _status )
	{
		case ST_INACTIVE:
			if ( Bounds_addIndex( THIS,Bounds_getFree( THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		case ST_LOWER:
			if ( Bounds_addIndex( THIS,Bounds_getFixed( THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		case ST_UPPER:
			if ( Bounds_addIndex( THIS,Bounds_getFixed( THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		default:
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A l l F r e e
 */
returnValue Bounds_setupAllFree( Bounds* THIS )
{
	return Bounds_setupAll( THIS,ST_INACTIVE );
}


/*
 *	s e t u p A l l L o w e r
 */
returnValue Bounds_setupAllLower( Bounds* THIS )
{
	return Bounds_setupAll( THIS,ST_LOWER );
}


/*
 *	s e t u p A l l U p p e r
 */
returnValue Bounds_setupAllUpper( Bounds* THIS )
{
	return Bounds_setupAll( THIS,ST_UPPER );
}


/*
 *	m o v e F i x e d T o F r e e
 */
returnValue Bounds_moveFixedToFree( Bounds* THIS, int number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of fixed variables to that of free ones. */
	if ( Bounds_removeIndex( THIS,Bounds_getFixed( THIS ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( Bounds_addIndex( THIS,Bounds_getFree( THIS ),number,ST_INACTIVE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	m o v e F r e e T o F i x e d
 */
returnValue Bounds_moveFreeToFixed(	Bounds* THIS, int number, SubjectToStatus _status
									)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of free variables to that of fixed ones. */
	if ( Bounds_removeIndex( THIS,Bounds_getFree( THIS ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( Bounds_addIndex( THIS,Bounds_getFixed( THIS ),number,_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	f l i p F i x e d
 */
returnValue Bounds_flipFixed( Bounds* THIS, int number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( THIS->status != 0 )
		switch (THIS->status[number])
		{
			case ST_LOWER: THIS->status[number] = ST_UPPER; break;
			case ST_UPPER: THIS->status[number] = ST_LOWER; break;
			default: return THROWERROR( RET_MOVING_BOUND_FAILED );
		}

	return SUCCESSFUL_RETURN;
}


/*
 *	s w a p F r e e
 */
returnValue Bounds_swapFree(	Bounds* THIS, int number1, int number2
								)
{
	/* consistency check */
	if ( ( number1 < 0 ) || ( number1 >= THIS->n ) || ( number2 < 0 ) || ( number2 >= THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Swap index within indexlist of free variables. */
	return Bounds_swapIndex( THIS,Bounds_getFree( THIS ),number1,number2 );
}


/*
 *	s h i f t
 */
returnValue Bounds_shift(	Bounds* THIS, int offset )
{
	int i;

	static Indexlist shiftedFreee;
	static Indexlist shiftedFixed;

	Indexlist_init( &shiftedFreee,THIS->n );
	Indexlist_init( &shiftedFixed,THIS->n );

	/* consistency check */
	if ( ( offset == 0 ) || ( THIS->n <= 1 ) )
		return SUCCESSFUL_RETURN;

	if ( ( offset < 0 ) || ( offset > THIS->n/2 ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( ( THIS->n % offset ) != 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Shift types and THIS->status. */
	for( i=0; i<THIS->n-offset; ++i )
	{
		Bounds_setType( THIS,i,Bounds_getType( THIS,i+offset ) );
		Bounds_setStatus( THIS,i,Bounds_getStatus( THIS,i+offset ) );
	}

	/* 2) Construct shifted index lists of free and fixed variables. */
	for( i=0; i<THIS->n; ++i )
	{
		switch ( Bounds_getStatus( THIS,i ) )
		{
			case ST_INACTIVE:
				if ( Indexlist_addNumber( &shiftedFreee,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_LOWER:
				if ( Indexlist_addNumber( &shiftedFixed,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_UPPER:
				if ( Indexlist_addNumber( &shiftedFixed,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			default:
				return THROWERROR( RET_SHIFTING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	IndexlistCPY( &shiftedFreee,&(THIS->freee) );
	IndexlistCPY( &shiftedFixed,&(THIS->fixed) );

	return SUCCESSFUL_RETURN;
}


/*
 *	r o t a t e
 */
returnValue Bounds_rotate( Bounds* THIS, int offset )
{
	int i;
	myStatic SubjectToType   typeTmp[NVMAX];
	myStatic SubjectToStatus statusTmp[NVMAX];
	
	myStatic Indexlist rotatedFreee;
	myStatic Indexlist rotatedFixed;

	Indexlist_init( &rotatedFreee,THIS->n );
	Indexlist_init( &rotatedFixed,THIS->n );

	/* consistency check */
	if ( ( offset == 0 ) || ( offset == THIS->n ) || ( THIS->n <= 1 ) )
		return SUCCESSFUL_RETURN;

	if ( ( offset < 0 ) || ( offset > THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );


	/* 1) Rotate types and status. */
	for( i=0; i<offset; ++i )
	{
		typeTmp[i] = Bounds_getType( THIS,i );
		statusTmp[i] = Bounds_getStatus( THIS,i );
	}

	for( i=0; i<THIS->n-offset; ++i )
	{
		Bounds_setType( THIS,i,Bounds_getType( THIS,i+offset ) );
		Bounds_setStatus( THIS,i,Bounds_getStatus( THIS,i+offset ) );
	}

	for( i=THIS->n-offset; i<THIS->n; ++i )
	{
		Bounds_setType( THIS,i,typeTmp[i-THIS->n+offset] );
		Bounds_setStatus( THIS,i,statusTmp[i-THIS->n+offset] );
	}

	/* 2) Construct shifted index lists of free and fixed variables. */
	for( i=0; i<THIS->n; ++i )
	{
		switch ( Bounds_getStatus( THIS,i ) )
		{
			case ST_INACTIVE:
				if ( Indexlist_addNumber( &rotatedFreee,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_LOWER:
				if ( Indexlist_addNumber( &rotatedFixed,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_UPPER:
				if ( Indexlist_addNumber( &rotatedFixed,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			default:
				return THROWERROR( RET_ROTATING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	IndexlistCPY( &rotatedFreee,&(THIS->freee) );
	IndexlistCPY( &rotatedFixed,&(THIS->fixed) );

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue Bounds_print( Bounds* THIS )
{
	#ifndef __XPCTARGET__
	#ifndef __DSPACE__
	static char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	int nFR = Bounds_getNFR( THIS );
	int nFX = Bounds_getNFX( THIS );

	int *FR_idx, *FX_idx;

	if ( THIS->n == 0 )
		return SUCCESSFUL_RETURN;

	Indexlist_getNumberArray( Bounds_getFree( THIS ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( THIS ),&FX_idx );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"Bounds object comprising %d variables (%d free, %d fixed):\n",THIS->n,nFR,nFX );
	qpOASES_myPrintf( myPrintfString );

	REFER_NAMESPACE_QPOASES qpOASES_printNI( FR_idx,nFR,"free " );
	REFER_NAMESPACE_QPOASES qpOASES_printNI( FX_idx,nFX,"fixed" );

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
returnValue Bounds_setupAll( Bounds* THIS, SubjectToStatus _status )
{
	int i;

	/* 1) Place unbounded variables at the beginning of the index list of free variables. */
	for( i=0; i<THIS->n; ++i )
	{
		if ( Bounds_getType( THIS,i ) == ST_UNBOUNDED )
		{
			if ( Bounds_setupBound( THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 2) Add remaining (i.e. bounded but possibly free) variables to the index list of free variables. */
	for( i=0; i<THIS->n; ++i )
	{
		if ( Bounds_getType( THIS,i ) == ST_BOUNDED )
		{
			if ( Bounds_setupBound( THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 3) Place implicitly fixed variables at the end of the index list of free variables. */
	for( i=0; i<THIS->n; ++i )
	{
		if ( Bounds_getType( THIS,i ) == ST_EQUALITY )
		{
			if ( Bounds_setupBound( THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 4) Moreover, add all bounds of unknown type. */
	for( i=0; i<THIS->n; ++i )
	{
		if ( Bounds_getType( THIS,i ) == ST_UNKNOWN || Bounds_getType( THIS,i ) == ST_DISABLED )
		{
			if ( Bounds_setupBound( THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	a d d I n d e x
 */
returnValue Bounds_addIndex(	Bounds* THIS, Indexlist* const indexlist,
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
returnValue Bounds_removeIndex(	Bounds* THIS, Indexlist* const indexlist,
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
returnValue Bounds_swapIndex(	Bounds* THIS, Indexlist* const indexlist,
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
