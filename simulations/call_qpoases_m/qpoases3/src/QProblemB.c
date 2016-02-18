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
 *	\file src/QProblemB.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0embedded
 *	\date 2007-2014
 *
 *	Implementation of the QProblemB class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming.
 */


#include <qpOASES/QProblemB.h>


BEGIN_NAMESPACE_QPOASES


/*
 *	Q P r o b l e m B
 */
void QProblemBCON(	QProblemB* THIS,
					int _nV, HessianType _hessianType )
{
	#ifdef __CODE_GENERATION__
	Options_setToFast( &(THIS->options) );
	#else
	Options_setToDefault( &(THIS->options) );
	#endif /* __CODE_GENERATION__ */

	/* print copyright notice */
	if (THIS->options.printLevel != PL_NONE)
		qpOASES_printCopyrightNotice( );

	/* consistency check */
	if ( ( _nV <= 0 ) || ( _nV > NVMAX ) )
	{
		_nV = 1;
		THROWERROR( RET_INVALID_ARGUMENTS );
		assert( 1 == 0 );
	}

	/* reset global message handler */
	MessageHandling_reset( qpOASES_getGlobalMessageHandler() );

	THIS->H = &(THIS->HH);

	Bounds_init( &(THIS->bounds),_nV );

	THIS->tau = 0.0;

	THIS->hessianType = _hessianType;
	THIS->regVal = 0.0;

	THIS->infeasible  = BT_FALSE;
	THIS->unbounded   = BT_FALSE;

	THIS->status = QPS_NOTINITIALISED;

	THIS->count = 0;

	THIS->ramp0 = THIS->options.initialRamping;
	THIS->ramp1 = THIS->options.finalRamping;
	THIS->rampOffset = 0;

	QProblemB_setPrintLevel( THIS,THIS->options.printLevel );
}


/*
 *	c o p y
 */
void QProblemBCPY(	QProblemB* FROM,
					QProblemB* TO
					)
{
	unsigned int _nV = (unsigned int)QProblemB_getNV( FROM );

	TO->bounds = FROM->bounds;

	TO->HH = FROM->HH;
	TO->H = &(TO->HH);

	QProblemB_setG( TO,FROM->g );
	QProblemB_setLB( TO,FROM->lb );
	QProblemB_setUB( TO,FROM->ub );

	memcpy( TO->R,FROM->R,NVMAX*NVMAX*sizeof(real_t) );
	
	TO->haveCholesky = FROM->haveCholesky;

	memcpy( TO->x,FROM->x,_nV*sizeof(real_t) );
	memcpy( TO->y,FROM->y,_nV*sizeof(real_t) );

	TO->tau = FROM->tau;

	TO->hessianType = FROM->hessianType;
	TO->regVal = FROM->regVal;

	TO->infeasible = FROM->infeasible;
	TO->unbounded = FROM->unbounded;

	TO->status = FROM->status;

	TO->count = FROM->count;

	TO->ramp0 = FROM->ramp0;
	TO->ramp1 = FROM->ramp1;

	OptionsCPY( &(FROM->options),&(TO->options) );
	QProblemB_setPrintLevel( TO,TO->options.printLevel );
}



/*
 *	r e s e t
 */
returnValue QProblemB_reset( QProblemB* THIS )
{
	int i;
	int nV = QProblemB_getNV( THIS );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Reset bounds. */
	Bounds_init( &(THIS->bounds),nV );

	/* 2) Reset Cholesky decomposition. */
	for( i=0; i<NVMAX*NVMAX; ++i )
		THIS->R[i] = 0.0;
	
	THIS->haveCholesky = BT_FALSE;

	/* 3) Reset steplength and status flags. */
	THIS->tau = 0.0;

	THIS->hessianType = HST_UNKNOWN;
	THIS->regVal = 0.0;

	THIS->infeasible  = BT_FALSE;
	THIS->unbounded   = BT_FALSE;

	THIS->status = QPS_NOTINITIALISED;

	THIS->ramp0 = THIS->options.initialRamping;
	THIS->ramp1 = THIS->options.finalRamping;
	THIS->rampOffset = 0;

	return SUCCESSFUL_RETURN;
}


/*
 *	i n i t
 */
returnValue QProblemB_initM(	QProblemB* THIS, DenseMatrix *_H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int* nWSR, real_t* const cputime
								)
{
	if ( QProblemB_getNV( THIS ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( QProblemB_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblemB_reset( THIS );
	}

	/* 2) Setup QP data. */
	if ( QProblemB_setupQPdataM( THIS,_H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine (without any additional information). */
	return QProblemB_solveInitialQP( THIS,0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB_init(	QProblemB* THIS, real_t* const _H, const real_t* const _g,
							const real_t* const _lb, const real_t* const _ub,
							int* nWSR, real_t* const cputime
							)
{
	if ( QProblemB_getNV( THIS ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( QProblemB_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblemB_reset( THIS );
	}

	/* 2) Setup QP data. */
	if ( QProblemB_setupQPdata( THIS,_H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine (without any additional information). */
	return QProblemB_solveInitialQP( THIS,0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB_initF(	QProblemB* THIS, const char* const H_file, const char* const g_file,
								const char* const lb_file, const char* const ub_file,
								int* nWSR, real_t* const cputime
								)
{
	if ( QProblemB_getNV( THIS ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( QProblemB_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblemB_reset( THIS );
	}

	/* 2) Setup QP data from files. */
	if ( QProblemB_setupQPdataFromFile( THIS,H_file,g_file,lb_file,ub_file ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* 3) Call to main initialisation routine (without any additional information). */
	return QProblemB_solveInitialQP( THIS,0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB_initMW( 	QProblemB* THIS, DenseMatrix *_H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int* nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								Bounds* const guessedBounds
								)
{
	int i;
	int nV = QProblemB_getNV( THIS );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( QProblemB_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblemB_reset( THIS );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( Bounds_getStatus( guessedBounds,i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	/* exclude THIS possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data. */
	if ( QProblemB_setupQPdataM( THIS,_H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return QProblemB_solveInitialQP( THIS,xOpt,yOpt,guessedBounds, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB_initW( 	QProblemB* THIS, real_t* const _H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int* nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								Bounds* const guessedBounds
								)
{
	int i;
	int nV = QProblemB_getNV( THIS );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( QProblemB_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblemB_reset( THIS );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( Bounds_getStatus( guessedBounds,i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	/* exclude THIS possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data. */
	if ( QProblemB_setupQPdata( THIS,_H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return QProblemB_solveInitialQP( THIS,xOpt,yOpt,guessedBounds, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB_initFW( 	QProblemB* THIS, const char* const H_file, const char* const g_file,
								const char* const lb_file, const char* const ub_file,
								int* nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								Bounds* const guessedBounds
								)
{
	int i;
	int nV = QProblemB_getNV( THIS );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( QProblemB_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblemB_reset( THIS );
	}

	for( i=0; i<nV; ++i )
	{
		if ( Bounds_getStatus( guessedBounds,i ) == ST_UNDEFINED )
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* exclude THIS possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data from files. */
	if ( QProblemB_setupQPdataFromFile( THIS,H_file,g_file,lb_file,ub_file ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* 3) Call to main initialisation routine. */
	return QProblemB_solveInitialQP( THIS,xOpt,yOpt,guessedBounds, nWSR,cputime );
}


/*
 * s e t u p I n i t i a l C h o l e s k y
 */
returnValue QProblemB_setupInitialCholesky( QProblemB* THIS )
{
	returnValue returnvalueCholesky;

	/* If regularisation shall be used, always regularise at beginning 
	 * if initial working set is not empty. */
	if ( ( QProblemB_getNV( THIS ) != QProblemB_getNFR( THIS ) - QProblemB_getNFV( THIS ) ) && ( THIS->options.enableRegularisation == BT_TRUE ) )
		if ( QProblemB_regulariseHessian( THIS ) != SUCCESSFUL_RETURN )
			return RET_INIT_FAILED_REGULARISATION;

	/* Factorise projected Hessian 
	 * now handles all special cases (no active bounds/constraints, no nullspace) */
	returnvalueCholesky = QProblemB_computeCholesky( THIS );

	/* If Hessian is not positive definite, regularise and try again. */
	if ( returnvalueCholesky == RET_HESSIAN_NOT_SPD )
	{
		if ( QProblemB_regulariseHessian( THIS ) != SUCCESSFUL_RETURN )
			return RET_INIT_FAILED_REGULARISATION;

		returnvalueCholesky = QProblemB_computeCholesky( THIS );
	}

	if ( returnvalueCholesky != SUCCESSFUL_RETURN )
		return RET_INIT_FAILED_CHOLESKY;

	THIS->haveCholesky = BT_TRUE;
	return SUCCESSFUL_RETURN;
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB_hotstart(	QProblemB* THIS, const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								int* nWSR, real_t* const cputime
								)
{
	returnValue returnvalue = SUCCESSFUL_RETURN;
	int i, nActiveFar;
	int nV = QProblemB_getNV( THIS );

	int nWSR_max = *nWSR;
	int nWSR_performed = 0;

	real_t cputime_remaining = QPOASES_INFTY;
	real_t cputime_needed = 0.0;

	real_t farbound = THIS->options.initialFarBounds;

	myStatic real_t ub_new_far[NVMAX];
	myStatic real_t lb_new_far[NVMAX];
	
	real_t tol;

	if ( QProblemB_getNV( THIS ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* Simple check for consistency of bounds */
	if ( QProblemB_areBoundsConsistent(THIS,lb_new, ub_new) != SUCCESSFUL_RETURN )
		return QProblemB_setInfeasibilityFlag(THIS,returnvalue,BT_TRUE);

	++(THIS->count);


	if ( THIS->haveCholesky == BT_FALSE )
	{
		returnvalue = QProblemB_setupInitialCholesky( THIS );
		if (returnvalue != SUCCESSFUL_RETURN)
			return THROWERROR(returnvalue);
	}

	if ( THIS->options.enableFarBounds == BT_FALSE )
	{
		/* Automatically call standard solveQP if regularisation is not active. */
		returnvalue = QProblemB_solveRegularisedQP( THIS,g_new,lb_new,ub_new, nWSR,cputime,0 );
	}
	else
	{
		/* possibly extend initial far bounds to largest bound/constraint data */
		if (ub_new)
			for (i = 0; i < nV; i++)
				if ((ub_new[i] < QPOASES_INFTY) && (ub_new[i] > farbound)) farbound = ub_new[i];
		if (lb_new)
			for (i = 0; i < nV; i++)
				if ((lb_new[i] > -QPOASES_INFTY) && (lb_new[i] < -farbound)) farbound = -lb_new[i];

		QProblemB_updateFarBounds(	THIS,farbound,nV,
									lb_new,lb_new_far, ub_new,ub_new_far
									);

		for ( ;; )
		{
			*nWSR = nWSR_max;
			if ( cputime != 0 )
				cputime_remaining = *cputime - cputime_needed;

			/* Automatically call standard solveQP if regularisation is not active. */
			returnvalue = QProblemB_solveRegularisedQP( THIS,g_new,lb_new_far,ub_new_far, nWSR,&cputime_remaining,nWSR_performed );
			
			nWSR_performed  = *nWSR;
			cputime_needed += cputime_remaining;

			/* Check for active far-bounds and move them away */
			nActiveFar = 0;
			farbound *= THIS->options.growFarBounds;

			if ( THIS->infeasible == BT_TRUE )
			{
				if ( farbound >= QPOASES_INFTY )
				{
					returnvalue = RET_HOTSTART_STOPPED_INFEASIBILITY;
					goto farewell;
				}

				QProblemB_updateFarBounds(	THIS,farbound,nV,
											lb_new,lb_new_far, ub_new,ub_new_far
											);
			}
			else if ( THIS->status == QPS_SOLVED )
			{
				tol = farbound/THIS->options.growFarBounds * THIS->options.boundTolerance;
				THIS->status = QPS_HOMOTOPYQPSOLVED;

				for ( i=0; i<nV; ++i )
				{
					if ( ( ( lb_new == 0 ) || ( lb_new_far[i] > lb_new[i] ) ) && ( qpOASES_getAbs ( lb_new_far[i] - THIS->x[i] ) < tol ) )
						++nActiveFar;
					if ( ( ( ub_new == 0 ) || ( ub_new_far[i] < ub_new[i] ) ) && ( qpOASES_getAbs ( ub_new_far[i] - THIS->x[i] ) < tol ) )
						++nActiveFar;
				}

				if ( nActiveFar == 0 )
					break;

				if ( farbound >= QPOASES_INFTY )
				{
					THIS->unbounded = BT_TRUE;
					returnvalue = RET_HOTSTART_STOPPED_UNBOUNDEDNESS;
					goto farewell;
				}

				QProblemB_updateFarBounds(	THIS,farbound,nV,
											lb_new,lb_new_far, ub_new,ub_new_far
											);
			}
			else
			{
				/* some other error when solving QP */
				break;
			}

			/* advance ramp offset to avoid Ramping cycles */
			(THIS->rampOffset)++;
		}

		farewell:
			if ( cputime != 0 )
				*cputime = cputime_needed;
	}

	return ( returnvalue != SUCCESSFUL_RETURN ) ? THROWERROR( returnvalue ) : returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB_hotstartF(	QProblemB* THIS, const char* const g_file,
									const char* const lb_file, const char* const ub_file,
									int* nWSR, real_t* const cputime
									)
{
	int nV  = QProblemB_getNV( THIS );
	returnValue returnvalue;

	/* 1) Allocate memory (if bounds exist). */
	myStatic real_t g_new[NVMAX];
	myStatic real_t lb_new[NVMAX];
	myStatic real_t ub_new[NVMAX];


	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* consistency check */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 2) Load new QP vectors from file. */
	returnvalue = QProblemB_loadQPvectorsFromFile(	THIS,g_file,lb_file,ub_file,
													g_new,lb_new,ub_new
													);
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* 3) Actually perform hotstart. */
	returnvalue = QProblemB_hotstart( THIS,g_new,lb_new,ub_new, nWSR,cputime );

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB_hotstartW(	QProblemB* THIS, const real_t* const g_new,
									const real_t* const lb_new, const real_t* const ub_new,
									int* nWSR, real_t* const cputime,
									Bounds* const guessedBounds
									)
{
	int nV = QProblemB_getNV( THIS );
	
	returnValue returnvalue;
	real_t starttime = 0.0;
	
	myStatic Bounds emptyBounds;
	BoundsCON( &emptyBounds,nV );


	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );


	/* start runtime measurement */
	if ( cputime != 0 )
		starttime = qpOASES_getCPUtime( );


	/* 1) Update working set according to guess for working set of bounds. */
	if ( guessedBounds != 0 )
	{
		if ( QProblemB_setupAuxiliaryQP( THIS,guessedBounds ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}
	else
	{
		/* create empty bounds for setting up auxiliary QP */
		if ( Bounds_setupAllFree( &emptyBounds ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( QProblemB_setupAuxiliaryQP( THIS,&emptyBounds ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}

	/* 2) Perform usual homotopy. */

	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= qpOASES_getCPUtime( ) - starttime;

	returnvalue = QProblemB_hotstart( THIS,g_new,lb_new,ub_new, nWSR,cputime );

	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = qpOASES_getCPUtime( ) - starttime;

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB_hotstartFW(	QProblemB* THIS, const char* const g_file,
									const char* const lb_file, const char* const ub_file,
									int* nWSR, real_t* const cputime,
									Bounds* const guessedBounds
									)
{
	int nV = QProblemB_getNV( THIS );
	returnValue returnvalue;

	/* 1) Allocate memory (if bounds exist). */
	myStatic real_t g_new[NVMAX];
	myStatic real_t lb_new[NVMAX];
	myStatic real_t ub_new[NVMAX];

	
	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* consistency check */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 2) Load new QP vectors from file. */
	returnvalue = QProblemB_loadQPvectorsFromFile(	THIS,g_file,lb_file,ub_file,
													g_new,lb_new,ub_new
													);
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* 3) Actually perform hotstart using initialised homotopy. */
	returnvalue = QProblemB_hotstartW(	THIS,g_new,lb_new,ub_new, nWSR,cputime,
										guessedBounds
										);

	return returnvalue;
}



/*
 *	g e t W o r k i n g S e t
 */
returnValue QProblemB_getWorkingSet( QProblemB* THIS, real_t* workingSet )
{
	int i;
	int nV = QProblemB_getNV( THIS );

	/* At which limit is the bound active? */
	for (i = 0; i < nV; i++) {
		switch ( Bounds_getStatus( &(THIS->bounds),i ) ) {
			case ST_LOWER: workingSet[i] = -1.0; break;
			case ST_UPPER: workingSet[i] = +1.0; break;
			default: workingSet[i] = 0.0; break;
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	g e t N Z
 */
int QProblemB_getNZ( QProblemB* THIS )
{
	/* if no constraints are present: nZ=nFR */
	return QProblemB_getNFR( THIS );
}


/*
 *	g e t O b j V a l
 */
real_t QProblemB_getObjVal( QProblemB* THIS )
{
	real_t objVal;

	/* calculated optimal objective function value
	 * only if current QP has been solved */
	if ( ( QProblemB_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblemB_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblemB_getStatus( THIS ) == QPS_SOLVED ) )
	{
		objVal = QProblemB_getObjValX( THIS,THIS->x );
	}
	else
	{
		objVal = QPOASES_INFTY;
	}

	return objVal;
}


/*
 *	g e t O b j V a l
 */
real_t QProblemB_getObjValX( QProblemB* THIS, const real_t* const _x )
{
	int i;
	int nV = QProblemB_getNV( THIS );

	real_t objVal = 0.0;
	myStatic real_t Hx[NVMAX];

	if ( nV == 0 )
		return 0.0;

	for( i=0; i<nV; ++i )
		objVal += _x[i]*THIS->g[i];

	switch ( THIS->hessianType )
	{
		case HST_ZERO:
			break;

		case HST_IDENTITY:
			for( i=0; i<nV; ++i )
				objVal += 0.5*_x[i]*_x[i];
			break;

		default:
			DenseMatrix_times(THIS->H,1, 1.0, _x, nV, 0.0, Hx, nV);
			for( i=0; i<nV; ++i )
				objVal += 0.5*_x[i]*Hx[i];
			break;
	}

	/* When using regularisation, the objective function value
	 * needs to be modified as follows:
	 * objVal = objVal - 0.5*_x*(Hmod-H)*_x - _x'*(gMod-g)
	 *        = objVal - 0.5*_x*eps*_x * - _x'*(-eps*_x)
	 *        = objVal + 0.5*_x*eps*_x */
	if ( QProblemB_usingRegularisation( THIS ) == BT_TRUE )
	{
		for( i=0; i<nV; ++i )
			objVal += 0.5*_x[i]*THIS->regVal*_x[i];
	}

	return objVal;
}


/*
 *	g e t P r i m a l S o l u t i o n
 */
returnValue QProblemB_getPrimalSolution( QProblemB* THIS, real_t* const xOpt )
{
	int i;

	/* return optimal primal solution vector
	 * only if current QP has been solved */
	if ( ( QProblemB_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblemB_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblemB_getStatus( THIS ) == QPS_SOLVED ) )
	{
		for( i=0; i<QProblemB_getNV( THIS ); ++i )
			xOpt[i] = THIS->x[i];

		return SUCCESSFUL_RETURN;
	}
	else
	{
		return RET_QP_NOT_SOLVED;
	}
}


/*
 *	g e t D u a l S o l u t i o n
 */
returnValue QProblemB_getDualSolution( QProblemB* THIS, real_t* const yOpt )
{
	int i;

	for( i=0; i<QProblemB_getNV( THIS ); ++i )
		yOpt[i] = THIS->y[i];

	/* return optimal dual solution vector
	 * only if current QP has been solved */
	if ( ( QProblemB_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblemB_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblemB_getStatus( THIS ) == QPS_SOLVED ) )
	{
		return SUCCESSFUL_RETURN;
	}
	else
	{
		return RET_QP_NOT_SOLVED;
	}
}



/*
 *	p r i n t P r o p e r t i e s
 */
returnValue QProblemB_printProperties( QProblemB* THIS )
{
	#ifndef __XPCTARGET__
	#ifndef __DSPACE__
	myStatic char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	/* Do not print properties if print level is set to none! */
	if ( THIS->options.printLevel == PL_NONE )
		return SUCCESSFUL_RETURN;

	qpOASES_myPrintf( "\n#################   qpOASES  --  QP PROPERTIES   #################\n" );
	qpOASES_myPrintf( "\n" );

	/* 1) Variables properties. */
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,  "Number of Variables: %4.1d\n",QProblemB_getNV( THIS ) );
	qpOASES_myPrintf( myPrintfString );

	if ( Bounds_hasNoLower( &(THIS->bounds) ) == BT_TRUE )
			qpOASES_myPrintf( "Variables are not bounded from below.\n" );
		else
			qpOASES_myPrintf( "Variables are bounded from below.\n" );

	if ( Bounds_hasNoUpper( &(THIS->bounds) ) == BT_TRUE )
			qpOASES_myPrintf( "Variables are not bounded from above.\n" );
		else
			qpOASES_myPrintf( "Variables are bounded from above.\n" );

	qpOASES_myPrintf( "\n" );


	/* 2) Further properties. */
	switch ( THIS->hessianType )
	{
		case HST_ZERO:
			qpOASES_myPrintf( "Hessian is zero matrix (i.e. actually an LP is solved).\n" );
			break;

		case HST_IDENTITY:
			qpOASES_myPrintf( "Hessian is identity matrix.\n" );
			break;

		case HST_POSDEF:
			qpOASES_myPrintf( "Hessian matrix is (strictly) positive definite.\n" );
			break;

		case HST_POSDEF_NULLSPACE:
			qpOASES_myPrintf( "Hessian matrix is positive definite on null space of active constraints.\n" );
			break;

		case HST_SEMIDEF:
			qpOASES_myPrintf( "Hessian matrix is positive semi-definite.\n" );
			break;

		case HST_INDEF:
			qpOASES_myPrintf( "Hessian matrix is indefinite.\n" );
			break;

		default:
			qpOASES_myPrintf( "Hessian matrix has unknown type.\n" );
			break;
	}

	if ( THIS->infeasible == BT_TRUE )
		qpOASES_myPrintf( "QP was found to be infeasible.\n" );
	else
		qpOASES_myPrintf( "QP seems to be feasible.\n" );

	if ( THIS->unbounded == BT_TRUE )
		qpOASES_myPrintf( "QP was found to be unbounded from below.\n" );
	else
		qpOASES_myPrintf( "QP seems to be bounded from below.\n" );

	qpOASES_myPrintf( "\n" );


	/* 3) QP object properties. */
	switch ( THIS->status )
	{
		case QPS_NOTINITIALISED:
			qpOASES_myPrintf( "Status of QP object: freshly instantiated or reset.\n" );
			break;

		case QPS_PREPARINGAUXILIARYQP:
			qpOASES_myPrintf( "Status of QP object: an auxiliary QP is currently setup.\n" );
			break;

		case QPS_AUXILIARYQPSOLVED:
			qpOASES_myPrintf( "Status of QP object: an auxilary QP was solved.\n" );
			break;

		case QPS_PERFORMINGHOMOTOPY:
			qpOASES_myPrintf( "Status of QP object: a homotopy step is performed.\n" );
			break;

		case QPS_HOMOTOPYQPSOLVED:
			qpOASES_myPrintf( "Status of QP object: an intermediate QP along the homotopy path was solved.\n" );
			break;

		case QPS_SOLVED:
			qpOASES_myPrintf( "Status of QP object: solution of the actual QP was found.\n" );
			break;
	}

	switch ( THIS->options.printLevel )
	{
		case PL_DEBUG_ITER:
			qpOASES_myPrintf( "Print level of QP object is set to display a tabular output for debugging.\n" );
			break;

		case PL_TABULAR:
			qpOASES_myPrintf( "Print level of QP object is set to display a tabular output.\n" );
			break;

		case PL_LOW:
					qpOASES_myPrintf( "Print level of QP object is low, i.e. only error are printed.\n" );
			break;

		case PL_MEDIUM:
			qpOASES_myPrintf( "Print level of QP object is medium, i.e. error and warnings are printed.\n" );
			break;

		case PL_HIGH:
			qpOASES_myPrintf( "Print level of QP object is high, i.e. all available output is printed.\n" );
			break;

		default:
			break;
	}

	qpOASES_myPrintf( "\n" );
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}


returnValue QProblemB_printOptions( QProblemB* THIS )
{
	return Options_print( &(THIS->options) );
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/




/*
 *	d e t e r m i n e H e s s i a n T y p e
 */
returnValue QProblemB_determineHessianType( QProblemB* THIS )
{
	int i;
	real_t curDiag;
	BooleanType isIdentity, isZero;

	int nV = QProblemB_getNV( THIS );
	

	/* if Hessian type has been set by user, do NOT change it! */
	if ( THIS->hessianType != HST_UNKNOWN )
		return SUCCESSFUL_RETURN;

	/* if Hessian has not been allocated, assume it to be all zeros! */
	if ( THIS->H == 0 )
	{
		THIS->hessianType = HST_ZERO;

		if ( THIS->options.enableRegularisation == BT_FALSE )
			return THROWERROR( RET_USE_REGULARISATION_FOR_LP );

		return SUCCESSFUL_RETURN;
	}

	/* 1) If Hessian has outer-diagonal elements,
	 *    Hessian is assumed to be positive definite. */
	THIS->hessianType = HST_POSDEF;
	if (DenseMatrix_isDiag(THIS->H) == BT_FALSE)
		return SUCCESSFUL_RETURN;

	/* 2) Otherwise it is diagonal and test for identity or zero matrix is performed. */
	isIdentity = BT_TRUE;
	isZero = BT_TRUE;

	for ( i=0; i<nV; ++i )
	{
		curDiag = DenseMatrix_diag( THIS->H,i );
        if ( curDiag >= QPOASES_INFTY )
            return RET_DIAGONAL_NOT_INITIALISED;

		if ( curDiag < -QPOASES_ZERO )
		{
			THIS->hessianType = HST_INDEF;
			if ( THIS->options.enableFlippingBounds == BT_FALSE )
				return THROWERROR( RET_HESSIAN_INDEFINITE );
			else
				return SUCCESSFUL_RETURN;
		}

		if ( qpOASES_getAbs( curDiag - 1.0 ) > QPOASES_EPS )
			isIdentity = BT_FALSE;

		if ( qpOASES_getAbs( curDiag ) > QPOASES_EPS )
			isZero = BT_FALSE;
	}

	if ( isIdentity == BT_TRUE )
		THIS->hessianType = HST_IDENTITY;

	if ( isZero == BT_TRUE )
		THIS->hessianType = HST_ZERO;

	if ( ( THIS->hessianType == HST_ZERO ) && ( THIS->options.enableRegularisation == BT_FALSE ) )
		return THROWERROR( RET_USE_REGULARISATION_FOR_LP );

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblemB_setupSubjectToType( QProblemB* THIS )
{
	return QProblemB_setupSubjectToTypeNew( THIS,THIS->lb,THIS->ub );
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblemB_setupSubjectToTypeNew( QProblemB* THIS, const real_t* const lb_new, const real_t* const ub_new )
{
	int i;
	int nV = QProblemB_getNV( THIS );


	/* 1) Check if lower bounds are present. */
	Bounds_setNoLower( &(THIS->bounds),BT_TRUE );
	if ( lb_new != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( lb_new[i] > -QPOASES_INFTY )
			{
				Bounds_setNoLower( &(THIS->bounds),BT_FALSE );
				break;
			}
		}
	}

	/* 2) Check if upper bounds are present. */
	Bounds_setNoUpper( &(THIS->bounds),BT_TRUE );
	if ( ub_new != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( ub_new[i] < QPOASES_INFTY )
			{
				Bounds_setNoUpper( &(THIS->bounds),BT_FALSE );
				break;
			}
		}
	}

	/* 3) Determine implicitly fixed and unbounded variables. */
	if ( ( lb_new != 0 ) && ( ub_new != 0 ) )
	{
		for( i=0; i<nV; ++i )
		{
			if ( ( lb_new[i] <= -QPOASES_INFTY ) && ( ub_new[i] >= QPOASES_INFTY )
					&& (THIS->options.enableFarBounds == BT_FALSE))
			{
				Bounds_setType( &(THIS->bounds),i,ST_UNBOUNDED );
			}
			else
			{
				if ( THIS->options.enableEqualities
						&& THIS->lb[i] > THIS->ub[i] - THIS->options.boundTolerance
						&& lb_new[i] > ub_new[i] - THIS->options.boundTolerance)
					Bounds_setType( &(THIS->bounds),i,ST_EQUALITY );
				else
					Bounds_setType( &(THIS->bounds),i,ST_BOUNDED );
			}
		}
	}
	else
	{
		if ( ( lb_new == 0 ) && ( ub_new == 0 ) )
		{
			for( i=0; i<nV; ++i )
				Bounds_setType( &(THIS->bounds),i,ST_UNBOUNDED );
		}
		else
		{
			for( i=0; i<nV; ++i )
				Bounds_setType( &(THIS->bounds),i,ST_BOUNDED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c o m p u t e C h o l e s k y
 */
returnValue QProblemB_computeCholesky( QProblemB* THIS )
{
	int i, j;
	int nV  = QProblemB_getNV( THIS );
	int nFR = QProblemB_getNFR( THIS );
	
	int* FR_idx;

	long info = 0;
	unsigned long _nFR = (unsigned long)nFR, _nV = NVMAX;
	
	/* 1) Initialises R with all zeros. */
	for( i=0; i<NVMAX*NVMAX; ++i )
		THIS->R[i] = 0.0;

	/* 2) Calculate Cholesky decomposition of H (projected to free variables). */
	if ( ( THIS->hessianType == HST_ZERO ) || ( THIS->hessianType == HST_IDENTITY ) )
	{
		if ( THIS->hessianType == HST_ZERO )
		{
			/* if Hessian is zero matrix and it has been regularised,
			 * its Cholesky factor is the identity matrix scaled by sqrt(eps). */
			if ( QProblemB_usingRegularisation( THIS ) == BT_TRUE )
			{
				for( i=0; i<nV; ++i )
					RR(i,i) = qpOASES_getSqrt( THIS->regVal );
			}
			/*else
				return THROWERROR( RET_CHOLESKY_OF_ZERO_HESSIAN );*/
		}
		else
		{
			/* if Hessian is identity, so is its Cholesky factor. */
			for( i=0; i<nV; ++i )
				RR(i,i) = 1.0;
		}
	}
	else
	{
		if ( nFR > 0 )
		{
			Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );

			/* get H */
			for ( j=0; j < nFR; ++j )
				DenseMatrix_getCol(THIS->H,FR_idx[j], Bounds_getFree( &(THIS->bounds) ), 1.0, &(THIS->R[j*NVMAX]));

			/* R'*R = H */
			POTRF( "U", &_nFR, THIS->R, &_nV, &info );

			/* <0 = invalid call, =0 ok, >0 not spd */
			if (info > 0) {
				if ( THIS->R[0] < 0.0 )
				{
					/* Cholesky decomposition has tunneled a negative
					 * diagonal element. */ 
					THIS->options.epsRegularisation = qpOASES_getMin( -THIS->R[0]+THIS->options.epsRegularisation,qpOASES_getSqrt(qpOASES_getAbs(THIS->options.epsRegularisation)) );
				}

				THIS->hessianType = HST_SEMIDEF;
				return RET_HESSIAN_NOT_SPD;
			}


			/* zero first subdiagonal to make givens updates work */
			for (i=0;i<nFR-1;++i)
				RR(i+1,i) = 0.0;

		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	o b t a i n A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblemB_obtainAuxiliaryWorkingSet(	QProblemB* THIS, const real_t* const xOpt, const real_t* const yOpt,
													Bounds* const guessedBounds, Bounds* auxiliaryBounds
													)
{
	int i = 0;
	int nV = QProblemB_getNV( THIS );


	/* 1) Ensure that desiredBounds is allocated (and different from guessedBounds). */
	if ( ( auxiliaryBounds == 0 ) || ( auxiliaryBounds == guessedBounds ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 2) Setup working set for auxiliary initial QP. */
	if ( guessedBounds != 0 )
	{
		/* If an initial working set is specific, use it!
		 * Moreover, add all implictly fixed variables if specified. */
		for( i=0; i<nV; ++i )
		{
			#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
			if ( Bounds_getType( &(THIS->bounds),i ) == ST_EQUALITY )
			{
				if ( Bounds_setupBound( auxiliaryBounds,i,ST_LOWER ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
			else
			#endif
			{
				if ( Bounds_setupBound( auxiliaryBounds,i,Bounds_getStatus( guessedBounds,i ) ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
		}
	}
	else	/* No initial working set specified. */
	{
		if ( ( xOpt != 0 ) && ( yOpt == 0 ) )
		{
			/* Obtain initial working set by "clipping". */
			for( i=0; i<nV; ++i )
			{
				if ( xOpt[i] <= THIS->lb[i] + THIS->options.boundTolerance )
				{
					if ( Bounds_setupBound( auxiliaryBounds,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( xOpt[i] >= THIS->ub[i] - THIS->options.boundTolerance )
				{
					if ( Bounds_setupBound( auxiliaryBounds,i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all implictly fixed variables if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( Bounds_getType( &(THIS->bounds),i ) == ST_EQUALITY )
				{
					if ( Bounds_setupBound( auxiliaryBounds,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( Bounds_setupBound( auxiliaryBounds,i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		if ( ( xOpt == 0 ) && ( yOpt != 0 ) )
		{
			/* Obtain initial working set in accordance to sign of dual solution vector. */
			for( i=0; i<nV; ++i )
			{
				if ( yOpt[i] > QPOASES_EPS )
				{
					if ( Bounds_setupBound( auxiliaryBounds,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( yOpt[i] < -QPOASES_EPS )
				{
					if ( Bounds_setupBound( auxiliaryBounds,i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all implictly fixed variables if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( Bounds_getType( &(THIS->bounds),i ) == ST_EQUALITY )
				{
					if ( Bounds_setupBound( auxiliaryBounds,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( Bounds_setupBound( auxiliaryBounds,i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		/* If xOpt and yOpt are null pointer and no initial working is specified,
		 * start with empty working set (or implicitly fixed bounds only)
		 * for auxiliary QP. */
		if ( ( xOpt == 0 ) && ( yOpt == 0 ) )
		{
			for( i=0; i<nV; ++i )
			{
				switch( Bounds_getType( &(THIS->bounds),i ) )
				{
					case ST_UNBOUNDED:
						if ( Bounds_setupBound( auxiliaryBounds,i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;

					/* Only add all implictly fixed variables if specified. */
					#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
					case ST_EQUALITY:
						if ( Bounds_setupBound( auxiliaryBounds,i,ST_LOWER ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;
					#endif

					default:
						if ( Bounds_setupBound( auxiliaryBounds,i,THIS->options.initialStatusBounds ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	b a c k s o l v e R
 */
returnValue QProblemB_backsolveR(	QProblemB* THIS, const real_t* const b, BooleanType transposed,
									real_t* const a
									)
{
	/* Call standard backsolve procedure (i.e. removingBound == BT_FALSE). */
	return QProblemB_backsolveRrem( THIS,b,transposed,BT_FALSE,a );
}


/*
 *	b a c k s o l v e R
 */
returnValue QProblemB_backsolveRrem(	QProblemB* THIS, const real_t* const b, BooleanType transposed,
										BooleanType removingBound,
										real_t* const a
										)
{
	int i, j;
	int nR = QProblemB_getNZ( THIS );

	real_t sum;

	/* if backsolve is called while removing a bound, reduce nZ by one. */
	if ( removingBound == BT_TRUE )
		--nR;

	/* nothing to do */
	if ( nR <= 0 )
		return SUCCESSFUL_RETURN;


	/* Solve Ra = b, where R might be transposed. */
	if ( transposed == BT_FALSE )
	{
		/* solve Ra = b */
		for( i=(nR-1); i>=0; --i )
		{
			sum = b[i];
			for( j=(i+1); j<nR; ++j )
				sum -= RR(i,j) * a[j];

			if ( qpOASES_getAbs( RR(i,i) ) >= QPOASES_ZERO*qpOASES_getAbs( sum ) )
				a[i] = sum / RR(i,i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}
	else
	{
		/* solve R^T*a = b */
		for( i=0; i<nR; ++i )
		{
			sum = b[i];
			for( j=0; j<i; ++j )
				sum -= RR(j,i) * a[j];

			if ( qpOASES_getAbs( RR(i,i) ) >= QPOASES_ZERO*qpOASES_getAbs( sum ) )
				a[i] = sum / RR(i,i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e D a t a S h i f t
 */
returnValue QProblemB_determineDataShift(	QProblemB* THIS, const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new,
											real_t* const delta_g, real_t* const delta_lb, real_t* const delta_ub,
											BooleanType* Delta_bB_isZero
											)
{
	int i, ii;
	int nV  = QProblemB_getNV( THIS );
	int nFX = QProblemB_getNFX( THIS );

	int* FX_idx;
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );


	/* 1) Calculate shift directions. */
	for( i=0; i<nV; ++i )
		delta_g[i]  = g_new[i]  - THIS->g[i];

	if ( lb_new != 0 )
	{
		for( i=0; i<nV; ++i )
			delta_lb[i] = lb_new[i] - THIS->lb[i];
	}
	else
	{
		/* if no lower bounds exist, assume the new lower bounds to be -infinity */
		for( i=0; i<nV; ++i )
			delta_lb[i] = -QPOASES_INFTY - THIS->lb[i];
	}

	if ( ub_new != 0 )
	{
		for( i=0; i<nV; ++i )
			delta_ub[i] = ub_new[i] - THIS->ub[i];
	}
	else
	{
		/* if no upper bounds exist, assume the new upper bounds to be infinity */
		for( i=0; i<nV; ++i )
			delta_ub[i] = QPOASES_INFTY - THIS->ub[i];
	}

	/* 2) Determine if active bounds are to be shifted. */
	*Delta_bB_isZero = BT_TRUE;

	for ( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];

		if ( ( qpOASES_getAbs( delta_lb[ii] ) > QPOASES_EPS ) || ( qpOASES_getAbs( delta_ub[ii] ) > QPOASES_EPS ) )
		{
			*Delta_bB_isZero = BT_FALSE;
			break;
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p Q P d a t a
 */
returnValue QProblemB_setupQPdataM(	QProblemB* THIS, DenseMatrix *_H, const real_t* const _g,
									const real_t* const _lb, const real_t* const _ub
									)
{
	if ( _H == 0 )
		return QProblemB_setupQPdata( THIS,(real_t*)0,_g,_lb,_ub );
	else
		return QProblemB_setupQPdata( THIS,DenseMatrix_getVal(_H),_g,_lb,_ub );
}


/*
 *	s e t u p Q P d a t a
 */
returnValue QProblemB_setupQPdata(	QProblemB* THIS, real_t* const _H, const real_t* const _g,
									const real_t* const _lb, const real_t* const _ub
									)
{
	int i;
	int nV = QProblemB_getNV( THIS );

	/* 1) Setup Hessian matrix. */
	QProblemB_setH( THIS,_H );

	/* 2) Setup gradient vector. */
	if ( _g == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );
	else
		QProblemB_setG( THIS,_g );

	/* 3) Setup lower bounds vector. */
	if ( _lb != 0 )
	{
		QProblemB_setLB( THIS,_lb );
	}
	else
	{
		/* if no lower bounds are specified, set them to -infinity */
		for( i=0; i<nV; ++i )
			THIS->lb[i] = -QPOASES_INFTY;
	}

	/* 4) Setup upper bounds vector. */
	if ( _ub != 0 )
	{
		QProblemB_setUB( THIS,_ub );
	}
	else
	{
		/* if no upper bounds are specified, set them to infinity */
		for( i=0; i<nV; ++i )
			THIS->ub[i] = QPOASES_INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p Q P d a t a F r o m F i l e
 */
returnValue QProblemB_setupQPdataFromFile(	QProblemB* THIS, const char* const H_file, const char* const g_file,
											const char* const lb_file, const char* const ub_file
											)
{
	int i;
	int nV = QProblemB_getNV( THIS );

	returnValue returnvalue;


	/* 1) Load Hessian matrix from file. */
	myStatic real_t _H[NVMAX*NVMAX];

	if ( H_file != 0 )
	{
		returnvalue = qpOASES_readFromFileM( _H, nV,nV, H_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
		QProblemB_setH( THIS,_H );
	}
	else
	{
		QProblemB_setH( THIS,(real_t*)0 );
	}

	/* 2) Load gradient vector from file. */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	returnvalue = qpOASES_readFromFileV( THIS->g, nV, g_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
		return THROWERROR( returnvalue );

	/* 3) Load lower bounds vector from file. */
	if ( lb_file != 0 )
	{
		returnvalue = qpOASES_readFromFileV( THIS->lb, nV, lb_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* if no lower bounds are specified, set them to -infinity */
		for( i=0; i<nV; ++i )
			THIS->lb[i] = -QPOASES_INFTY;
	}

	/* 4) Load upper bounds vector from file. */
	if ( ub_file != 0 )
	{
		returnvalue = qpOASES_readFromFileV( THIS->ub, nV, ub_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* if no upper bounds are specified, set them to infinity */
		for( i=0; i<nV; ++i )
			THIS->ub[i] = QPOASES_INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	l o a d Q P v e c t o r s F r o m F i l e
 */
returnValue QProblemB_loadQPvectorsFromFile(	QProblemB* THIS, const char* const g_file, const char* const lb_file, const char* const ub_file,
												real_t* const g_new, real_t* const lb_new, real_t* const ub_new
												)
{
	int nV = QProblemB_getNV( THIS );

	returnValue returnvalue;


	/* 1) Load gradient vector from file. */
	if ( ( g_file != 0 ) && ( g_new != 0 ) )
	{
		returnvalue = qpOASES_readFromFileV( g_new, nV, g_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* At least gradient vector needs to be specified! */
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* 2) Load lower bounds vector from file. */
	if ( lb_file != 0 )
	{
		if ( lb_new != 0 )
		{
			returnvalue = qpOASES_readFromFileV( lb_new, nV, lb_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* If filename is given, storage must be provided! */
			return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	/* 3) Load upper bounds vector from file. */
	if ( ub_file != 0 )
	{
		if ( ub_new != 0 )
		{
			returnvalue = qpOASES_readFromFileV( ub_new, nV, ub_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* If filename is given, storage must be provided! */
			return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t I n f e a s i b i l i t y F l a g
 */
returnValue QProblemB_setInfeasibilityFlag(	QProblemB* THIS,
											returnValue returnvalue, BooleanType doThrowError
											)
{
	THIS->infeasible = BT_TRUE;

	if ( ( doThrowError == BT_TRUE ) || ( THIS->options.enableFarBounds == BT_FALSE ) )
		THROWERROR( returnvalue );

	return returnvalue;
}


/*
 *	a r e B o u n d s C o n s i s t e n t
 */
returnValue QProblemB_areBoundsConsistent(	QProblemB* THIS,
											const real_t* const lb_new, const real_t* const ub_new )
{
	int i;

	if (lb_new && ub_new) {
		for (i = 0; i < QProblemB_getNV(THIS); ++i) {
			if (lb_new[i] > ub_new[i]+QPOASES_EPS) {
				return RET_QP_INFEASIBLE;
			}
		}
	}
	return SUCCESSFUL_RETURN;
}


/*
 *	i s C P U t i m e L i m i t E x c e e d e d
 */
BooleanType QProblemB_isCPUtimeLimitExceeded(	QProblemB* THIS, const real_t* const cputime,
												real_t starttime,
												int nWSR
												)
{
	real_t elapsedTime, timePerIteration;

	/* Always perform next QP iteration if no CPU time limit is given. */
	if ( cputime == 0 )
		return BT_FALSE;

	/* Always perform first QP iteration. */
	if ( nWSR <= 0 )
		return BT_FALSE;

	elapsedTime = qpOASES_getCPUtime( ) - starttime;
	timePerIteration = elapsedTime / ((real_t) nWSR);

	/* Determine if next QP iteration exceed CPU time limit
	 * considering the (current) average CPU time per iteration. */
	if ( ( elapsedTime + timePerIteration*1.25 ) <= ( *cputime ) )
		return BT_FALSE;
	else
		return BT_TRUE;
}


/*
 *	r e g u l a r i s e H e s s i a n
 */
returnValue QProblemB_regulariseHessian( QProblemB* THIS )
{
	/* Do nothing if Hessian regularisation is disbaled! */
	if ( THIS->options.enableRegularisation == BT_FALSE )
		return SUCCESSFUL_RETURN;

	/* Regularisation of identity Hessian not possible. */
	if ( THIS->hessianType == HST_IDENTITY )
		return THROWERROR( RET_CANNOT_REGULARISE_IDENTITY );

	/* Determine regularisation parameter. */
	if ( QProblemB_usingRegularisation( THIS ) == BT_TRUE )
		return THROWERROR( RET_HESSIAN_ALREADY_REGULARISED );
	else
	{
		/* Regularisation of zero Hessian is done implicitly. */
		if ( THIS->hessianType == HST_ZERO )
		{
			THIS->regVal = qpOASES_getNorm( THIS->g,QProblemB_getNV( THIS ),2 ) * THIS->options.epsRegularisation;
		}
		else
		{
			THIS->regVal = DenseMatrix_getNorm( THIS->H,2 ) * THIS->options.epsRegularisation;

			if ( DenseMatrix_addToDiag( THIS->H,THIS->regVal ) == RET_NO_DIAGONAL_AVAILABLE )
				return THROWERROR( RET_CANNOT_REGULARISE_SPARSE );
		}

		THROWINFO( RET_USING_REGULARISATION );
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	p e r f o r m R a t i o T e s t
 */
returnValue QProblemB_performRatioTestB(	QProblemB* THIS, 
											const int nIdx,
											const int* const idxList,
											Bounds* const subjectTo,
											const real_t* const num,
											const real_t* const den,
											real_t epsNum,
											real_t epsDen,
											real_t* t,
											int* BC_idx
											)
{
	int i, ii;

	*BC_idx = -1;

	for( i=0; i<nIdx; ++i )
	{
		ii = idxList[i];

		if ( Bounds_getType( subjectTo,ii ) != ST_EQUALITY )
		{
			if ( ( Bounds_getStatus( subjectTo,ii ) == ST_LOWER ) || ( Bounds_getStatus( subjectTo,ii ) == ST_INACTIVE ) )
			{
				if ( QProblemB_isBlocking( THIS,num[i],den[i],epsNum,epsDen,t ) == BT_TRUE )
				{
					*t = num[i] / den[i];
					*BC_idx = ii;
				}
			}
			else
			if ( Bounds_getStatus( subjectTo,ii ) == ST_UPPER )
			{
				if ( QProblemB_isBlocking( THIS,-num[i],-den[i],epsNum,epsDen,t ) == BT_TRUE )
				{
					*t = num[i] / den[i];
					*BC_idx = ii;
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 * g e t R e l a t i v e H o m o t o p y L e n g t h
 */
real_t QProblemB_getRelativeHomotopyLength(	QProblemB* THIS,
											const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new
											)
{
	int i;
	int nV = QProblemB_getNV( THIS );
	real_t d, s, len = 0.0;

	/* gradient */
	for (i = 0; i < nV; i++)
	{
		s = qpOASES_getAbs(g_new[i]);
		if (s < 1.0) s = 1.0;
		d = qpOASES_getAbs(g_new[i] - THIS->g[i]) / s;
		if (d > len) len = d;
	}

	/* lower bounds */
	for (i = 0; i < nV && lb_new; i++)
	{
		s = qpOASES_getAbs(lb_new[i]);
		if (s < 1.0) s = 1.0;
		d = qpOASES_getAbs(lb_new[i] - THIS->lb[i]) / s;
		if (d > len) len = d;
	}

	/* upper bounds */
	for (i = 0; i < nV && ub_new; i++)
	{
		s = qpOASES_getAbs(ub_new[i]);
		if (s < 1.0) s = 1.0;
		d = qpOASES_getAbs(ub_new[i] - THIS->ub[i]) / s;
		if (d > len) len = d;
	}

	return len;
}


/*
 * u p d a t e F a r B o u n d s
 */
returnValue QProblemB_updateFarBounds(	QProblemB* THIS,
										real_t curFarBound, int nRamp,
                                        const real_t* const lb_new, real_t* const lb_new_far,
                                        const real_t* const ub_new, real_t* const ub_new_far
                                        )
{
	int i;
	real_t rampVal, t;
	int nV = QProblemB_getNV( THIS );

	if ( THIS->options.enableRamping == BT_TRUE )
	{
		for ( i=0; i<nV; ++i )
		{
			t = (real_t)((i + THIS->rampOffset) % nRamp) / (real_t)(nRamp-1);
			rampVal = curFarBound * (1.0 + (1.0-t)*THIS->ramp0 + t*THIS->ramp1);

			if ( lb_new == 0 )
				lb_new_far[i] = -rampVal;
			else
				lb_new_far[i] = qpOASES_getMax( -rampVal,lb_new[i] );

			if ( ub_new == 0 )
				ub_new_far[i] = rampVal;
			else
				ub_new_far[i] = qpOASES_getMin( rampVal,ub_new[i] );
		}
	}
	else
	{
		for ( i=0; i<nV; ++i )
		{
			if ( lb_new == 0 )
				lb_new_far[i] = -curFarBound;
			else
				lb_new_far[i] = qpOASES_getMax( -curFarBound,lb_new[i] );

			if ( ub_new == 0 )
				ub_new_far[i] = curFarBound;
			else
				ub_new_far[i] = qpOASES_getMin( curFarBound,ub_new[i] );
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 * p e r f o r m R a m p i n g
 */
returnValue QProblemB_performRamping( QProblemB* THIS )
{
	int nV = QProblemB_getNV( THIS ), bstat, i;
	real_t t, rampVal;

	/* ramp inactive bounds and active dual variables */
	for (i = 0; i < nV; i++)
	{
		switch (Bounds_getType( &(THIS->bounds),i))
		{
			case ST_EQUALITY: THIS->lb[i] = THIS->x[i]; THIS->ub[i] = THIS->x[i]; continue; /* reestablish exact feasibility */
			case ST_UNBOUNDED: continue;
			case ST_DISABLED: continue;
			default: break;
		}

		t = (real_t)((i + THIS->rampOffset) % nV) / (real_t)(nV-1);
		rampVal = (1.0-t) * THIS->ramp0 + t * THIS->ramp1;
		bstat = Bounds_getStatus(&(THIS->bounds),i);
		if (bstat != ST_LOWER) { THIS->lb[i] = THIS->x[i] - rampVal; }
		if (bstat != ST_UPPER) { THIS->ub[i] = THIS->x[i] + rampVal; }
		if (bstat == ST_LOWER) { THIS->lb[i] = THIS->x[i]; THIS->y[i] = +rampVal; }
		if (bstat == ST_UPPER) { THIS->ub[i] = THIS->x[i]; THIS->y[i] = -rampVal; }
		if (bstat == ST_INACTIVE) THIS->y[i] = 0.0; /* reestablish exact complementarity */
	}

	/* reestablish exact stationarity */
	QProblemB_setupAuxiliaryQPgradient( THIS );

	/* advance ramp offset to avoid Ramping cycles */
	(THIS->rampOffset)++;

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R I V A T E                                                            *
 *****************************************************************************/

/*
 *	s o l v e I n i t i a l Q P
 */
returnValue QProblemB_solveInitialQP(	QProblemB* THIS, const real_t* const xOpt, const real_t* const yOpt,
										Bounds* const guessedBounds,
										int* nWSR, real_t* const cputime
										)
{
	int i;
	int nV = QProblemB_getNV( THIS );
	
	myStatic Bounds auxiliaryBounds;
	
	returnValue returnvalue;
	
	myStatic real_t g_original[NVMAX];
	myStatic real_t lb_original[NVMAX];
	myStatic real_t ub_original[NVMAX];

	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = qpOASES_getCPUtime( );

	BoundsCON( &auxiliaryBounds,nV );


	THIS->status = QPS_NOTINITIALISED;

	/* I) ANALYSE QP DATA: */
	/* 1) Check if Hessian happens to be the identity matrix. */
	if ( QProblemB_determineHessianType( THIS ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup type of bounds (i.e. unbounded, implicitly fixed etc.). */
	if ( QProblemB_setupSubjectToType( THIS ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	THIS->status = QPS_PREPARINGAUXILIARYQP;


	/* II) SETUP AUXILIARY QP WITH GIVEN OPTIMAL SOLUTION: */
	/* 1) Setup bounds data structure. */
	if ( Bounds_setupAllFree( &(THIS->bounds) ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup optimal primal/dual solution for auxiliary QP. */
	if ( QProblemB_setupAuxiliaryQPsolution( THIS,xOpt,yOpt ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 3) Obtain linear independent working set for auxiliary QP. */
	if ( QProblemB_obtainAuxiliaryWorkingSet( THIS,xOpt,yOpt,guessedBounds, &auxiliaryBounds ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 4) Setup working set of auxiliary QP and setup cholesky decomposition. */
	/* a) Working set of auxiliary QP. */
	if ( QProblemB_setupAuxiliaryWorkingSet( THIS,&auxiliaryBounds,BT_TRUE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* b) Regularise Hessian if necessary. */
	if ( ( THIS->hessianType == HST_ZERO ) || ( THIS->hessianType == HST_SEMIDEF ) )
	{
		if ( QProblemB_regulariseHessian( THIS ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_INIT_FAILED_REGULARISATION );
	}

	THIS->haveCholesky = BT_FALSE;
	
	/* 5) Store original QP formulation... */
	for( i=0; i<nV; ++i )
	{
		g_original[i]  = THIS->g[i];
		lb_original[i] = THIS->lb[i];
		ub_original[i] = THIS->ub[i];
	}

	/* ... and setup QP data of an auxiliary QP having an optimal solution
	 * as specified by the user (or xOpt = yOpt = 0, by default). */
	if ( QProblemB_setupAuxiliaryQPgradient( THIS ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	if ( QProblemB_setupAuxiliaryQPbounds( THIS,BT_TRUE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	THIS->status = QPS_AUXILIARYQPSOLVED;


	/* III) SOLVE ACTUAL INITIAL QP: */

	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= qpOASES_getCPUtime( ) - starttime;

	/* Use hotstart method to find the solution of the original initial QP,... */
	returnvalue = QProblemB_hotstart( THIS,g_original,lb_original,ub_original, nWSR,cputime );


	/* ... check for infeasibility and unboundedness... */
	if ( QProblemB_isInfeasible( THIS ) == BT_TRUE )
		return THROWERROR( RET_INIT_FAILED_INFEASIBILITY );

	if ( QProblemB_isUnbounded( THIS ) == BT_TRUE )
		return THROWERROR( RET_INIT_FAILED_UNBOUNDEDNESS );

	/* ... and internal errors. */
	if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
		return THROWERROR( RET_INIT_FAILED_HOTSTART );


	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = qpOASES_getCPUtime( ) - starttime;

	THROWINFO( RET_INIT_SUCCESSFUL );

	return returnvalue;
}


/*
 *	s o l v e Q P
 */
returnValue QProblemB_solveQP(	QProblemB* THIS, const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								int* nWSR, real_t* const cputime, int nWSRperformed
								)
{
	int iter;
	
	/* I) PREPARATIONS */
	/* 1) Allocate delta vectors of gradient and bounds,
	 *    index arrays and step direction arrays. */
	myStatic real_t delta_xFR[NVMAX];
	myStatic real_t delta_xFX[NVMAX];
	myStatic real_t delta_yFX[NVMAX];

	myStatic real_t delta_g[NVMAX];
	myStatic real_t delta_lb[NVMAX];
	myStatic real_t delta_ub[NVMAX];

	returnValue returnvalue;
	BooleanType Delta_bB_isZero;

	int BC_idx;
	SubjectToStatus BC_status;

	real_t homotopyLength;

	#ifndef __XPCTARGET__
	myStatic char messageString[QPOASES_MAX_STRING_LENGTH];
	#endif

	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = qpOASES_getCPUtime( );

	/* consistency check */
	if ( ( QProblemB_getStatus( THIS ) == QPS_NOTINITIALISED )       ||
		 ( QProblemB_getStatus( THIS ) == QPS_PREPARINGAUXILIARYQP ) ||
		 ( QProblemB_getStatus( THIS ) == QPS_PERFORMINGHOMOTOPY )   )
	{
		return THROWERROR( RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED );
	}

	/* 2) Update type of bounds,e.g. a formerly implicitly fixed
	 *    variable might have become a normal one etc. */
	if ( QProblemB_setupSubjectToTypeNew( THIS,lb_new,ub_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_HOTSTART_FAILED );

	/* 3) Reset status flags. */
	THIS->infeasible = BT_FALSE;
	THIS->unbounded  = BT_FALSE;


	/* II) MAIN HOMOTOPY LOOP */
	for( iter=nWSRperformed; iter<*nWSR; ++iter )
	{
		THIS->tabularOutput.idxAddB = THIS->tabularOutput.idxRemB = THIS->tabularOutput.idxAddC = THIS->tabularOutput.idxRemC = -1;
		THIS->tabularOutput.excAddB = THIS->tabularOutput.excRemB = THIS->tabularOutput.excAddC = THIS->tabularOutput.excRemC = 0;

		if ( QProblemB_isCPUtimeLimitExceeded( THIS,cputime,starttime,iter-nWSRperformed ) == BT_TRUE )
		{
			/* Assign number of working set recalculations and stop runtime measurement. */
			*nWSR = iter;
			if ( cputime != 0 )
				*cputime = qpOASES_getCPUtime( ) - starttime;

			break;
		}

		THIS->status = QPS_PERFORMINGHOMOTOPY;

		#ifndef __XPCTARGET__
		snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"%d ...",iter );
		MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_ITERATION_STARTED,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		#endif

		/* 2) Initialise shift direction of the gradient and the bounds. */
		returnvalue = QProblemB_determineDataShift(	THIS,g_new,lb_new,ub_new,
													delta_g,delta_lb,delta_ub,
													&Delta_bB_isZero
													);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			/* Assign number of working set recalculations and stop runtime measurement. */
			*nWSR = iter;
			if ( cputime != 0 )
				*cputime = qpOASES_getCPUtime( ) - starttime;

			THROWERROR( RET_SHIFT_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 3) Determination of step direction of X and Y. */
		returnvalue = QProblemB_determineStepDirection(	THIS,delta_g,delta_lb,delta_ub,
														Delta_bB_isZero,
														delta_xFX,delta_xFR,delta_yFX
														);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			/* Assign number of working set recalculations and stop runtime measurement. */
			*nWSR = iter;
			if ( cputime != 0 )
				*cputime = qpOASES_getCPUtime( ) - starttime;

			THROWERROR( RET_STEPDIRECTION_DETERMINATION_FAILED );
			return returnvalue;
		}


		/* 4) Determination of step length TAU.
		 *    This step along the homotopy path is also taken (without changing working set). */
		returnvalue = QProblemB_performStep(	THIS,delta_g,delta_lb,delta_ub,
												delta_xFX,delta_xFR,delta_yFX,
												&BC_idx,&BC_status
												);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			/* Assign number of working set recalculations and stop runtime measurement. */
			*nWSR = iter;
			if ( cputime != 0 )
				*cputime = qpOASES_getCPUtime( ) - starttime;

			THROWERROR( RET_STEPLENGTH_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 5) Termination criterion. */
		homotopyLength = QProblemB_getRelativeHomotopyLength( THIS,g_new, lb_new, ub_new );
		if ( homotopyLength <= THIS->options.terminationTolerance )
		{
			THIS->status = QPS_SOLVED;

			THROWINFO( RET_OPTIMAL_SOLUTION_FOUND );

			if ( QProblemB_printIteration( THIS,iter,BC_idx,BC_status,homotopyLength ) != SUCCESSFUL_RETURN )
				THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass THIS as return value! */

			*nWSR = iter;

			if ( cputime != 0 )
				*cputime = qpOASES_getCPUtime( ) - starttime;

			return SUCCESSFUL_RETURN;
		}


		/* 6) Change active set. */
		returnvalue = QProblemB_changeActiveSet( THIS,BC_idx,BC_status );

		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			/* Assign number of working set recalculations and stop runtime measurement. */
			*nWSR = iter;
			if ( cputime != 0 )
				*cputime = qpOASES_getCPUtime( ) - starttime;

			/* checks for infeasibility... */
			if ( THIS->infeasible == BT_TRUE )
			{
				THIS->status = QPS_HOMOTOPYQPSOLVED;
				return QProblemB_setInfeasibilityFlag( THIS,RET_HOTSTART_STOPPED_INFEASIBILITY,BT_FALSE );
			}

			/* ...unboundedness... */
			if ( THIS->unbounded == BT_TRUE ) /* not necessary since objective function convex! */
				return THROWERROR( RET_HOTSTART_STOPPED_UNBOUNDEDNESS );

			/* ... and throw unspecific error otherwise */
			THROWERROR( RET_HOMOTOPY_STEP_FAILED );
			return returnvalue;
		}

		/* 6a) Possibly refactorise projected Hessian from scratch. */
		if (THIS->options.enableCholeskyRefactorisation > 0 && iter % THIS->options.enableCholeskyRefactorisation == 0)
		{
			returnvalue = QProblemB_computeCholesky( THIS );
			if (returnvalue != SUCCESSFUL_RETURN)
				return returnvalue;
		}


		/* 7) Perform Ramping Strategy on zero homotopy step or drift correction (if desired). */
		 if ( ( THIS->tau <= QPOASES_EPS ) && ( THIS->options.enableRamping == BT_TRUE ) )
			QProblemB_performRamping( THIS );
		else
		if ( (THIS->options.enableDriftCorrection > 0)
		     && ((iter+1) % THIS->options.enableDriftCorrection == 0) )
			QProblemB_performDriftCorrection( THIS );  /* always returns SUCCESSFUL_RETURN */

		/* 8) Output information of successful QP iteration. */
		THIS->status = QPS_HOMOTOPYQPSOLVED;

		if ( QProblemB_printIteration( THIS,iter,BC_idx,BC_status,homotopyLength ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass THIS as return value! */
	}

	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = qpOASES_getCPUtime( ) - starttime;


	/* if programm gets to here, output information that QP could not be solved
	 * within the given maximum numbers of working set changes */
	if ( THIS->options.printLevel == PL_HIGH )
	{
		#ifndef __XPCTARGET__
		snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"(nWSR = %d)",iter );
		return MessageHandling_throwWarning( qpOASES_getGlobalMessageHandler(),RET_MAX_NWSR_REACHED,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		#else
		return RET_MAX_NWSR_REACHED;
		#endif
	}
	else
	{
		return RET_MAX_NWSR_REACHED;
	}
}


/*
 *	s o l v e R e g u l a r i s e d Q P
 */
returnValue QProblemB_solveRegularisedQP(	QProblemB* THIS, const real_t* const g_new,
											const real_t* const lb_new, const real_t* const ub_new,
											int* nWSR, real_t* const cputime, int nWSRperformed
											)
{
	int i, step;
	int nV = QProblemB_getNV( THIS );

	returnValue returnvalue;

	int nWSR_max   = *nWSR;
	int nWSR_total = nWSRperformed;

	real_t cputime_total = 0.0;
	real_t cputime_cur   = 0.0;

	myStatic real_t gMod[NVMAX];


	/* Perform normal QP solution if QP has not been regularised. */
	if ( QProblemB_usingRegularisation( THIS ) == BT_FALSE )
		return QProblemB_solveQP( THIS,g_new,lb_new,ub_new, nWSR,cputime,nWSRperformed );


	/* I) SOLVE USUAL REGULARISED QP */
	if ( cputime == 0 )
	{
		returnvalue = QProblemB_solveQP( THIS,g_new,lb_new,ub_new, nWSR,0,nWSRperformed );
	}
	else
	{
		cputime_cur = *cputime;
		returnvalue = QProblemB_solveQP( THIS,g_new,lb_new,ub_new, nWSR,&cputime_cur,nWSRperformed );
	}
	nWSR_total     = *nWSR;
	cputime_total += cputime_cur;


	/* Only continue if QP solution has been successful. */
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		if ( cputime != 0 )
			*cputime = cputime_total;

		if ( returnvalue == RET_MAX_NWSR_REACHED )
			THROWWARNING( RET_NO_REGSTEP_NWSR );

		return returnvalue;
	}


	/* II) PERFORM SUCCESSIVE REGULARISATION STEPS */
	for( step=0; step<THIS->options.numRegularisationSteps; ++step )
	{
		/* 1) Modify gradient: gMod = g - eps*xOpt
		 *    (assuming regularisation matrix to be regVal*Id). */
		for( i=0; i<nV; ++i )
			gMod[i] = g_new[i] - THIS->regVal * THIS->x[i];

		/* 2) Solve regularised QP with modified gradient allowing
		 *    only as many working set recalculations and CPU time
		 *    as have been left from previous QP solutions. */
		if ( cputime == 0 )
		{
			*nWSR = nWSR_max;
			returnvalue = QProblemB_solveQP( THIS,gMod,lb_new,ub_new, nWSR,0,nWSR_total );
		}
		else
		{
			*nWSR = nWSR_max;
			cputime_cur = *cputime - cputime_total;
			returnvalue = QProblemB_solveQP( THIS,gMod,lb_new,ub_new, nWSR,&cputime_cur,nWSR_total );
		}

		nWSR_total     = *nWSR;
		cputime_total += cputime_cur;

		/* Only continue if QP solution has been successful. */
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			if ( cputime != 0 )
				*cputime = cputime_total;

			if ( returnvalue == RET_MAX_NWSR_REACHED )
				THROWWARNING( RET_FEWER_REGSTEPS_NWSR );

			return returnvalue;
		}
	}

	if ( cputime != 0 )
		*cputime = cputime_total;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblemB_setupAuxiliaryWorkingSet( 	QProblemB* THIS, 
													Bounds* const auxiliaryBounds,
													BooleanType setupAfresh
													)
{
	int i;
	int nV = QProblemB_getNV( THIS );
	
	BooleanType updateCholesky;

	/* consistency checks */
	if ( auxiliaryBounds != 0 )
	{
		for( i=0; i<nV; ++i )
			if ( ( Bounds_getStatus(&(THIS->bounds),i ) == ST_UNDEFINED ) || ( Bounds_getStatus( auxiliaryBounds,i ) == ST_UNDEFINED ) )
				return THROWERROR( RET_UNKNOWN_BUG );
	}
	else
	{
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}


	/* I) SETUP CHOLESKY FLAG:
	 *    Cholesky decomposition shall only be updated if working set
	 *    shall be updated (i.e. NOT setup afresh!) */
	if ( setupAfresh == BT_TRUE )
		updateCholesky = BT_FALSE;
	else
		updateCholesky = BT_TRUE;


	/* II) REMOVE FORMERLY ACTIVE BOUNDS (IF NECESSARY): */
	if ( setupAfresh == BT_FALSE )
	{
		/* Remove all active bounds that shall be inactive AND
		*  all active bounds that are active at the wrong bound. */
		for( i=0; i<nV; ++i )
		{
			if ( ( Bounds_getStatus(&(THIS->bounds),i ) == ST_LOWER ) && ( Bounds_getStatus( auxiliaryBounds,i ) != ST_LOWER ) )
				if ( QProblemB_removeBound( THIS,i,updateCholesky ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

			if ( ( Bounds_getStatus(&(THIS->bounds),i ) == ST_UPPER ) && ( Bounds_getStatus( auxiliaryBounds,i ) != ST_UPPER ) )
				if ( QProblemB_removeBound( THIS,i,updateCholesky ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}


	/* III) ADD NEWLY ACTIVE BOUNDS: */
	/*      Add all inactive bounds that shall be active AND
	 *      all formerly active bounds that have been active at the wrong bound. */
	for( i=0; i<nV; ++i )
	{
		if ( ( Bounds_getStatus(&(THIS->bounds),i ) == ST_INACTIVE ) && ( Bounds_getStatus( auxiliaryBounds,i ) != ST_INACTIVE ) )
		{
			if ( QProblemB_addBound( THIS,i,Bounds_getStatus( auxiliaryBounds,i ),updateCholesky ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P s o l u t i o n
 */
returnValue QProblemB_setupAuxiliaryQPsolution(	QProblemB* THIS, const real_t* const xOpt, const real_t* const yOpt
												)
{
	int i;
	int nV = QProblemB_getNV( THIS );


	/* Setup primal/dual solution vectors for auxiliary initial QP:
	 * if a null pointer is passed, a zero vector is assigned;
	 * old solution vector is kept if pointer to internal solution vector is passed. */
	if ( xOpt != 0 )
	{
		if ( xOpt != THIS->x )
			for( i=0; i<nV; ++i )
				THIS->x[i] = xOpt[i];
	}
	else
	{
		for( i=0; i<nV; ++i )
			THIS->x[i] = 0.0;
	}

	if ( yOpt != 0 )
	{
		if ( yOpt != THIS->y )
			for( i=0; i<nV; ++i )
				THIS->y[i] = yOpt[i];
	}
	else
	{
		for( i=0; i<nV; ++i )
			THIS->y[i] = 0.0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P g r a d i e n t
 */
returnValue QProblemB_setupAuxiliaryQPgradient( QProblemB* THIS )
{
	int i;
	int nV = QProblemB_getNV( THIS );

	/* Setup gradient vector: g = -H*x + y'*Id. */
	switch ( THIS->hessianType )
	{
		case HST_ZERO:
			if ( QProblemB_usingRegularisation( THIS ) == BT_FALSE )
				for ( i=0; i<nV; ++i )
					THIS->g[i] = THIS->y[i];
			else
				for ( i=0; i<nV; ++i )
					THIS->g[i] = THIS->y[i] - THIS->regVal * THIS->x[i];
			break;

		case HST_IDENTITY:
			for ( i=0; i<nV; ++i )
				THIS->g[i] = THIS->y[i] - THIS->x[i];
			break;

		default:
			/* y'*Id */
			for ( i=0; i<nV; ++i )
				THIS->g[i] = THIS->y[i];

			/* -H*x */
			DenseMatrix_times(THIS->H,1, -1.0, THIS->x, nV, 1.0, THIS->g, nV);

			break;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P b o u n d s
 */
returnValue QProblemB_setupAuxiliaryQPbounds( QProblemB* THIS, BooleanType useRelaxation )
{
	int i;
	int nV = QProblemB_getNV( THIS );


	/* Setup bound vectors. */
	for ( i=0; i<nV; ++i )
	{
		switch ( Bounds_getStatus(&(THIS->bounds),i ) )
		{
			case ST_INACTIVE:
				if ( useRelaxation == BT_TRUE )
				{
					if ( Bounds_getType( &(THIS->bounds),i ) == ST_EQUALITY )
					{
						THIS->lb[i] = THIS->x[i];
						THIS->ub[i] = THIS->x[i];
					}
					else
					{
						THIS->lb[i] = THIS->x[i] - THIS->options.boundRelaxation;
						THIS->ub[i] = THIS->x[i] + THIS->options.boundRelaxation;
					}
				}
				break;

			case ST_LOWER:
				THIS->lb[i] = THIS->x[i];
				if ( Bounds_getType( &(THIS->bounds),i ) == ST_EQUALITY )
				{
					THIS->ub[i] = THIS->x[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						THIS->ub[i] = THIS->x[i] + THIS->options.boundRelaxation;
				}
				break;

			case ST_UPPER:
				THIS->ub[i] = THIS->x[i];
				if ( Bounds_getType( &(THIS->bounds),i ) == ST_EQUALITY )
				{
					THIS->lb[i] = THIS->x[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						THIS->lb[i] = THIS->x[i] - THIS->options.boundRelaxation;
				}
				break;

            case ST_DISABLED:
                break;
                
			default:
				return THROWERROR( RET_UNKNOWN_BUG );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P
 */
returnValue QProblemB_setupAuxiliaryQP( QProblemB* THIS, Bounds* const guessedBounds )
{
	int i;
	int nV = QProblemB_getNV( THIS );

	/* nothing to do */
	if ( guessedBounds == &(THIS->bounds) )
		return SUCCESSFUL_RETURN;

	THIS->status = QPS_PREPARINGAUXILIARYQP;


	/* I) SETUP WORKING SET ... */
	if ( QProblemB_shallRefactorise( THIS,guessedBounds ) == BT_TRUE )
	{
		/* ... WITH REFACTORISATION: */
		/* 1) Reset bounds ... */
		Bounds_init( &(THIS->bounds),nV );

		/*    ... and set them up afresh. */
		if ( QProblemB_setupSubjectToType( THIS ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( Bounds_setupAllFree( &(THIS->bounds) ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 2) Setup guessed working set afresh. */
		if ( QProblemB_setupAuxiliaryWorkingSet( THIS,guessedBounds,BT_TRUE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 3) Calculate Cholesky decomposition. */
		if ( QProblemB_computeCholesky( THIS ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}
	else
	{
		/* ... WITHOUT REFACTORISATION: */
		if ( QProblemB_setupAuxiliaryWorkingSet( THIS,guessedBounds,BT_FALSE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}


	/* II) SETUP AUXILIARY QP DATA: */
	/* 1) Ensure that dual variable is zero for free bounds. */
	for ( i=0; i<nV; ++i )
		if ( Bounds_getStatus(&(THIS->bounds),i ) == ST_INACTIVE )
			THIS->y[i] = 0.0;

	/* 2) Setup gradient and bound vectors. */
	if ( QProblemB_setupAuxiliaryQPgradient( THIS ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	if ( QProblemB_setupAuxiliaryQPbounds( THIS,BT_FALSE ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e S t e p D i r e c t i o n
 */
returnValue QProblemB_determineStepDirection(	QProblemB* THIS, 
												const real_t* const delta_g, const real_t* const delta_lb, const real_t* const delta_ub,
												BooleanType Delta_bB_isZero,
												real_t* const delta_xFX, real_t* const delta_xFR,
												real_t* const delta_yFX
												)
{
	int i, ii;
	int r;
	int nFR = QProblemB_getNFR( THIS );
	int nFX = QProblemB_getNFX( THIS );
	
	int* FR_idx;
	int* FX_idx;

	real_t rnrm;

	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
	
	/* This routine computes
	 * delta_xFX := delta_b
	 * delta_xFR := R \ R' \ -( delta_g + HMX*delta_xFX )
	 * delta_yFX := HMX'*delta_xFR + HFX*delta_xFX  { + eps*delta_xFX }
	 */

	/* I) DETERMINE delta_xFX := delta_{l|u}b */
	if ( Delta_bB_isZero == BT_FALSE )
	{
		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];

			if ( Bounds_getStatus(&(THIS->bounds),ii ) == ST_LOWER )
				delta_xFX[i] = delta_lb[ii];
			else
				delta_xFX[i] = delta_ub[ii];
		}
	}
	else
	{
		for( i=0; i<nFX; ++i )
			delta_xFX[i] = 0.0;
	}


	/* delta_xFR_TMP holds the residual, initialized with right hand side
	 * delta_xFR holds the step that gets refined incrementally */
	for ( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		THIS->delta_xFR_TMP[i] = - delta_g[ii];
		delta_xFR[i] = 0.0;
	}


	/* Iterative refinement loop for delta_xFR */
	for ( r=0; r<=THIS->options.numRefinementSteps; ++r )
	{
		/* II) DETERMINE delta_xFR */
		if ( nFR > 0 )
		{
			/* Add - HMX*delta_xFX
			 * This is skipped if delta_b=0 or mixed part HM=0 (H=0 or H=Id) */
			if ( ( THIS->hessianType != HST_ZERO ) && ( THIS->hessianType != HST_IDENTITY ) && ( Delta_bB_isZero == BT_FALSE ) && ( r == 0 ) )
				DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFixed( &(THIS->bounds) ), 1, -1.0, delta_xFX, nFX, 1.0, THIS->delta_xFR_TMP, nFR, BT_TRUE);

			/* Determine R' \ ( - HMX*delta_xFX - delta_gFR ) where R'R = HFR */
			if ( QProblemB_backsolveR( THIS,THIS->delta_xFR_TMP,BT_TRUE,THIS->delta_xFR_TMP ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );

			/* Determine HFR \ ( - HMX*delta_xFX - delta_gFR ) */
			if ( QProblemB_backsolveR( THIS,THIS->delta_xFR_TMP,BT_FALSE,THIS->delta_xFR_TMP ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );
		}

		/* refine solution found for delta_xFR so far */
		for ( i=0; i<nFR; ++i )
			delta_xFR[i] += THIS->delta_xFR_TMP[i];

		if ( THIS->options.numRefinementSteps > 0 )
		{
			rnrm = 0.0;
			/* compute new residual in delta_xFR_TMP:
			 * residual := - HFR*delta_xFR - HMX*delta_xFX - delta_gFR
			 * set to -delta_gFR */
			for ( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				THIS->delta_xFR_TMP[i] = -delta_g[ii];
			}
			/* add - HFR*delta_xFR */
			switch ( THIS->hessianType )
			{
				case HST_ZERO:
					break;

				case HST_IDENTITY:
					for ( i=0; i<nFR; ++i )
					{
						THIS->delta_xFR_TMP[i] -= delta_xFR[i];

						/* compute max norm */
						if (rnrm < qpOASES_getAbs (THIS->delta_xFR_TMP[i]))
							rnrm = qpOASES_getAbs (THIS->delta_xFR_TMP[i]);
					}
					break;

				default:
					DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFree( &(THIS->bounds) ),  1, -1.0, delta_xFR, nFR, 1.0, THIS->delta_xFR_TMP, nFR, BT_TRUE);
					DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFixed( &(THIS->bounds) ), 1, -1.0, delta_xFX, nFX, 1.0, THIS->delta_xFR_TMP, nFR, BT_TRUE);

					/* compute max norm */
					for ( i=0; i<nFR; ++i )
						if (rnrm < qpOASES_getAbs (THIS->delta_xFR_TMP[i]))
							rnrm = qpOASES_getAbs (THIS->delta_xFR_TMP[i]);

					break;
			}
			
			/* early termination of residual norm small enough */
			if ( rnrm < THIS->options.epsIterRef )
				break;
		}

	} /* end of refinement loop for delta_xFR */

	/* III) DETERMINE delta_yFX */
	if ( nFX > 0 )
	{
		if ( ( THIS->hessianType == HST_ZERO ) || ( THIS->hessianType == HST_IDENTITY ) )
		{
			for( i=0; i<nFX; ++i )
			{
				/* set to delta_g */
				ii = FX_idx[i];
				delta_yFX[i] = delta_g[ii];

				/* add HFX*delta_xFX = {0|I}*delta_xFX */
				if ( THIS->hessianType == HST_ZERO )
				{
					if ( QProblemB_usingRegularisation( THIS ) == BT_TRUE )
						delta_yFX[i] += THIS->regVal*delta_xFX[i];
				}
				else
					delta_yFX[i] += delta_xFX[i];
			}
		}
		else
		{
			for( i=0; i<nFX; ++i )
			{
				/* set to delta_g */
				ii = FX_idx[i];
				delta_yFX[i] = delta_g[ii];
			}
			DenseMatrix_subTimes(THIS->H,Bounds_getFixed( &(THIS->bounds) ), Bounds_getFree( &(THIS->bounds) ), 1, 1.0, delta_xFR, nFR, 1.0, delta_yFX, nFX, BT_TRUE);
			if (Delta_bB_isZero == BT_FALSE)
				DenseMatrix_subTimes(THIS->H,Bounds_getFixed( &(THIS->bounds) ), Bounds_getFixed( &(THIS->bounds) ), 1, 1.0, delta_xFX, nFX, 1.0, delta_yFX, nFX, BT_TRUE);
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p e r f o r m S t e p
 */
returnValue QProblemB_performStep(	QProblemB* THIS, 
									const real_t* const delta_g,
									const real_t* const delta_lb, const real_t* const delta_ub,
									const real_t* const delta_xFX,
									const real_t* const delta_xFR,
									const real_t* const delta_yFX,
									int* BC_idx, SubjectToStatus* BC_status
									)
{
	int i, ii;
	int nV = QProblemB_getNV( THIS );
	int nFR = QProblemB_getNFR( THIS );
	int nFX = QProblemB_getNFX( THIS );
	
	myStatic char messageString[QPOASES_MAX_STRING_LENGTH];

	int BC_idx_tmp = -1;

	myStatic real_t num[NVMAX];
	myStatic real_t den[NVMAX];
	
	int* FR_idx;
	int* FX_idx;

	THIS->tau = 1.0;
	*BC_idx = -1;
	*BC_status = ST_UNDEFINED;

	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );


	/* I) DETERMINE MAXIMUM DUAL STEPLENGTH, i.e. ensure that
	 *    active dual bounds remain valid (ignoring implicitly fixed variables): */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = THIS->y[ii];
		den[i] = -delta_yFX[i];
	}

	QProblemB_performRatioTestB( THIS,nFX,FX_idx,&(THIS->bounds),num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

	if ( BC_idx_tmp >= 0 )
	{
		*BC_idx = BC_idx_tmp;
		*BC_status = ST_INACTIVE;
	}


	/* II) DETERMINE MAXIMUM PRIMAL STEPLENGTH, i.e. ensure that
	 *     inactive bounds remain valid (ignoring unbounded variables). */
	/* 1) Inactive lower bounds. */
	if ( Bounds_hasNoLower( &(THIS->bounds) ) == BT_FALSE )
	{
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			num[i] = qpOASES_getMax( THIS->x[ii] - THIS->lb[ii],0.0 );
			den[i] = delta_lb[ii] - delta_xFR[i];
		}

		QProblemB_performRatioTestB( THIS,nFR,FR_idx,&(THIS->bounds),num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			*BC_idx = BC_idx_tmp;
			*BC_status = ST_LOWER;
		}
	}

	/* 2) Inactive upper bounds. */
	if ( Bounds_hasNoUpper( &(THIS->bounds) ) == BT_FALSE )
	{
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			num[i] = qpOASES_getMax( THIS->ub[ii] - THIS->x[ii],0.0 );
			den[i] = delta_xFR[i] - delta_ub[ii];
		}

		QProblemB_performRatioTestB( THIS,nFR,FR_idx,&(THIS->bounds),num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			*BC_idx = BC_idx_tmp;
			*BC_status = ST_UPPER;
		}
	}


	#ifndef __XPCTARGET__
	if ( *BC_status == ST_UNDEFINED )
		snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"Stepsize is %.15e!",THIS->tau );
	else
		snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"Stepsize is %.15e! (idx = %d, status = %d)",THIS->tau,*BC_idx,*BC_status );

	MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_STEPSIZE_NONPOSITIVE,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
	#endif


	/* III) PERFORM STEP ALONG HOMOTOPY PATH */
	if ( THIS->tau > QPOASES_ZERO )
	{
		/* 1) Perform step in primal und dual space. */
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			THIS->x[ii] += THIS->tau * delta_xFR[i];
		}

		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];
			THIS->x[ii] += THIS->tau * delta_xFX[i];
			THIS->y[ii] += THIS->tau * delta_yFX[i];
		}

		/* 2) Shift QP data. */
		for( i=0; i<nV; ++i )
		{
			THIS->g[i]  += THIS->tau * delta_g[i];
			THIS->lb[i] += THIS->tau * delta_lb[i];
			THIS->ub[i] += THIS->tau * delta_ub[i];
		}
	}
	else
	{
		/* print a warning if stepsize is zero */
		#ifndef __XPCTARGET__
		snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"Stepsize is %.15e",THIS->tau );
		MessageHandling_throwWarning( qpOASES_getGlobalMessageHandler(),RET_STEPSIZE,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		#endif
	}


	return SUCCESSFUL_RETURN;
}


/*
 *	c h a n g e A c t i v e S e t
 */
returnValue QProblemB_changeActiveSet( QProblemB* THIS, int BC_idx, SubjectToStatus BC_status )
{
	#ifndef __XPCTARGET__
	myStatic char messageString[QPOASES_MAX_STRING_LENGTH];
	#endif

	/* IV) UPDATE ACTIVE SET */
	switch ( BC_status )
	{
		/* Optimal solution found as no working set change detected. */
		case ST_UNDEFINED:
			return RET_OPTIMAL_SOLUTION_FOUND;


		/* Remove one variable from active set. */
		case ST_INACTIVE:
			#ifndef __XPCTARGET__
			snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"bound no. %d.", BC_idx );
			MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( QProblemB_removeBound( THIS,BC_idx,BT_TRUE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_REMOVE_FROM_ACTIVESET_FAILED );

			THIS->y[BC_idx] = 0.0;
			break;


		/* Add one variable to active set. */
		default:
			#ifndef __XPCTARGET__
			if ( BC_status == ST_LOWER )
				snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"lower bound no. %d.", BC_idx );
			else
				snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"upper bound no. %d.", BC_idx );
				MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_ADD_TO_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( QProblemB_addBound( THIS,BC_idx,BC_status,BT_TRUE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_ADD_TO_ACTIVESET_FAILED );
			break;
	}

	return SUCCESSFUL_RETURN;
}



/*
 * p e r f o r m D r i f t C o r r e c t i o n
 */
returnValue QProblemB_performDriftCorrection( QProblemB* THIS )
{
	int i;
	int nV = QProblemB_getNV( THIS );

	for ( i=0; i<nV; ++i )
	{
		switch ( Bounds_getType ( &(THIS->bounds),i ) )
		{
			case ST_BOUNDED:
				switch ( Bounds_getStatus( &(THIS->bounds),i ) )
				{
					case ST_LOWER:
						THIS->lb[i] = THIS->x[i];
						THIS->ub[i] = qpOASES_getMax (THIS->ub[i], THIS->x[i]);
						THIS->y[i] = qpOASES_getMax (THIS->y[i], 0.0);
						break;
					case ST_UPPER:
						THIS->lb[i] = qpOASES_getMin (THIS->lb[i], THIS->x[i]);
						THIS->ub[i] = THIS->x[i];
						THIS->y[i] = qpOASES_getMin (THIS->y[i], 0.0);
						break;
					case ST_INACTIVE:
						THIS->lb[i] = qpOASES_getMin (THIS->lb[i], THIS->x[i]);
						THIS->ub[i] = qpOASES_getMax (THIS->ub[i], THIS->x[i]);
						THIS->y[i] = 0.0;
						break;
					case ST_UNDEFINED:
					case ST_INFEASIBLE_LOWER:
					case ST_INFEASIBLE_UPPER:
						break;
				}
				break;
			case ST_EQUALITY:
				THIS->lb[i] = THIS->x[i];
				THIS->ub[i] = THIS->x[i];
				break;
			case ST_UNBOUNDED:
			case ST_UNKNOWN:
            case ST_DISABLED:
				break;
		}
	}

	return QProblemB_setupAuxiliaryQPgradient( THIS );
}


/*
 *	s h a l l R e f a c t o r i s e
 */
BooleanType QProblemB_shallRefactorise( QProblemB* THIS, Bounds* const guessedBounds )
{
	int i;
	int nV = QProblemB_getNV( THIS );
	
	int differenceNumber = 0;

	/* always refactorise if Hessian is not known to be positive definite */
	if ( ( THIS->hessianType == HST_SEMIDEF ) || ( THIS->hessianType == HST_INDEF ) )
		return BT_TRUE;

	/* 1) Determine number of bounds that have same status
	 *    in guessed AND current bounds.*/
	for( i=0; i<nV; ++i )
		if ( Bounds_getStatus( guessedBounds,i ) != Bounds_getStatus(&(THIS->bounds),i ) )
			++differenceNumber;

	/* 2) Decide wheter to refactorise or not. */
	if ( 2*differenceNumber > Bounds_getNFX( guessedBounds ) )
		return BT_TRUE;
	else
		return BT_FALSE;
}


/*
 *	a d d B o u n d
 */
returnValue QProblemB_addBound(	QProblemB* THIS, 
								int number, SubjectToStatus B_status,
								BooleanType updateCholesky
								)
{
	int i, j;
	int nFR = QProblemB_getNFR( THIS );

	int number_idx;
	real_t c, s, nu;

	/* consistency check */
	if ( ( QProblemB_getStatus( THIS ) == QPS_NOTINITIALISED )    ||
		 ( QProblemB_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblemB_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblemB_getStatus( THIS ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* Perform cholesky updates only if QProblemB has been initialised! */
	if ( QProblemB_getStatus( THIS ) == QPS_PREPARINGAUXILIARYQP )
	{
		/* UPDATE INDICES */
		if ( Bounds_moveFreeToFixed( &(THIS->bounds),number,B_status ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_ADDBOUND_FAILED );

		return SUCCESSFUL_RETURN;
	}


	/* I) PERFORM CHOLESKY UPDATE: */
	if ( ( updateCholesky == BT_TRUE ) &&
		 ( THIS->hessianType != HST_ZERO )   && ( THIS->hessianType != HST_IDENTITY ) )
	{
		/* 1) Index of variable to be added within the list of free variables. */
		number_idx = Indexlist_getIndex( Bounds_getFree( &(THIS->bounds) ),number );

		/* 2) Use row-wise Givens rotations to restore upper triangular form of R. */
		for( i=number_idx+1; i<nFR; ++i )
		{
			QProblemB_computeGivens( RR(i-1,i),RR(i,i), &RR(i-1,i),&RR(i,i),&c,&s );
			nu = s/(1.0+c);

			for( j=(1+i); j<nFR; ++j ) /* last column of R is thrown away */
				QProblemB_applyGivens( c,s,nu,RR(i-1,j),RR(i,j), &RR(i-1,j),&RR(i,j) );
		}

		/* 3) Delete <number_idx>th column and ... */
		for( i=0; i<nFR-1; ++i )
			for( j=number_idx+1; j<nFR; ++j )
				RR(i,j-1) = RR(i,j);
		/* ... last column of R. */
		for( i=0; i<nFR; ++i )
			RR(i,nFR-1) = 0.0;
	}

	/* II) UPDATE INDICES */
	THIS->tabularOutput.idxAddB = number;
	if ( Bounds_moveFreeToFixed( &(THIS->bounds),number,B_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDBOUND_FAILED );


	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e B o u n d
 */
returnValue QProblemB_removeBound(	QProblemB* THIS, 
									int number,
									BooleanType updateCholesky
									)
{
	int i;
	int nFR = QProblemB_getNFR( THIS );
	
	int* FR_idx;

	myStatic real_t rhs[NVMAX+1];
	myStatic real_t r[NVMAX];

	real_t r0;


	/* consistency check */
	if ( ( QProblemB_getStatus( THIS ) == QPS_NOTINITIALISED )    ||
		 ( QProblemB_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblemB_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblemB_getStatus( THIS ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* I) UPDATE INDICES */
	THIS->tabularOutput.idxRemB = number;
	if ( Bounds_moveFixedToFree( &(THIS->bounds),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_REMOVEBOUND_FAILED );

	/* Perform cholesky updates only if QProblemB has been initialised! */
	if ( QProblemB_getStatus( THIS ) == QPS_PREPARINGAUXILIARYQP )
		return SUCCESSFUL_RETURN;


	/* II) PERFORM CHOLESKY UPDATE */
	if ( ( updateCholesky == BT_TRUE ) &&
		 ( THIS->hessianType != HST_ZERO )   && ( THIS->hessianType != HST_IDENTITY ) )
	{
		Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );

		/* 1) Calculate new column of cholesky decomposition. */
		switch ( THIS->hessianType )
		{
			case HST_ZERO:
				if ( QProblemB_usingRegularisation( THIS ) == BT_FALSE )
					r0 = 0.0;
				else
					r0 = THIS->regVal;
				for( i=0; i<nFR; ++i )
					rhs[i] = 0.0;
				break;

			case HST_IDENTITY:
				r0 = 1.0;
				for( i=0; i<nFR; ++i )
					rhs[i] = 0.0;
				break;

			default:
				DenseMatrix_getRow(THIS->H,number, Bounds_getFree( &(THIS->bounds) ), 1.0, rhs);
				r0 = DenseMatrix_diag(THIS->H,number);
				break;
		}

		if ( QProblemB_backsolveRrem( THIS,rhs,BT_TRUE,BT_TRUE,r ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_REMOVEBOUND_FAILED );

		for( i=0; i<nFR; ++i )
			r0 -= r[i]*r[i];

		/* 2) Store new column into R. */
		for( i=0; i<nFR; ++i )
			RR(i,nFR) = r[i];

		if ( r0 > QPOASES_ZERO )
			RR(nFR,nFR) = qpOASES_getSqrt( r0 );
		else
		{
			THIS->hessianType = HST_SEMIDEF;
			return THROWERROR( RET_HESSIAN_NOT_SPD );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t I t e r a t i o n
 */
returnValue QProblemB_printIteration( 	QProblemB* THIS, 
										int iter,
										int BC_idx,	SubjectToStatus BC_status,
										real_t homotopyLength
		  								)
{
	#ifndef __XPCTARGET__
	myStatic char myPrintfString[QPOASES_MAX_STRING_LENGTH];
	myStatic char info[QPOASES_MAX_STRING_LENGTH];

	/* consistency check */
	if ( iter < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* nothing to do */
	if ( THIS->options.printLevel != PL_MEDIUM )
		return SUCCESSFUL_RETURN;


	/* 1) Print header at first iteration. */
 	if ( iter == 0 )
	{
		snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"\n\n#################   qpOASES  --  QP NO. %3.0d   ##################\n\n", THIS->count );
		qpOASES_myPrintf( myPrintfString );

		qpOASES_myPrintf( "    Iter   |    StepLength    |       Info       |   nFX    \n" );
		qpOASES_myPrintf( " ----------+------------------+------------------+--------- \n" );
	}

	/* 2) Print iteration line. */
	if ( BC_status == ST_UNDEFINED )
	{
		if ( THIS->hessianType == HST_ZERO )
			snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |    LP SOLVED     |  %4.1d   \n", iter,THIS->tau,QProblemB_getNFX( THIS ) );
		else
			snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |    QP SOLVED     |  %4.1d   \n", iter,THIS->tau,QProblemB_getNFX( THIS ) );
		qpOASES_myPrintf( myPrintfString );
	}
	else
	{
		if ( BC_status == ST_INACTIVE )
			snprintf( info,8,"REM BND" );
		else
			snprintf( info,8,"ADD BND" );

		snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |   %s %4.1d   |  %4.1d   \n", iter,THIS->tau,info,BC_idx,QProblemB_getNFX( THIS ) );
		qpOASES_myPrintf( myPrintfString );
	}
	#endif

	return SUCCESSFUL_RETURN;
}



END_NAMESPACE_QPOASES


/*
 *	end of file
 */
