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
 *	\file src/QProblem.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0embedded
 *	\date 2007-2014
 *
 *	Implementation of the QProblem class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming.
 */


#include <qpOASES/QProblem.h>
#include <qpOASES/QProblemB.h>



BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

/*
 *	Q P r o b l e m
 */
void QProblemCON(	QProblem* THIS,
					int _nV, int _nC, HessianType _hessianType )
{
	#ifdef __CODE_GENERATION__
	Options_setToFast( &(THIS->options) );
	#else
	Options_setToDefault( &(THIS->options) );
	#endif /* __CODE_GENERATION__ */

	/* print copyright notice */
	if (THIS->options.printLevel != PL_NONE)
		qpOASES_printCopyrightNotice( );

	/* consistency checks */
	if ( ( _nV <= 0 ) || ( _nV > NVMAX ) )
	{
		_nV = 1;
		THROWERROR( RET_INVALID_ARGUMENTS );
		assert( 1 == 0 );
	}

	if ( ( _nC < 0 ) || ( _nC > NCMAX ) )
	{
		_nC = 0;
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

	QProblem_setPrintLevel( THIS,THIS->options.printLevel );

	THIS->A = &(THIS->AA);

	Constraints_init( &(THIS->constraints),_nC );

	THIS->sizeT = qpOASES_getMinI( _nV,_nC );

	THIS->constraintProduct = 0;
}



/*
 *	c o p y
 */
void QProblemCPY(	QProblem* FROM,
					QProblem* TO
					)
{
	unsigned int _nV = (unsigned int)QProblem_getNV( FROM );
	unsigned int _nC = (unsigned int)QProblem_getNC( FROM );

	TO->bounds = FROM->bounds;

	TO->HH = FROM->HH;
	TO->H = &(TO->HH);

	QProblem_setG( TO,FROM->g );
	QProblem_setLB( TO,FROM->lb );
	QProblem_setUB( TO,FROM->ub );

	memcpy( TO->R,FROM->R,NVMAX*NVMAX*sizeof(real_t) );
	
	TO->haveCholesky = FROM->haveCholesky;

	memcpy( TO->x,FROM->x,_nV*sizeof(real_t) );

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
	QProblem_setPrintLevel( TO,TO->options.printLevel );
	
	ConstraintsCPY( &(FROM->constraints),&(TO->constraints) );

	TO->AA = FROM->AA;
	TO->A = &(TO->AA);

	QProblem_setLBA( TO,FROM->lbA );
	QProblem_setUBA( TO,FROM->ubA );

	memcpy( TO->y,FROM->y,(_nV+_nC)*sizeof(real_t) );

	TO->sizeT = FROM->sizeT;

	memcpy( TO->T,FROM->T,NVCMIN*NVCMIN*sizeof(real_t) );
	memcpy( TO->Q,FROM->Q,NVMAX*NVMAX*sizeof(real_t) );

	memcpy( TO->Ax,FROM->Ax,_nC*sizeof(real_t) );
	memcpy( TO->Ax_l,FROM->Ax_l,_nC*sizeof(real_t) );
	memcpy( TO->Ax_u,FROM->Ax_u,_nC*sizeof(real_t) );

	if ( FROM->constraintProduct != 0 )
		TO->constraintProduct = FROM->constraintProduct;
	else
		TO->constraintProduct = 0;
}



/*
 *	r e s e t
 */
returnValue QProblem_reset( QProblem* THIS )
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

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

	/* 2) Reset constraints. */
	Constraints_init( &(THIS->constraints),nC );

	/* 3) Reset TQ factorisation. */
	for( i=0; i<NVCMIN*NVCMIN; ++i )
		THIS->T[i] = 0.0;

	for( i=0; i<NVMAX*NVMAX; ++i )
		THIS->Q[i] = 0.0;

	/* 4) Reset constraint product pointer. */
	THIS->constraintProduct = 0;

	return SUCCESSFUL_RETURN;
}


/*
 *	i n i t
 */
returnValue QProblem_initM(	QProblem* THIS, DenseMatrix *_H, const real_t* const _g, DenseMatrix *_A,
							const real_t* const _lb, const real_t* const _ub,
							const real_t* const _lbA, const real_t* const _ubA,
							int* nWSR, real_t* const cputime
							)
{
	int nV,nC;

	if ( QProblem_getNV( THIS ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( QProblem_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblem_reset( THIS );
	}

	nV = QProblem_getNV( THIS );
	nC = QProblem_getNC( THIS );

	/*DenseMatrix_print( _H );
	DenseMatrix_print( _A );

	qpOASES_printV( _g,nV );
	qpOASES_printV( _lb,nV );
	qpOASES_printV( _ub,nV );
	qpOASES_printV( _lbA,nC );
	qpOASES_printV( _ubA,nC );*/


	/* 2) Setup QP data. */
	if ( QProblem_setupQPdataM( THIS,_H,_g,_A,_lb,_ub,_lbA,_ubA ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine (without any additional information). */
	return QProblem_solveInitialQP( THIS,0,0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblem_init(	QProblem* THIS, real_t* const _H, const real_t* const _g, real_t* const _A,
							const real_t* const _lb, const real_t* const _ub,
							const real_t* const _lbA, const real_t* const _ubA,
							int* nWSR, real_t* const cputime
							)
{
	if ( QProblem_getNV( THIS ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( QProblem_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblem_reset( THIS );
	}

	/* 2) Setup QP data. */
	if ( QProblem_setupQPdata( THIS,_H,_g,_A,_lb,_ub,_lbA,_ubA ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine (without any additional information). */
	return QProblem_solveInitialQP( THIS,0,0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblem_initF(	QProblem* THIS, const char* const H_file, const char* const g_file, const char* const A_file,
							const char* const lb_file, const char* const ub_file,
							const char* const lbA_file, const char* const ubA_file,
							int* nWSR, real_t* const cputime
							)
{
	if ( QProblem_getNV( THIS ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( QProblem_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblem_reset( THIS );
	}

	/* 2) Setup QP data from files. */
	if ( QProblem_setupQPdataFromFile( THIS,H_file,g_file,A_file,lb_file,ub_file,lbA_file,ubA_file ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* 3) Call to main initialisation routine (without any additional information). */
	return QProblem_solveInitialQP( THIS,0,0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblem_initMW(	QProblem* THIS, DenseMatrix *_H, const real_t* const _g, DenseMatrix *_A,
								const real_t* const _lb, const real_t* const _ub,
								const real_t* const _lbA, const real_t* const _ubA,
								int* nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								Bounds* const guessedBounds, Constraints* const guessedConstraints
								)
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( QProblem_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblem_reset( THIS );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( Bounds_getStatus( guessedBounds,i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	if ( guessedConstraints != 0 )
	{
		for( i=0; i<nC; ++i )
			if ( Constraints_getStatus( guessedConstraints,i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* exclude these possibilities in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data. */
	if ( QProblem_setupQPdataM( THIS,_H,_g,_A,_lb,_ub,_lbA,_ubA ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return QProblem_solveInitialQP( THIS,xOpt,yOpt,guessedBounds,guessedConstraints, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblem_initW(	QProblem* THIS, real_t* const _H, const real_t* const _g, real_t* const _A,
							const real_t* const _lb, const real_t* const _ub,
							const real_t* const _lbA, const real_t* const _ubA,
							int* nWSR, real_t* const cputime,
							const real_t* const xOpt, const real_t* const yOpt,
							Bounds* const guessedBounds, Constraints* const guessedConstraints
							)
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( QProblem_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblem_reset( THIS );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( Bounds_getStatus( guessedBounds,i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	if ( guessedConstraints != 0 )
	{
		for( i=0; i<nC; ++i )
			if ( Constraints_getStatus( guessedConstraints,i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* exclude these possibilities in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data. */
	if ( QProblem_setupQPdata( THIS,_H,_g,_A,_lb,_ub,_lbA,_ubA ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return QProblem_solveInitialQP( THIS,xOpt,yOpt,guessedBounds,guessedConstraints, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblem_initFW(	QProblem* THIS, const char* const H_file, const char* const g_file, const char* const A_file,
								const char* const lb_file, const char* const ub_file,
								const char* const lbA_file, const char* const ubA_file,
								int* nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								Bounds* const guessedBounds, Constraints* const guessedConstraints
								)
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( QProblem_isInitialised( THIS ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		QProblem_reset( THIS );
	}

	for( i=0; i<nV; ++i )
	{
		if ( Bounds_getStatus( guessedBounds,i ) == ST_UNDEFINED )
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	for( i=0; i<nC; ++i )
		if ( Constraints_getStatus( guessedConstraints,i ) == ST_UNDEFINED )
			return THROWERROR( RET_INVALID_ARGUMENTS );

	/* exclude these possibilities in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data from files. */
	if ( QProblem_setupQPdataFromFile( THIS,H_file,g_file,A_file,lb_file,ub_file,lbA_file,ubA_file ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* 3) Call to main initialisation routine. */
	return QProblem_solveInitialQP( THIS,xOpt,yOpt,guessedBounds,guessedConstraints, nWSR,cputime );
}


/*
 * s e t u p I n i t i a l C h o l e s k y
 */
returnValue QProblem_setupInitialCholesky( QProblem* THIS )
{
	returnValue returnvalueCholesky;

	/* If regularisation shall be used, always regularise at beginning 
	 * if initial working set is not empty. */
	if ( ( QProblem_getNV(THIS) != QProblem_getNFR(THIS) - QProblem_getNFV(THIS) ) && ( THIS->options.enableRegularisation == BT_TRUE ) )
		if ( QProblem_regulariseHessian( THIS ) != SUCCESSFUL_RETURN )
			return RET_INIT_FAILED_REGULARISATION;
	
	/* Factorise projected Hessian 
	 * now handles all special cases (no active bounds/constraints, no nullspace) */
	returnvalueCholesky = QProblem_computeProjectedCholesky( THIS );

	/* If Hessian is not positive definite, regularise and try again. */
	if ( returnvalueCholesky == RET_HESSIAN_NOT_SPD )
	{
		if ( QProblem_regulariseHessian( THIS ) != SUCCESSFUL_RETURN )
			return RET_INIT_FAILED_REGULARISATION;

		returnvalueCholesky = QProblem_computeProjectedCholesky( THIS );
	}

	if ( returnvalueCholesky != SUCCESSFUL_RETURN )
		return RET_INIT_FAILED_CHOLESKY;

	THIS->haveCholesky = BT_TRUE;
	return SUCCESSFUL_RETURN;
}


/*
 *	h o t s t a r t
 */
returnValue QProblem_hotstart(	QProblem* THIS,
								const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								const real_t* const lbA_new, const real_t* const ubA_new,
								int* nWSR, real_t* const cputime
								)
{
	returnValue returnvalue = SUCCESSFUL_RETURN;
	int i, nActiveFar;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	int nWSR_max = *nWSR;
	int nWSR_performed = 0;

	real_t cputime_remaining = QPOASES_INFTY, *pcputime_rem;
	real_t cputime_needed = 0.0;

	real_t farbound = THIS->options.initialFarBounds;

	myStatic real_t ub_new_far[NVMAX];
	myStatic real_t lb_new_far[NVMAX];
	myStatic real_t ubA_new_far[NCMAX];
	myStatic real_t lbA_new_far[NCMAX];
	
	real_t tol;
	
	if ( QProblem_getNV( THIS ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* Simple check for consistency of bounds and constraints. */
	if ( QProblem_areBoundsConsistent(THIS,lb_new, ub_new, lbA_new, ubA_new) != SUCCESSFUL_RETURN )
		return QProblem_setInfeasibilityFlag(THIS,returnvalue,BT_TRUE);

	++(THIS->count);

    /*qpOASES_printNV( g_new,nV,"g" );
    qpOASES_printNV( lb_new,nV,"lb" );
    qpOASES_printNV( ub_new,nV,"ub" );
    qpOASES_printNV( lbA_new,nC,"lbA" );
    qpOASES_printNV( ubA_new,nC,"ubA" );*/
    
	/*QProblem_writeQpDataIntoMatFile( THIS,"qpData.mat" );*/
	/*QProblem_writeQpWorkspaceIntoMatFile( THIS,"qpWorkspace.mat" );*/


	if ( THIS->haveCholesky == BT_FALSE )
	{
		returnvalue = QProblem_setupInitialCholesky( THIS );
		if (returnvalue != SUCCESSFUL_RETURN)
			return THROWERROR(returnvalue);
	}


	if ( THIS->options.enableFarBounds == BT_FALSE )
	{
		/* Automatically call standard solveQP if regularisation is not active. */
		returnvalue = QProblem_solveRegularisedQP( THIS,g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,cputime,0 );
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
		if (ubA_new)
			for (i = 0; i < nC; i++)
				if ((ubA_new[i] < QPOASES_INFTY) && (ubA_new[i] > farbound)) farbound = ubA_new[i];
		if (lbA_new)
			for (i = 0; i < nC; i++)
				if ((lbA_new[i] > -QPOASES_INFTY) && (lbA_new[i] < -farbound)) farbound = -lbA_new[i];

		QProblem_updateFarBounds(	THIS,farbound,nV+nC,
									lb_new,lb_new_far, ub_new,ub_new_far,
									lbA_new,lbA_new_far, ubA_new,ubA_new_far
									);

		for ( ;; )
		{
			*nWSR = nWSR_max;
			if ( cputime != 0 )
			{
				cputime_remaining = *cputime - cputime_needed;
				pcputime_rem = &cputime_remaining;
			}
			else
				pcputime_rem = 0;

			/* Automatically call standard solveQP if regularisation is not active. */
			returnvalue = QProblem_solveRegularisedQP( THIS,g_new,lb_new_far,ub_new_far,lbA_new_far,ubA_new_far, nWSR,pcputime_rem,nWSR_performed );

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
					break; /* goto farewell; */
				}

				QProblem_updateFarBounds(	THIS,farbound,nV+nC,
											lb_new,lb_new_far, ub_new,ub_new_far,
											lbA_new,lbA_new_far, ubA_new,ubA_new_far
											);
			}
			else if ( THIS->status == QPS_SOLVED )
			{
				tol = farbound/THIS->options.growFarBounds * THIS->options.boundTolerance;
				
				for ( i=0; i<nV; ++i )
				{
					if ( ( ( lb_new == 0 ) || ( lb_new_far[i] > lb_new[i] ) ) && ( qpOASES_getAbs ( lb_new_far[i] - THIS->x[i] ) < tol ) )
						++nActiveFar;
					if ( ( ( ub_new == 0 ) || ( ub_new_far[i] < ub_new[i] ) ) && ( qpOASES_getAbs ( ub_new_far[i] - THIS->x[i] ) < tol ) )
						++nActiveFar;
				}
				for ( i=0; i<nC; ++i )
				{
					if ( ( ( lbA_new == 0 ) || ( lbA_new_far[i] > lbA_new[i] ) ) && ( qpOASES_getAbs ( lbA_new_far[i] - THIS->Ax[i] ) < tol ) )
						++nActiveFar;
					if ( ( ( ubA_new == 0 ) || ( ubA_new_far[i] < ubA_new[i] ) ) && ( qpOASES_getAbs ( ubA_new_far[i] - THIS->Ax[i] ) < tol ) )
						++nActiveFar;
				}

				if ( nActiveFar == 0 )
					break;

				THIS->status = QPS_HOMOTOPYQPSOLVED;

				if ( farbound >= QPOASES_INFTY )
				{
					THIS->unbounded = BT_TRUE;
					returnvalue = RET_HOTSTART_STOPPED_UNBOUNDEDNESS;
					goto farewell;
				}

				QProblem_updateFarBounds(	THIS,farbound,nV+nC,
											lb_new,lb_new_far, ub_new,ub_new_far,
											lbA_new,lbA_new_far, ubA_new,ubA_new_far
											);
			}
			else
			{
				/* some other error */
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
returnValue QProblem_hotstartF(	QProblem* THIS, const char* const g_file,
								const char* const lb_file, const char* const ub_file,
								const char* const lbA_file, const char* const ubA_file,
								int* nWSR, real_t* const cputime
								)
{
	int nV = QProblem_getNV( THIS );
	returnValue returnvalue;

	/* 1) Allocate memory (if bounds exist). */
	myStatic real_t g_new[NVMAX];
	myStatic real_t lb_new[NVMAX];
	myStatic real_t ub_new[NVMAX];
	myStatic real_t lbA_new[NCMAX];
	myStatic real_t ubA_new[NCMAX];


	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* consistency check */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 2) Load new QP vectors from file. */
	returnvalue = QProblem_loadQPvectorsFromFile(	THIS,g_file,lb_file,ub_file,lbA_file,ubA_file,
													g_new,lb_new,ub_new,lbA_new,ubA_new
													);
	if ( returnvalue != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* 3) Actually perform hotstart. */
	returnvalue = QProblem_hotstart( THIS,g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,cputime );

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblem_hotstartW(	QProblem* THIS, const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								const real_t* const lbA_new, const real_t* const ubA_new,
								int* nWSR, real_t* const cputime,
								Bounds* const guessedBounds, Constraints* const guessedConstraints
								)
{
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );
	
	real_t starttime = 0.0;
	myStatic Bounds emptyBounds;
	myStatic Constraints emptyConstraints;
	
	returnValue returnvalue;

	BoundsCON( &emptyBounds,nV );
	ConstraintsCON( &emptyConstraints,nC );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );


	/* start runtime measurement */
	if ( cputime != 0 )
		starttime = qpOASES_getCPUtime( );


	/* 1) Update working sets according to guesses for working sets of bounds and constraints. */
	if ( ( guessedBounds != 0 ) && ( guessedConstraints != 0 ) )
	{
		if ( QProblem_setupAuxiliaryQP( THIS,guessedBounds,guessedConstraints ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}

	if ( ( guessedBounds == 0 ) && ( guessedConstraints != 0 ) )
	{
		/* create empty bounds for setting up auxiliary QP */
		if ( Bounds_setupAllFree( &emptyBounds ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( QProblem_setupAuxiliaryQP( THIS,&emptyBounds,guessedConstraints ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}

	if ( ( guessedBounds != 0 ) && ( guessedConstraints == 0 ) )
	{
		/* create empty constraints for setting up auxiliary QP */
		if ( Constraints_setupAllInactive( &emptyConstraints ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( QProblem_setupAuxiliaryQP( THIS,guessedBounds,&emptyConstraints ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}

	if ( ( guessedBounds == 0 ) && ( guessedConstraints == 0 ) )
	{
		/* create empty bounds and constraints for setting up auxiliary QP */
		if ( Bounds_setupAllFree( &emptyBounds ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
		if ( Constraints_setupAllInactive( &emptyConstraints ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( QProblem_setupAuxiliaryQP( THIS,&emptyBounds,&emptyConstraints ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}

	THIS->status = QPS_AUXILIARYQPSOLVED;


	/* 2) Perform usual homotopy. */

	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= qpOASES_getCPUtime( ) - starttime;

	returnvalue = QProblem_hotstart( THIS,g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,cputime );

	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = qpOASES_getCPUtime( ) - starttime;

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblem_hotstartFW(	QProblem* THIS, const char* const g_file,
									const char* const lb_file, const char* const ub_file,
									const char* const lbA_file, const char* const ubA_file,
									int* nWSR, real_t* const cputime,
									Bounds* const guessedBounds, Constraints* const guessedConstraints
									)
{
	int nV = QProblem_getNV( THIS );
	returnValue returnvalue;
	
	/* 1) Allocate memory (if bounds exist). */
	myStatic real_t g_new[NVMAX];
	myStatic real_t lb_new[NVMAX];
	myStatic real_t ub_new[NVMAX];
	myStatic real_t lbA_new[NCMAX];
	myStatic real_t ubA_new[NCMAX];


	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* consistency check */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Load new QP vectors from file. */
	returnvalue = QProblem_loadQPvectorsFromFile(	THIS,g_file,lb_file,ub_file,lbA_file,ubA_file,
													g_new,lb_new,ub_new,lbA_new,ubA_new
													);
	if (returnvalue != SUCCESSFUL_RETURN) {
		returnvalue = RET_UNABLE_TO_READ_FILE;
		goto farewell;
	}
	
	/* 3) Actually perform hotstart using initialised homotopy. */
	returnvalue = QProblem_hotstartW(	THIS,g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,cputime,
										guessedBounds,guessedConstraints
										);
farewell:
	return ( returnvalue != SUCCESSFUL_RETURN ) ? THROWERROR ( returnvalue ) : returnvalue;
}


/*
 *
 */
returnValue QProblem_solveCurrentEQP(	QProblem* THIS,
										const int n_rhs,
										const real_t* g_in,
										const real_t* lb_in,
										const real_t* ub_in,
										const real_t* lbA_in,
										const real_t* ubA_in,
										real_t* x_out,
										real_t* y_out
										)
{
	returnValue returnvalue = SUCCESSFUL_RETURN;
	int ii, jj;
	int nV  = QProblem_getNV( THIS );
	int nC  = QProblem_getNC( THIS );
	int nFR = QProblem_getNFR( THIS );
	int nFX = QProblem_getNFX( THIS );
	int nAC = QProblem_getNAC( THIS );

	myStatic real_t delta_xFX[NVMAX];
	myStatic real_t delta_xFR[NVMAX];
	myStatic real_t delta_yAC[NCMAX];
	myStatic real_t delta_yFX[NVMAX];

	/* 1) Determine index arrays. */
	int* FR_idx;
	int* FX_idx;
	int* AC_idx;

	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );

	for ( ii = 0 ; ii < (nV+nC)*n_rhs; ++ii )
		y_out[ii] = 0.0;

	for ( ii = 0 ; ii < n_rhs; ++ii )
	{
		returnvalue = QProblem_determineStepDirection( THIS,
						g_in, lbA_in, ubA_in, lb_in, ub_in, BT_FALSE, BT_FALSE,
						delta_xFX, delta_xFR, delta_yAC, delta_yFX );

		for ( jj = 0; jj < nFX; ++jj )
			x_out[FX_idx[jj]] = delta_xFX[jj];
		for ( jj = 0; jj < nFR; ++jj )
			x_out[FR_idx[jj]] = delta_xFR[jj];
		for ( jj = 0; jj < nFX; ++jj )
			y_out[FX_idx[jj]] = delta_yFX[jj];
		for ( jj = 0; jj < nAC; ++jj )
			y_out[nV+AC_idx[jj]] = delta_yAC[jj];

		g_in += nV;
		lb_in += nV;
		ub_in += nV;
		lbA_in += nC;
		ubA_in += nC;
		x_out += nV;
		y_out += nV+nC;
	}

	return returnvalue;
}



/*
 *	g e t W o r k i n g S e t
 */
returnValue QProblem_getWorkingSet( QProblem* THIS, real_t* workingSet )
{
	int i;
	int nC = QProblem_getNC( THIS );
	int nV = QProblem_getNV( THIS );

	/* At which limit is the bounds/contraints active? */
	for (i = 0; i < nV; i++) {
		switch ( Bounds_getStatus( &(THIS->bounds),i ) ) {
		case ST_LOWER:
			workingSet[i] = -1.0;
			break;
		case ST_UPPER:
			workingSet[i] = +1.0;
			break;
		default:
			workingSet[i] =  0.0;
			break;
		}
	}
	for (i = 0; i < nC; i++) {
		switch ( Constraints_getStatus( &(THIS->constraints),i ) ) {
		case ST_LOWER:
			workingSet[i + nV] = -1.0;
			break;
		case ST_UPPER:
			workingSet[i + nV] = +1.0;
			break;
		default:
			workingSet[i + nV] =  0.0;
			break;
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	g e t N Z
 */
int QProblem_getNZ( QProblem* THIS )
{
	/* nZ = nFR - nAC */
	return QProblem_getNFR( THIS ) - QProblem_getNAC( THIS );
}


/*
 *	g e t D u a l S o l u t i o n
 */
returnValue QProblem_getDualSolution( QProblem* THIS, real_t* const yOpt )
{
	int i;
	
	for( i=0; i<QProblem_getNV( THIS )+QProblem_getNC( THIS ); ++i )
		yOpt[i] = THIS->y[i];

	/* return optimal dual solution vector
	 * only if current QP has been solved */
	if ( ( QProblem_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblem_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblem_getStatus( THIS ) == QPS_SOLVED ) )
	{
		return SUCCESSFUL_RETURN;
	}
	else
	{
		return RET_QP_NOT_SOLVED;
	}
}



/*
 *	s e t C o n s t r a i n t P r o d u c t
 */
returnValue QProblem_setConstraintProduct( QProblem* THIS, ConstraintProduct* const _constraintProduct )
{
	THIS->constraintProduct = _constraintProduct;

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t O b j V a l
 */
real_t QProblem_getObjVal( QProblem* THIS )
{
	real_t objVal;

	/* calculated optimal objective function value
	 * only if current QP has been solved */
	if ( ( QProblem_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblem_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblem_getStatus( THIS ) == QPS_SOLVED ) )
	{
		objVal = QProblem_getObjValX( THIS,THIS->x );
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
real_t QProblem_getObjValX( QProblem* THIS, const real_t* const _x )
{
	int i;
	int nV = QProblem_getNV( THIS );

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
	if ( QProblem_usingRegularisation( THIS ) == BT_TRUE )
	{
		for( i=0; i<nV; ++i )
			objVal += 0.5*_x[i]*THIS->regVal*_x[i];
	}

	return objVal;
}


/*
 *	g e t P r i m a l S o l u t i o n
 */
returnValue QProblem_getPrimalSolution( QProblem* THIS, real_t* const xOpt )
{
	int i;

	/* return optimal primal solution vector
	 * only if current QP has been solved */
	if ( ( QProblem_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblem_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblem_getStatus( THIS ) == QPS_SOLVED ) )
	{
		for( i=0; i<QProblem_getNV( THIS ); ++i )
			xOpt[i] = THIS->x[i];

		return SUCCESSFUL_RETURN;
	}
	else
	{
		return RET_QP_NOT_SOLVED;
	}
}


/*
 *	s e t P r i n t L e v e l
 */
returnValue QProblem_setPrintLevel( QProblem* THIS, PrintLevel _printLevel )
{
	#ifndef __XPCTARGET__
		#ifndef __MATLAB__
		if ( ( THIS->options.printLevel == PL_HIGH ) && ( THIS->options.printLevel != _printLevel ) )
			THROWINFO( RET_PRINTLEVEL_CHANGED );
		#endif /* __MATLAB__ */
		THIS->options.printLevel = _printLevel;
	#else
	THIS->options.printLevel = PL_NONE;
	#endif /* __XPCTARGET__ */

	/* update message handler preferences */
 	switch ( THIS->options.printLevel )
 	{
 		case PL_NONE:
 			MessageHandling_setErrorVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_HIDDEN );
			MessageHandling_setWarningVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_HIDDEN );
			MessageHandling_setInfoVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_HIDDEN );
			break;

		case PL_DEBUG_ITER:
		case PL_TABULAR:
		case PL_LOW:
			MessageHandling_setErrorVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_VISIBLE );
			MessageHandling_setWarningVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_HIDDEN );
			MessageHandling_setInfoVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_HIDDEN );
			break;

		case PL_MEDIUM:
			MessageHandling_setErrorVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_VISIBLE );
			MessageHandling_setWarningVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_VISIBLE );
			MessageHandling_setInfoVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_HIDDEN );
			break;

		default: /* PL_HIGH */
			MessageHandling_setErrorVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_VISIBLE );
			MessageHandling_setWarningVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_VISIBLE );
			MessageHandling_setInfoVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_VISIBLE );
			break;
 	}

	return SUCCESSFUL_RETURN;
}



returnValue QProblem_printOptions( QProblem* THIS )
{
	return Options_print( &(THIS->options) );
}


/*
 *	p r i n t P r o p e r t i e s
 */
returnValue QProblem_printProperties( QProblem* THIS )
{
	#ifndef __XPCTARGET__
	myStatic char myPrintfString[QPOASES_MAX_STRING_LENGTH];
	
	/* Do not print properties if print level is set to none! */
	if ( THIS->options.printLevel == PL_NONE )
		return SUCCESSFUL_RETURN;

	qpOASES_myPrintf( "\n#################   qpOASES  --  QP PROPERTIES   #################\n" );
	qpOASES_myPrintf( "\n" );

	/* 1) Variables properties. */
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,  "Number of Variables: %4.1d\n",QProblem_getNV( THIS ) );
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


	/* 2) Constraints properties. */
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,  "Total number of Constraints:      %4.1d\n",QProblem_getNC( THIS ) );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,  "Number of Equality Constraints:   %4.1d\n",QProblem_getNEC( THIS ) );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,  "Number of Inequality Constraints: %4.1d\n",QProblem_getNC( THIS )-QProblem_getNEC( THIS ) );
	qpOASES_myPrintf( myPrintfString );

	if ( QProblem_getNC( THIS ) > 0 )
	{
		if ( Constraints_hasNoLower( &(THIS->constraints) ) == BT_TRUE )
				qpOASES_myPrintf( "Constraints are not bounded from below.\n" );
			else
				qpOASES_myPrintf( "Constraints are bounded from below.\n" );

		if ( Constraints_hasNoUpper( &(THIS->constraints) ) == BT_TRUE )
				qpOASES_myPrintf( "Constraints are not bounded from above.\n" );
			else
				qpOASES_myPrintf( "Constraints are bounded from above.\n" );
	}

	qpOASES_myPrintf( "\n" );


	/* 3) Further properties. */
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


	/* 4) QP object properties. */
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

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/


/*
 *	d e t e r m i n e H e s s i a n T y p e
 */
returnValue QProblem_determineHessianType( QProblem* THIS )
{
	int i;
	int nV = QProblem_getNV( THIS );
	
	BooleanType isIdentity, isZero;
	

	/* if Hessian type has been set by user, do NOT change it! */
	if ( THIS->hessianType != HST_UNKNOWN )
		return SUCCESSFUL_RETURN;

	/* if Hessian has not been allocated, assume it to be all zeros! */
	if ( THIS->H == 0 )
	{
		THIS->hessianType = HST_ZERO;
		return SUCCESSFUL_RETURN;
	}


	/* 1) If Hessian has outer-diagonal elements,
	 *    Hessian is assumed to be positive definite. */
	THIS->hessianType = HST_POSDEF;
	if (DenseMatrix_isDiag(THIS->H) == BT_FALSE)
		return SUCCESSFUL_RETURN;

	/* 2) Otherwise it is diagonal and test for identity or zero matrix is performed. */
	/* THIS->hessianType = HST_DIAGONAL; */

	isIdentity = BT_TRUE;
	isZero = BT_TRUE;

	for ( i=0; i<nV; ++i )
	{
		if ( THIS->options.enableFlippingBounds == BT_FALSE )
			if ( DenseMatrix_diag(THIS->H,i) < -QPOASES_ZERO ) 
				return THROWERROR( RET_HESSIAN_INDEFINITE );

		if ( qpOASES_getAbs( DenseMatrix_diag(THIS->H,i) - 1.0 ) > QPOASES_EPS )
			isIdentity = BT_FALSE;

		if ( qpOASES_getAbs( DenseMatrix_diag(THIS->H,i) ) > QPOASES_EPS )
			isZero = BT_FALSE;
	}

	if ( isIdentity == BT_TRUE )
		THIS->hessianType = HST_IDENTITY;

	if ( isZero == BT_TRUE )
		THIS->hessianType = HST_ZERO;

	return SUCCESSFUL_RETURN;
}


/*
 *	c o m p u t e C h o l e s k y
 */
returnValue QProblemBCPY_computeCholesky( QProblem* THIS )
{
	int i, j;
	int nV  = QProblem_getNV( THIS );
	int nFR = QProblem_getNFR( THIS );
	
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
			/* if Hessian is zero matrix, it is assumed that it has been
			 * regularised and thus its Cholesky factor is the identity
			 * matrix scaled by sqrt(eps). */
			if ( QProblem_usingRegularisation( THIS ) == BT_TRUE )
			{
				for( i=0; i<nV; ++i )
					RR(i,i) = qpOASES_getSqrt( THIS->regVal );
			}
			else
				return THROWERROR( RET_CHOLESKY_OF_ZERO_HESSIAN );
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
returnValue QProblemBCPY_obtainAuxiliaryWorkingSet(	QProblem* THIS, const real_t* const xOpt, const real_t* const yOpt,
													Bounds* const guessedBounds, Bounds* auxiliaryBounds
													)
{
	int i = 0;
	int nV = QProblem_getNV( THIS );


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
returnValue QProblem_backsolveR(	QProblem* THIS, const real_t* const b, BooleanType transposed,
									real_t* const a
									)
{
	/* Call standard backsolve procedure (i.e. removingBound == BT_FALSE). */
	return QProblem_backsolveRrem( THIS,b,transposed,BT_FALSE,a );
}


/*
 *	b a c k s o l v e R
 */
returnValue QProblem_backsolveRrem(	QProblem* THIS, const real_t* const b, BooleanType transposed,
									BooleanType removingBound,
									real_t* const a
									)
{
	int i, j;
	int nR = QProblem_getNZ( THIS );

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
returnValue QProblemBCPY_determineDataShift(	QProblem* THIS, const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new,
												real_t* const delta_g, real_t* const delta_lb, real_t* const delta_ub,
												BooleanType* Delta_bB_isZero
												)
{
	int i, ii;
	int nV  = QProblem_getNV( THIS );
	int nFX = QProblem_getNFX( THIS );

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
returnValue QProblemBCPY_setupQPdataM(	QProblem* THIS, DenseMatrix *_H, const real_t* const _g,
										const real_t* const _lb, const real_t* const _ub
										)
{
	if ( _H == 0 )
		return QProblemBCPY_setupQPdata( THIS,(real_t*)0,_g,_lb,_ub );
	else
		return QProblemBCPY_setupQPdata( THIS,DenseMatrix_getVal(_H),_g,_lb,_ub );
}


/*
 *	s e t u p Q P d a t a
 */
returnValue QProblemBCPY_setupQPdata(	QProblem* THIS, real_t* const _H, const real_t* const _g,
										const real_t* const _lb, const real_t* const _ub
										)
{
	int i;
	int nV = QProblem_getNV( THIS );

	/* 1) Setup Hessian matrix. */
	QProblem_setH( THIS,_H );

	/* 2) Setup gradient vector. */
	if ( _g == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );
	else
		QProblem_setG( THIS,_g );

	/* 3) Setup lower bounds vector. */
	if ( _lb != 0 )
	{
		QProblem_setLB( THIS,_lb );
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
		QProblem_setUB( THIS,_ub );
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
returnValue QProblemBCPY_setupQPdataFromFile(	QProblem* THIS, const char* const H_file, const char* const g_file,
												const char* const lb_file, const char* const ub_file
												)
{
	int i;
	int nV = QProblem_getNV( THIS );

	returnValue returnvalue;


	/* 1) Load Hessian matrix from file. */
	myStatic real_t _H[NVMAX*NVMAX];

	if ( H_file != 0 )
	{
		returnvalue = qpOASES_readFromFileM( _H, nV,nV, H_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
		QProblem_setH( THIS,_H );
	}
	else
	{
		QProblem_setH( THIS,(real_t*)0 );
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
returnValue QProblemBCPY_loadQPvectorsFromFile(	QProblem* THIS, const char* const g_file, const char* const lb_file, const char* const ub_file,
												real_t* const g_new, real_t* const lb_new, real_t* const ub_new
												)
{
	int nV = QProblem_getNV( THIS );

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
returnValue QProblem_setInfeasibilityFlag(	QProblem* THIS,
											returnValue returnvalue, BooleanType doThrowError
											)
{
	THIS->infeasible = BT_TRUE;

	if ( THIS->options.enableFarBounds == BT_FALSE )
		THROWERROR( returnvalue );

	return returnvalue;
}


/*
 *	i s C P U t i m e L i m i t E x c e e d e d
 */
BooleanType QProblem_isCPUtimeLimitExceeded(	QProblem* THIS, const real_t* const cputime,
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
returnValue QProblem_regulariseHessian( QProblem* THIS )
{
	/* Do nothing if Hessian regularisation is disbaled! */
	if ( THIS->options.enableRegularisation == BT_FALSE )
		return SUCCESSFUL_RETURN;

	/* Regularisation of identity Hessian not possible. */
	if ( THIS->hessianType == HST_IDENTITY )
		return THROWERROR( RET_CANNOT_REGULARISE_IDENTITY );

	/* Determine regularisation parameter. */
	if ( QProblem_usingRegularisation( THIS ) == BT_TRUE )
		return SUCCESSFUL_RETURN; /*THROWERROR( RET_HESSIAN_ALREADY_REGULARISED );*/
	else
	{
		/* Regularisation of zero Hessian is done implicitly. */
		if ( THIS->hessianType == HST_ZERO )
		{
			THIS->regVal = qpOASES_getNorm( THIS->g,QProblem_getNV( THIS ),2 ) * THIS->options.epsRegularisation;
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
returnValue QProblem_performRatioTestB(	QProblem* THIS, 
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
				if ( QProblem_isBlocking( THIS,num[i],den[i],epsNum,epsDen,t ) == BT_TRUE )
				{
					*t = num[i] / den[i];
					*BC_idx = ii;
				}
			}
			else
			if ( Bounds_getStatus( subjectTo,ii ) == ST_UPPER )
			{
				if ( QProblem_isBlocking( THIS,-num[i],-den[i],epsNum,epsDen,t ) == BT_TRUE )
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
 * r e l a t i v e H o m o t o p y L e n g t h
 */
real_t QProblemBCPY_getRelativeHomotopyLength(	QProblem* THIS,
												const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new
												)
{
	int i;
	int nV = QProblem_getNV( THIS );
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
 *	s o l v e I n i t i a l Q P
 */
returnValue QProblem_solveInitialQP(	QProblem* THIS,
										const real_t* const xOpt, const real_t* const yOpt,
										Bounds* const guessedBounds, Constraints* const guessedConstraints,
										int* nWSR, real_t* const cputime
										)
{
	int i;
	#ifdef __MANY_CONSTRAINTS__
	int j;
	real_t l1;
	#endif

	/* some definitions */
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );
	
	myStatic Bounds auxiliaryBounds;
	myStatic Constraints auxiliaryConstraints;
	
	myStatic real_t g_original[NVMAX];
	myStatic real_t lb_original[NVMAX];
	myStatic real_t ub_original[NVMAX];
	myStatic real_t lbA_original[NCMAX];
	myStatic real_t ubA_original[NCMAX];
	
	returnValue returnvalue;


	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = qpOASES_getCPUtime( );


	/*DenseMatrix_print( THIS->H );
	DenseMatrix_print( THIS->A );

	qpOASES_printV( THIS->g,nV );
	qpOASES_printV( THIS->lb,nV );
	qpOASES_printV( THIS->ub,nV );
	qpOASES_printV( THIS->lbA,nC );
	qpOASES_printV( THIS->ubA,nC );*/


	THIS->status = QPS_NOTINITIALISED;

	BoundsCON( &auxiliaryBounds,nV );
	ConstraintsCON( &auxiliaryConstraints,nC );

	/* I) ANALYSE QP DATA: */
	/* 1) Check if Hessian happens to be the identity matrix. */
	if ( QProblem_determineHessianType( THIS ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup type of bounds and constraints (i.e. unbounded, implicitly fixed etc.). */
	if ( QProblem_setupSubjectToType( THIS ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	THIS->status = QPS_PREPARINGAUXILIARYQP;


	/* II) SETUP AUXILIARY QP WITH GIVEN OPTIMAL SOLUTION: */
	/* 1) Setup bounds and constraints data structure. */
	if ( Bounds_setupAllFree( &(THIS->bounds) ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	if ( Constraints_setupAllInactive( &(THIS->constraints) ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup optimal primal/dual solution for auxiliary QP. */
	if ( QProblem_setupAuxiliaryQPsolution( THIS,xOpt,yOpt ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 3) Obtain linear independent working set for auxiliary QP. */
	if ( QProblem_obtainAuxiliaryWorkingSet(	THIS,xOpt,yOpt,guessedBounds,guessedConstraints,
												&auxiliaryBounds,&auxiliaryConstraints ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 4) Setup working set of auxiliary QP and setup matrix factorisations. */
	/* a) Regularise Hessian if necessary. */
	if ( ( THIS->hessianType == HST_ZERO ) || ( THIS->hessianType == HST_SEMIDEF ) )
	{
		if ( QProblem_regulariseHessian( THIS ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_INIT_FAILED_REGULARISATION );
	}

	/* b) TQ factorisation. */
	if ( QProblem_setupTQfactorisation( THIS ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED_TQ );

	/* c) Working set of auxiliary QP. */
	if ( QProblem_setupAuxiliaryWorkingSet( THIS,&auxiliaryBounds,&auxiliaryConstraints,BT_TRUE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	THIS->haveCholesky = BT_FALSE;

	/* 5) Store original QP formulation... */
	for( i=0; i<nV; ++i )
	{
		g_original[i] = THIS->g[i];
		lb_original[i] = THIS->lb[i];
		ub_original[i] = THIS->ub[i];
	}

	for( i=0; i<nC; ++i )
	{
		lbA_original[i] = THIS->lbA[i];
		ubA_original[i] = THIS->ubA[i];
	}

	/* ... and setup QP data of an auxiliary QP having an optimal solution
	 * as specified by the user (or xOpt = yOpt = 0, by default). */
	if ( QProblem_setupAuxiliaryQPgradient( THIS ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	if ( QProblem_setupAuxiliaryQPbounds( THIS,&auxiliaryBounds,&auxiliaryConstraints,BT_TRUE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	THIS->status = QPS_AUXILIARYQPSOLVED;


	if ( THIS->options.enableRamping == BT_TRUE )
		QProblem_performRamping( THIS );


	/* III) SOLVE ACTUAL INITIAL QP: */
	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= qpOASES_getCPUtime( ) - starttime;

	/* Use hotstart method to find the solution of the original initial QP,... */
	returnvalue = QProblem_hotstart( THIS,g_original,lb_original,ub_original,lbA_original,ubA_original, nWSR,cputime );

	/* ... check for infeasibility and unboundedness... */
	if ( QProblem_isInfeasible( THIS ) == BT_TRUE )
		return THROWERROR( RET_INIT_FAILED_INFEASIBILITY );

	if ( QProblem_isUnbounded( THIS ) == BT_TRUE )
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
returnValue QProblem_solveQP(	QProblem* THIS,
								const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								const real_t* const lbA_new, const real_t* const ubA_new,
								int* nWSR, real_t* const cputime, int nWSRperformed
								)
{
	int iter;

	/* I) PREPARATIONS */
	/* 1) Allocate delta vectors of gradient and (constraints') bounds,
	 *    index arrays and step direction arrays. */

	myStatic real_t delta_xFR[NVMAX];
	myStatic real_t delta_xFX[NVMAX];
	myStatic real_t delta_yAC[NCMAX];
	myStatic real_t delta_yFX[NVMAX];

	myStatic real_t delta_g[NVMAX];
	myStatic real_t delta_lb[NVMAX];
	myStatic real_t delta_ub[NVMAX];
	myStatic real_t delta_lbA[NCMAX];
	myStatic real_t delta_ubA[NCMAX];

	returnValue returnvalue;
	BooleanType Delta_bC_isZero, Delta_bB_isZero;

	int BC_idx;
	SubjectToStatus BC_status;
	BooleanType BC_isBound;

	real_t homotopyLength;

	#ifndef __XPCTARGET__
	myStatic char messageString[QPOASES_MAX_STRING_LENGTH];
	#endif

	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = qpOASES_getCPUtime( );


	/* consistency check */
	if ( ( QProblem_getStatus( THIS ) == QPS_NOTINITIALISED )       ||
		 ( QProblem_getStatus( THIS ) == QPS_PREPARINGAUXILIARYQP ) ||
		 ( QProblem_getStatus( THIS ) == QPS_PERFORMINGHOMOTOPY )   )
	{
		return THROWERROR( RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED );
	}

   
	/* 2) Update type of bounds and constraints, e.g.
	 *    a former equality constraint might have become a normal one etc. */
    
/*  (ckirches) disabled THIS, as inactive but tight bounds may become inactive equalities
               which would then never become active again!*/
/*    
	if ( setupSubjectToType( THIS,lb_new,ub_new,lbA_new,ubA_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_HOTSTART_FAILED );
*/
	/* 3) Reset status flags. */
	THIS->infeasible = BT_FALSE;
	THIS->unbounded  = BT_FALSE;


	/* II) MAIN HOMOTOPY LOOP */
	for( iter=nWSRperformed; iter<*nWSR; ++iter )
	{
		THIS->tabularOutput.idxAddB = THIS->tabularOutput.idxRemB = THIS->tabularOutput.idxAddC = THIS->tabularOutput.idxRemC = -1;
		THIS->tabularOutput.excAddB = THIS->tabularOutput.excRemB = THIS->tabularOutput.excAddC = THIS->tabularOutput.excRemC = 0;

		if ( QProblem_isCPUtimeLimitExceeded( THIS,cputime,starttime,iter-nWSRperformed ) == BT_TRUE )
		{
			/* If CPU time limit is exceeded, stop homotopy loop immediately!
			* Assign number of working set recalculations (runtime measurement is stopped later). */
			*nWSR = iter;
			break;
		}

		/*qpOASES_printV( THIS->x,QProblem_getNV(THIS) );*/

		THIS->status = QPS_PERFORMINGHOMOTOPY;

		#ifndef __XPCTARGET__
		snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"%d ...",iter );
		MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_ITERATION_STARTED,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		#endif


		/* 2) Determination of shift direction of the gradient and the (constraints') bounds. */
		returnvalue = QProblem_determineDataShift(	THIS,g_new,lbA_new,ubA_new,lb_new,ub_new,
													delta_g,delta_lbA,delta_ubA,delta_lb,delta_ub,
													&Delta_bC_isZero,&Delta_bB_isZero
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
		returnvalue = QProblem_determineStepDirection(	THIS,delta_g,delta_lbA,delta_ubA,delta_lb,delta_ub,
														Delta_bC_isZero, Delta_bB_isZero,
														delta_xFX,delta_xFR,delta_yAC,delta_yFX
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
		returnvalue = QProblem_performStep(	THIS,delta_g, delta_lbA,delta_ubA,delta_lb,delta_ub,
											delta_xFX,delta_xFR,delta_yAC,delta_yFX,
											&BC_idx,&BC_status,&BC_isBound
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
		homotopyLength = QProblem_getRelativeHomotopyLength( THIS,g_new,lb_new,ub_new,lbA_new,ubA_new );
		/* fprintf( stdFile, "*homotopyLength* = %.3e\n",homotopyLength ); */

		if ( homotopyLength <= THIS->options.terminationTolerance )
		{
			THIS->status = QPS_SOLVED;

			THROWINFO( RET_OPTIMAL_SOLUTION_FOUND );

			if ( QProblem_printIteration( THIS,iter,BC_idx,BC_status,BC_isBound,homotopyLength ) != SUCCESSFUL_RETURN )
				THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass THIS as return value! */

			*nWSR = iter;
			if ( cputime != 0 )
				*cputime = qpOASES_getCPUtime( ) - starttime;

			return SUCCESSFUL_RETURN;
		}


		/* 6) Change active set. */
		returnvalue = QProblem_changeActiveSet( THIS,BC_idx,BC_status,BC_isBound );
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			/* Assign number of working set recalculations and stop runtime measurement. */
			*nWSR = iter;
			if ( cputime != 0 )
				*cputime = qpOASES_getCPUtime( ) - starttime;

			/* Checks for infeasibility... */
			if ( QProblem_isInfeasible( THIS ) == BT_TRUE )
			{				
				THIS->status = QPS_HOMOTOPYQPSOLVED;
				return QProblem_setInfeasibilityFlag( THIS,RET_HOTSTART_STOPPED_INFEASIBILITY,BT_FALSE );
			}

			/* ...unboundedness... */
			if ( QProblem_isUnbounded( THIS ) == BT_TRUE ) /* not necessary since objective function convex! */
				return THROWERROR( RET_HOTSTART_STOPPED_UNBOUNDEDNESS );

			/* ... and throw unspecific error otherwise */
			THROWERROR( RET_HOMOTOPY_STEP_FAILED );
			return returnvalue;
		}

		/* 6a) Possibly refactorise projected Hessian from scratch. */
		if (THIS->options.enableCholeskyRefactorisation > 0 && iter % THIS->options.enableCholeskyRefactorisation == 0)
		{
			returnvalue = QProblem_computeProjectedCholesky( THIS );
			if (returnvalue != SUCCESSFUL_RETURN)
				return returnvalue;
		}

		/* 7) Output information of successful QP iteration. */
		THIS->status = QPS_HOMOTOPYQPSOLVED;

		if ( QProblem_printIteration( THIS,iter,BC_idx,BC_status,BC_isBound,homotopyLength ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass THIS as return value! */

		/* 8) Perform Ramping Strategy on zero homotopy step or drift correction (if desired). */
		if (BC_status != ST_UNDEFINED)
		{
			if ( ( THIS->tau <= QPOASES_EPS ) && ( THIS->options.enableRamping == BT_TRUE ) )
				QProblem_performRamping( THIS );
			else
			if ( (THIS->options.enableDriftCorrection > 0)
			  && ((iter+1) % THIS->options.enableDriftCorrection == 0) )
				QProblem_performDriftCorrection( THIS );  /* always returns SUCCESSFUL_RETURN */
		}
	}

	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = qpOASES_getCPUtime( ) - starttime;


	/* if program gets to here, output information that QP could not be solved
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
returnValue QProblem_solveRegularisedQP(	QProblem* THIS,
											const real_t* const g_new,
											const real_t* const lb_new, const real_t* const ub_new,
											const real_t* const lbA_new, const real_t* const ubA_new,
											int* nWSR, real_t* const cputime, int nWSRperformed
											)
{
	int i, step;
	int nV = QProblem_getNV( THIS );

	returnValue returnvalue;

	int nWSR_max   = *nWSR;
	int nWSR_total = nWSRperformed;

	real_t cputime_total = 0.0;
	real_t cputime_cur   = 0.0;
	
	myStatic real_t gMod[NVMAX];

	/* Perform normal QP solution if QP has not been regularised. */
	if ( QProblem_usingRegularisation( THIS ) == BT_FALSE )
		return QProblem_solveQP( THIS,g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,cputime,nWSRperformed );


	/* I) SOLVE USUAL REGULARISED QP */
	if ( cputime == 0 )
	{
		returnvalue = QProblem_solveQP( THIS,g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,0,nWSRperformed );
	}
	else
	{
		cputime_cur = *cputime;
		returnvalue = QProblem_solveQP( THIS,g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,&cputime_cur,nWSRperformed );
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
			gMod[i] = g_new[i] - THIS->regVal*THIS->x[i];

		/* 2) Solve regularised QP with modified gradient allowing
		 *    only as many working set recalculations and CPU time
		 *    as have been left from previous QP solutions. */
		*nWSR = nWSR_max;

		if ( cputime == 0 )
		{
			returnvalue = QProblem_solveQP( THIS,gMod,lb_new,ub_new,lbA_new,ubA_new, nWSR,0,nWSR_total );
		}
		else
		{
			cputime_cur = *cputime - cputime_total;
			returnvalue = QProblem_solveQP( THIS,gMod,lb_new,ub_new,lbA_new,ubA_new, nWSR,&cputime_cur,nWSR_total );
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
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblem_setupSubjectToType( QProblem* THIS )
{
	return QProblem_setupSubjectToTypeNew( THIS,THIS->lb,THIS->ub,THIS->lbA,THIS->ubA );
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblem_setupSubjectToTypeNew(	QProblem* THIS,
											const real_t* const lb_new, const real_t* const ub_new,
											const real_t* const lbA_new, const real_t* const ubA_new
											)
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );
	

	/* I) SETUP SUBJECTTOTYPE FOR BOUNDS */
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
			if ( ( lb_new[i] < -QPOASES_INFTY + THIS->options.boundTolerance ) && ( ub_new[i] > QPOASES_INFTY - THIS->options.boundTolerance )
					&& (THIS->options.enableFarBounds == BT_FALSE))
			{
				Bounds_setType( &(THIS->bounds),i,ST_UNBOUNDED );
			}
			else
			{
				if ( THIS->options.enableEqualities && lb_new[i] > ub_new[i] - THIS->options.boundTolerance )
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


	/* II) SETUP SUBJECTTOTYPE FOR CONSTRAINTS */
	/* 1) Check if lower constraints' bounds are present. */
	Constraints_setNoLower( &(THIS->constraints),BT_TRUE );
	if ( lbA_new != 0 )
	{
		for( i=0; i<nC; ++i )
		{
			if ( lbA_new[i] > -QPOASES_INFTY )
			{
				Constraints_setNoLower( &(THIS->constraints),BT_FALSE );
				break;
			}
		}
	}

	/* 2) Check if upper constraints' bounds are present. */
	Constraints_setNoUpper( &(THIS->constraints),BT_TRUE );
	if ( ubA_new != 0 )
	{
		for( i=0; i<nC; ++i )
		{
			if ( ubA_new[i] < QPOASES_INFTY )
			{
				Constraints_setNoUpper( &(THIS->constraints),BT_FALSE );
				break;
			}
		}
	}

	/* 3) Determine implicit equality constraints and unbounded constraints. */
	if ( ( lbA_new != 0 ) && ( ubA_new != 0 ) )
	{
		for( i=0; i<nC; ++i )
		{
			if (Constraints_getType(&(THIS->constraints),i) == ST_DISABLED)
				continue;
			
			if ( ( lbA_new[i] <= -QPOASES_INFTY ) && ( ubA_new[i] >= QPOASES_INFTY )
					&& (THIS->options.enableFarBounds == BT_FALSE))
			{
				Constraints_setType( &(THIS->constraints),i,ST_UNBOUNDED );
			}
			else
			{
				if ( THIS->options.enableEqualities && THIS->lbA[i] > THIS->ubA[i] - THIS->options.boundTolerance 
				                              && lbA_new[i] > ubA_new[i] - THIS->options.boundTolerance)
					Constraints_setType( &(THIS->constraints),i,ST_EQUALITY );
				else
					Constraints_setType( &(THIS->constraints),i,ST_BOUNDED );
			}
		}
	}
	else
	{
		if ( ( lbA_new == 0 ) && ( ubA_new == 0 ) )
		{
			for( i=0; i<nC; ++i )
				Constraints_setType( &(THIS->constraints),i,ST_UNBOUNDED );
		}
		else
		{
			for( i=0; i<nC; ++i )
				Constraints_setType( &(THIS->constraints),i,ST_BOUNDED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c h o l e s k y D e c o m p o s i t i o n P r o j e c t e d
 */
returnValue QProblem_computeProjectedCholesky( QProblem* THIS )
{
	int i, j;
	int nV  = QProblem_getNV( THIS );
	int nZ  = QProblem_getNZ( THIS );
	int nFR = QProblem_getNFR( THIS );
	
	int *FR_idx, *AC_idx;
	
	long info = 0;
	unsigned long _nZ = (unsigned long)nZ, _nV = NVMAX;

	/* Revert to unprotected Cholesky decomposition */
	if ( QProblem_getNFX( THIS ) + QProblem_getNAC( THIS ) == 0 )
		return QProblemBCPY_computeCholesky( THIS );
	
	/* 1) Initialises R with all zeros. */
	for( i=0; i<NVMAX*NVMAX; ++i )
		THIS->R[i] = 0.0;

	/* Do not do anything for empty null spaces (important for LP case, HST_ZERO !)*/
	if ( nZ == 0 ) /* nZ == nV - QProblem_getNFX( THIS ) - QProblem_getNAC( THIS ) */
		return SUCCESSFUL_RETURN;
	
	/* 2) Calculate Cholesky decomposition of projected Hessian Z'*H*Z. */
	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );

	/* calculate Z'*H*Z */
	switch ( THIS->hessianType )
	{
		case HST_ZERO:
			if ( QProblem_usingRegularisation( THIS ) == BT_TRUE )
			{
				/*Id = createDiagSparseMat( nV, THIS->regVal );
				DenseMatrix_bilinear( Id,Bounds_getFree( &(THIS->bounds) ), nZ, THIS->Q, NVMAX, THIS->R, NVMAX );*/
				/*fprintf( stderr,"\n\n!!!!!!! NOT YET IMPLEMENTED !!!!!!!!!!!!\n\n" );*/
			}
			else
			{
				return THROWERROR( RET_CHOLESKY_OF_ZERO_HESSIAN );
			}
			break;
		
		case HST_IDENTITY:
			/*Id = createDiagSparseMat( nV, 1.0 );
			DenseMatrix_bilinear( Id,Bounds_getFree( &(THIS->bounds) ), nZ, THIS->Q, NVMAX, THIS->R, NVMAX );*/
			/*fprintf( stderr,"\n\n!!!!!!! NOT YET IMPLEMENTED !!!!!!!!!!!!\n\n" );*/
			break;
	
		default:			
			if ( Indexlist_getLength( Constraints_getActive( &(THIS->constraints)) ) == 0 ) {
				/* make Z trivial */
				for ( j=0; j < nZ; ++j ) {
					for ( i=0; i < nV; ++i )
						QQ(i,j) = 0.0;
					QQ(FR_idx[j],j) = 1.0;
				}
				/* now Z is trivial, and so is Z'HZ */
				for ( j=0; j < nFR; ++j )
					DenseMatrix_getCol(THIS->H,FR_idx[j], Bounds_getFree( &(THIS->bounds) ), 1.0, &(THIS->R[j*NVMAX]));
			} else {
				/* THIS is expensive if Z is large! */
				DenseMatrix_bilinear(THIS->H,Bounds_getFree( &(THIS->bounds) ), nZ, THIS->Q, NVMAX, THIS->R, NVMAX);
			}
	}

	/* R'*R = Z'*H*Z */
	POTRF( "U", &_nZ, THIS->R, &_nV, &info );

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
	for (i=0;i<nZ-1;++i)
		RR(i+1,i) = 0.0;
	
	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p T Q f a c t o r i s a t i o n
 */
returnValue QProblem_setupTQfactorisation( QProblem* THIS )
{
	int i, ii;
	int nFR = QProblem_getNFR( THIS );

	int* FR_idx;
	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );

	/* 1) Set Q to unity matrix. */
	for( i=0; i<NVMAX*NVMAX; ++i )
		THIS->Q[i] = 0.0;

	for( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		QQ(ii,i) = 1.0;
	}

 	/* 2) Set T to zero matrix. */
	for( i=0; i<NVCMIN*NVCMIN; ++i )
		THIS->T[i] = 0.0;

	return SUCCESSFUL_RETURN;
}


/*
 *	o b t a i n A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblem_obtainAuxiliaryWorkingSet(	QProblem* THIS,
												const real_t* const xOpt, const real_t* const yOpt,
												Bounds* const guessedBounds, Constraints* const guessedConstraints,
												Bounds* auxiliaryBounds,Constraints* auxiliaryConstraints
												)
{
	int i = 0;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	SubjectToStatus guessedStatus;

	
	/* 1) Ensure that desiredBounds is allocated (and different from guessedBounds). */
	if ( ( auxiliaryBounds == 0 ) || ( auxiliaryBounds == guessedBounds ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( auxiliaryConstraints == 0 ) || ( auxiliaryConstraints == guessedConstraints ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup working set of bounds for auxiliary initial QP. */
	if ( QProblemBCPY_obtainAuxiliaryWorkingSet( THIS,xOpt,yOpt,guessedBounds, auxiliaryBounds ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
	
	/* 3) Setup working set of constraints for auxiliary initial QP. */
	if ( guessedConstraints != 0 )
	{
		/* If an initial working set is specific, use it!
		 * Moreover, add all equality constraints if specified. */
		for( i=0; i<nC; ++i )
		{
			/* Add constraint only if it is not (goint to be) disabled! */
			guessedStatus = Constraints_getStatus( guessedConstraints,i );

			#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
			if ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY )
			{
				if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_LOWER ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
			else
			#endif
			{
				if ( Constraints_setupConstraint( auxiliaryConstraints,i,guessedStatus ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
		}
	}
	else	/* No initial working set specified. */
	{
		/* Obtain initial working set by "clipping". */
		if ( ( xOpt != 0 ) && ( yOpt == 0 ) )
		{
			for( i=0; i<nC; ++i )
			{
				if ( THIS->Ax[i] - THIS->lbA[i] <= THIS->options.boundTolerance )
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( THIS->ubA[i] - THIS->Ax_u[i] <= THIS->options.boundTolerance )
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all equality constraints if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY )
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		/* Obtain initial working set in accordance to sign of dual solution vector. */
		if ( ( xOpt == 0 ) && ( yOpt != 0 ) )
		{
			for( i=0; i<nC; ++i )
			{
				if ( yOpt[nV+i] > QPOASES_EPS )
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( yOpt[nV+i] < -QPOASES_EPS )
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all equality constraints if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY )
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		/* If xOpt and yOpt are null pointer and no initial working is specified,
		 * start with empty working set (or implicitly fixed bounds and equality constraints only)
		 * for auxiliary QP. */
		if ( ( xOpt == 0 ) && ( yOpt == 0 ) )
		{
			for( i=0; i<nC; ++i )
			{
				/* Only add all equality constraints if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY )
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( Constraints_setupConstraint( auxiliaryConstraints,i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblem_setupAuxiliaryWorkingSet(	QProblem* THIS,
												Bounds* const auxiliaryBounds,
												Constraints* const auxiliaryConstraints,
												BooleanType setupAfresh
												)
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );
	BooleanType WSisTrivial = BT_TRUE;
	
	BooleanType updateCholesky;
	BooleanType was_fulli = THIS->options.enableFullLITests;
	real_t backupEpsLITests = THIS->options.epsLITests;

	/* consistency checks */
	if ( auxiliaryBounds != 0 )
	{
		for( i=0; i<nV; ++i )
			if ( ( Bounds_getStatus( &(THIS->bounds),i ) == ST_UNDEFINED ) || ( Bounds_getStatus( auxiliaryBounds,i ) == ST_UNDEFINED ) )
				return THROWERROR( RET_UNKNOWN_BUG );
	}
	else
	{
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	if ( auxiliaryConstraints != 0 )
	{
		for( i=0; i<nC; ++i )
			if ( ( Constraints_getStatus( &(THIS->constraints),i ) == ST_UNDEFINED ) || ( Constraints_getStatus( auxiliaryConstraints,i ) == ST_UNDEFINED ) )
				return THROWERROR( RET_UNKNOWN_BUG );
	}
	else
	{
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* Check for trivial working set (all and only bounds active) */
	for (i = 0; i < nV; i++)
		if (Bounds_getStatus( auxiliaryBounds,i) == ST_INACTIVE)
		{
			WSisTrivial = BT_FALSE;
			break;
		}
	for (i = 0; i < nC; i++)
/*		(ckirches) here we chose to ignore an invalid ST_INACTIVE on 
		           constraints that are ST_EQUALITies or may just have become equalities*/
		if ( (Constraints_getType( &(THIS->constraints),i) == ST_EQUALITY) /* NOT auxiliaryConstraints here*/
			|| (Constraints_getStatus(auxiliaryConstraints,i) != ST_INACTIVE) )
		{
			WSisTrivial = BT_FALSE;
			break;
		}

	if (WSisTrivial == BT_TRUE)
	{
		for (i = 0; i < nV; i++)
			if (Bounds_getStatus( &(THIS->bounds),i) == ST_INACTIVE)
				Bounds_moveFreeToFixed( &(THIS->bounds),i, Bounds_getStatus( auxiliaryBounds,i));

		return SUCCESSFUL_RETURN;
	}


	/* I) SETUP CHOLESKY FLAG:
	 *    Cholesky decomposition shall only be updated if working set
	 *    shall be updated (i.e. NOT setup afresh!) */
	if ( setupAfresh == BT_TRUE )
		updateCholesky = BT_FALSE;
	else
		updateCholesky = BT_TRUE;

	THIS->options.enableFullLITests = BT_FALSE;
	/* THIS->options.epsLITests = 1e-1; */

	/* II) REMOVE FORMERLY ACTIVE (CONSTRAINTS') BOUNDS (IF NECESSARY): */
	if ( setupAfresh == BT_FALSE )
	{
		/* 1) Remove all active constraints that shall be inactive or disabled AND
		*    all active constraints that are active at the wrong bound. */
		for( i=0; i<nC; ++i )
		{
			if ( ( Constraints_getStatus( &(THIS->constraints),i ) == ST_LOWER ) && ( Constraints_getStatus( auxiliaryConstraints,i ) != ST_LOWER ) )
				if ( QProblem_removeConstraint( THIS,i,updateCholesky,BT_FALSE,THIS->options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

			if ( ( Constraints_getStatus( &(THIS->constraints),i ) == ST_UPPER ) && ( Constraints_getStatus( auxiliaryConstraints,i ) != ST_UPPER ) )
				if ( QProblem_removeConstraint( THIS,i,updateCholesky,BT_FALSE,THIS->options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}

		/* 2) Remove all active bounds that shall be inactive AND
		*    all active bounds that are active at the wrong bound. */
		for( i=0; i<nV; ++i )
		{
			if ( ( Bounds_getStatus( &(THIS->bounds),i ) == ST_LOWER ) && ( Bounds_getStatus( auxiliaryBounds,i ) != ST_LOWER ) )
				if ( QProblem_removeBound( THIS,i,updateCholesky,BT_FALSE,THIS->options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

			if ( ( Bounds_getStatus( &(THIS->bounds),i ) == ST_UPPER ) && ( Bounds_getStatus( auxiliaryBounds,i ) != ST_UPPER ) )
				if ( QProblem_removeBound( THIS,i,updateCholesky,BT_FALSE,THIS->options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}


	/* III) ADD NEWLY ACTIVE (CONSTRAINTS') BOUNDS: */

	/* 1) Add all equality bounds. */
	for( i=0; i<nV; ++i )
	{
		/*if ( ( Bounds_getType( &(THIS->bounds),i ) == ST_EQUALITY ) && ( ( Bounds_getStatus( &(THIS->bounds),i ) == ST_INACTIVE ) && ( Bounds_getStatus( auxiliaryBounds,i ) != ST_INACTIVE ) ) )
		
		(ckirches) force equalities active*/
		
		if ( ( Bounds_getType( &(THIS->bounds),i ) == ST_EQUALITY ) && ( Bounds_getStatus( &(THIS->bounds),i ) == ST_INACTIVE ) )
		{
            /* assert ( Bounds_getStatus( auxiliaryBounds,i ) != ST_INACTIVE ); */
			/* No check for linear independence necessary. */
			if ( QProblem_addBound( THIS,i,ST_LOWER,updateCholesky,BT_TRUE ) != SUCCESSFUL_RETURN ) /* was Bounds_getStatus( auxiliaryBounds,i ) */
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}

	/* 2) Add all equality constraints. */
	for( i=0; i<nC; ++i )
	{
        /*if ( ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY ) && ( ( Constraints_getStatus( &(THIS->constraints),i ) == ST_INACTIVE ) && ( Constraints_getStatus( auxiliaryConstraints,i ) != ST_INACTIVE ) ) )

		(ckirches) force equalities active */
		
		if ( ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY ) && ( Constraints_getStatus( &(THIS->constraints),i ) == ST_INACTIVE ) )
		{
            /* assert ( Constraints_getStatus( auxiliaryConstraints,i ) != ST_INACTIVE ); */
			/* Add constraint only if it is linearly independent from the current working set. */
			if ( QProblem_addConstraint_checkLI( THIS,i ) == RET_LINEARLY_INDEPENDENT )
			{
				if ( QProblem_addConstraint( THIS,i,ST_LOWER,updateCholesky,BT_TRUE ) != SUCCESSFUL_RETURN )  /* was Constraints_getStatus( auxiliaryConstraints,i )*/
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
			}
			else
			{
				/* Equalities are not linearly independent! */
				Constraints_setType( &(THIS->constraints),i, ST_BOUNDED );
			}
		}
	}


	/* 3) Add all inactive bounds that shall be active AND
	 *    all formerly active bounds that have been active at the wrong bound. */
	for( i=0; i<nV; ++i )
	{
		if ( ( Bounds_getType( &(THIS->bounds),i ) != ST_EQUALITY ) && ( ( Bounds_getStatus( &(THIS->bounds),i ) == ST_INACTIVE ) && ( Bounds_getStatus( auxiliaryBounds,i ) != ST_INACTIVE ) ) )
		{
			/* Add bound only if it is linearly independent from the current working set. */
			if ( QProblem_addBound_checkLI( THIS,i ) == RET_LINEARLY_INDEPENDENT )
			{
				if ( QProblem_addBound( THIS,i,Bounds_getStatus( auxiliaryBounds,i ),updateCholesky,BT_TRUE ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
			}
		}
	}

	/* 4) Add all inactive constraints that shall be active AND
	 *    all formerly active constraints that have been active at the wrong bound. */
	for( i=0; i<nC; ++i )
	{
		if ( ( Constraints_getType( &(THIS->constraints),i ) != ST_EQUALITY ) && ( Constraints_getStatus( auxiliaryConstraints,i ) != ST_INACTIVE ) )
		{
			/* formerly inactive */
			if ( Constraints_getStatus( &(THIS->constraints),i ) == ST_INACTIVE )
			{
				/* Add constraint only if it is linearly independent from the current working set. */
				if ( QProblem_addConstraint_checkLI( THIS,i ) == RET_LINEARLY_INDEPENDENT )
				{
					if ( QProblem_addConstraint( THIS,i,Constraints_getStatus( auxiliaryConstraints,i ),updateCholesky,BT_TRUE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
				}
			}
		}
	}

	THIS->options.enableFullLITests = was_fulli;
	THIS->options.epsLITests = backupEpsLITests;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P s o l u t i o n
 */
returnValue QProblem_setupAuxiliaryQPsolution(	QProblem* THIS,
												const real_t* const xOpt, const real_t* const yOpt
												)
{
	int i, j;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );


	/* Setup primal/dual solution vector for auxiliary initial QP:
	 * if a null pointer is passed, a zero vector is assigned;
	 *  old solution vector is kept if pointer to internal solution vevtor is passed. */
	if ( xOpt != 0 )
	{
		if ( xOpt != THIS->x )
			for( i=0; i<nV; ++i )
				THIS->x[i] = xOpt[i];

		DenseMatrix_times(THIS->A,1, 1.0, THIS->x, nV, 0.0, THIS->Ax, nC);

		for ( j=0; j<nC; ++j )
		{
			THIS->Ax_l[j] = THIS->Ax[j];
			THIS->Ax_u[j] = THIS->Ax[j];
		}
	}
	else
	{
		for( i=0; i<nV; ++i )
			THIS->x[i] = 0.0;

		for ( j=0; j<nC; ++j )
		{
			THIS->Ax[j] = 0.0;
			THIS->Ax_l[j] = 0.0;
			THIS->Ax_u[j] = 0.0;
		}
	}

	if ( yOpt != 0 )
	{
		if ( yOpt != THIS->y )
			for( i=0; i<nV+nC; ++i )
				THIS->y[i] = yOpt[i];
	}
	else
	{
		for( i=0; i<nV+nC; ++i )
			THIS->y[i] = 0.0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P g r a d i e n t
 */
returnValue QProblem_setupAuxiliaryQPgradient( QProblem* THIS )
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );


	/* Setup gradient vector: g = -H*x + [Id A]'*[yB yC]
	 *                          = yB - H*x + A'*yC. */
	switch ( THIS->hessianType )
	{
		case HST_ZERO:
			if ( QProblem_usingRegularisation( THIS ) == BT_FALSE )
				for ( i=0; i<nV; ++i )
					THIS->g[i] = THIS->y[i];
			else
				for ( i=0; i<nV; ++i )
					THIS->g[i] = THIS->y[i] - THIS->regVal*THIS->x[i];
			break;

		case HST_IDENTITY:
			for ( i=0; i<nV; ++i )
				THIS->g[i] = THIS->y[i] - THIS->x[i];
			break;

		default:
			/* y'*Id */
			for ( i=0; i<nV; ++i )
				THIS->g[i] = THIS->y[i];

			/* - H*x */
			DenseMatrix_times(THIS->H,1, -1.0, THIS->x, nV, 1.0, THIS->g, nV);
			break;
	}

	/* + A'*yC */
	DenseMatrix_transTimes(THIS->A,1, 1.0, &(THIS->y[nV]), nC, 1.0, THIS->g, nV);
	
	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P b o u n d s
 */
returnValue QProblem_setupAuxiliaryQPbounds(	QProblem* THIS,
												Bounds* const auxiliaryBounds,
												Constraints* const auxiliaryConstraints,
												BooleanType useRelaxation
												)
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );


	/* 1) Setup bound vectors. */
	for ( i=0; i<nV; ++i )
	{
		switch ( Bounds_getStatus( &(THIS->bounds),i ) )
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
						/* If a bound is inactive although it was supposed to be
						* active by the auxiliaryBounds it could not be added
						* due to linear dependence. Thus set it "strongly inactive". */
						if ( Bounds_getStatus( auxiliaryBounds,i ) == ST_LOWER )
							THIS->lb[i] = THIS->x[i];
						else
							THIS->lb[i] = THIS->x[i] - THIS->options.boundRelaxation;

						if ( Bounds_getStatus( auxiliaryBounds,i ) == ST_UPPER )
							THIS->ub[i] = THIS->x[i];
						else
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

	/* 2) Setup constraints vectors. */
	for ( i=0; i<nC; ++i )
	{
		switch ( Constraints_getStatus( &(THIS->constraints),i ) )
		{
			case ST_INACTIVE:
				if ( useRelaxation == BT_TRUE )
				{
					if ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY )
					{
						THIS->lbA[i] = THIS->Ax_l[i];
						THIS->ubA[i] = THIS->Ax_u[i];
					}
					else
					{
						/* If a constraint is inactive although it was supposed to be
						* active by the auxiliaryConstraints, it could not be added
						* due to linear dependence. Thus set it "strongly inactive". */
						if ( Constraints_getStatus( auxiliaryConstraints,i ) == ST_LOWER )
							THIS->lbA[i] = THIS->Ax_l[i];
						else
							THIS->lbA[i] = THIS->Ax_l[i] - THIS->options.boundRelaxation;

						if ( Constraints_getStatus( auxiliaryConstraints,i ) == ST_UPPER )
							THIS->ubA[i] = THIS->Ax_u[i];
						else
							THIS->ubA[i] = THIS->Ax_u[i] + THIS->options.boundRelaxation;
					}
				}
				break;

			case ST_LOWER:
				THIS->lbA[i] = THIS->Ax_l[i];
				if ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY )
				{
					THIS->ubA[i] = THIS->Ax_l[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						THIS->ubA[i] = THIS->Ax_l[i] + THIS->options.boundRelaxation;
				}
				break;

			case ST_UPPER:
				THIS->ubA[i] = THIS->Ax_u[i];
				if ( Constraints_getType( &(THIS->constraints),i ) == ST_EQUALITY )
				{
					THIS->lbA[i] = THIS->Ax_u[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						THIS->lbA[i] = THIS->Ax_u[i] - THIS->options.boundRelaxation;
				}
				break;

            case ST_DISABLED:
                break;
                
			default:
				return THROWERROR( RET_UNKNOWN_BUG );
		}
		THIS->Ax_l[i] = THIS->Ax_l[i] - THIS->lbA[i];
		THIS->Ax_u[i] = THIS->ubA[i] - THIS->Ax_u[i];
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	a d d C o n s t r a i n t
 */
returnValue QProblem_addConstraint(	QProblem* THIS,
									int number, SubjectToStatus C_status,
									BooleanType updateCholesky,
									BooleanType ensureLI
									)
{
	int i, j, ii;

	returnValue ensureLIreturnvalue;
	
	int nFR, nAC, nZ, tcol;
	int* FR_idx;
	
	myStatic real_t aFR[NVMAX];
	myStatic real_t wZ[NVMAX];
	
	real_t c, s, nu;


	/* consistency checks */
	if ( Constraints_getStatus( &(THIS->constraints),number ) != ST_INACTIVE )
		return THROWERROR( RET_CONSTRAINT_ALREADY_ACTIVE );

	if ( ( Constraints_getNC( &(THIS->constraints) ) - QProblem_getNAC( THIS ) ) == Constraints_getNUC( &(THIS->constraints) ) )
		return THROWERROR( RET_ALL_CONSTRAINTS_ACTIVE );

	if ( ( QProblem_getStatus( THIS ) == QPS_NOTINITIALISED )    ||
		 ( QProblem_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblem_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( QProblem_getStatus( THIS ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}


	/* I) ENSURE LINEAR INDEPENDENCE OF THE WORKING SET,
	 *    i.e. remove a constraint or bound if linear dependence occurs. */
	/* check for LI only if Cholesky decomposition shall be updated! */
	if ( updateCholesky == BT_TRUE && ensureLI == BT_TRUE )
	{
		ensureLIreturnvalue = QProblem_addConstraint_ensureLI( THIS,number,C_status );

		switch ( ensureLIreturnvalue )
		{
			case SUCCESSFUL_RETURN:
				break;

			case RET_LI_RESOLVED:
				break;

			case RET_ENSURELI_FAILED_NOINDEX:
				return RET_ADDCONSTRAINT_FAILED_INFEASIBILITY;

			case RET_ENSURELI_FAILED_CYCLING:
				return RET_ADDCONSTRAINT_FAILED_INFEASIBILITY;
				
			case RET_ENSURELI_DROPPED:
				return SUCCESSFUL_RETURN;

			default:
				return THROWERROR( RET_ENSURELI_FAILED );
		}
	}

	/* some definitions */
	nFR = QProblem_getNFR( THIS );
	nAC = QProblem_getNAC( THIS );
	nZ  = QProblem_getNZ( THIS );

	tcol = THIS->sizeT - nAC;

	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );

	for( i=0; i<nZ; ++i )
		wZ[i] = 0.0;


	/* II) ADD NEW ACTIVE CONSTRAINT TO MATRIX T: */
	/* 1) Add row [wZ wY] = aFR'*[Z Y] to the end of T: assign aFR. */
	DenseMatrix_getRow(THIS->A,number, Bounds_getFree( &(THIS->bounds) ), 1.0, aFR);

	/* calculate wZ */
	for( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		for( j=0; j<nZ; ++j )
			wZ[j] += aFR[i] * QQ(ii,j);
	}

	/* 2) Calculate wY and store it directly into T. */
	if ( nAC > 0 )
	{
		for( j=0; j<nAC; ++j )
			TT(nAC,tcol+j) = 0.0;
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			for( j=0; j<nAC; ++j )
				TT(nAC,tcol+j) += aFR[i] * QQ(ii,nZ+j);
		}
	}

	if ( nZ > 0 )
	{
		/* II) RESTORE TRIANGULAR FORM OF T: */
		/*     Use column-wise Givens rotations to restore reverse triangular form
		*      of T, simultanenous change of Q (i.e. Z) and R. */
		for( j=0; j<nZ-1; ++j )
		{
			QProblemB_computeGivens( wZ[j+1],wZ[j], &wZ[j+1],&wZ[j],&c,&s );
			nu = s/(1.0+c);

			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				QProblemB_applyGivens( c,s,nu,QQ(ii,1+j),QQ(ii,j), &QQ(ii,1+j),&QQ(ii,j) );
			}

			if ( ( updateCholesky == BT_TRUE ) &&
				 ( THIS->hessianType != HST_ZERO )   && ( THIS->hessianType != HST_IDENTITY ) )
			{
				for( i=0; i<=j+1; ++i )
					QProblemB_applyGivens( c,s,nu,RR(i,1+j),RR(i,j), &RR(i,1+j),&RR(i,j) );
			}
		}

		TT(nAC,tcol-1) = wZ[nZ-1];


		if ( ( updateCholesky == BT_TRUE ) &&
			 ( THIS->hessianType != HST_ZERO )   && ( THIS->hessianType != HST_IDENTITY ) )
		{
			/* III) RESTORE TRIANGULAR FORM OF R:
			 *      Use row-wise Givens rotations to restore upper triangular form of R. */
			for( i=0; i<nZ-1; ++i )
			{
				QProblemB_computeGivens( RR(i,i),RR(1+i,i), &RR(i,i),&RR(1+i,i),&c,&s );
				nu = s/(1.0+c);

				for( j=(1+i); j<(nZ-1); ++j ) /* last column of R is thrown away */
					QProblemB_applyGivens( c,s,nu,RR(i,j),RR(1+i,j), &RR(i,j),&RR(1+i,j) );
			}
			/* last column of R is thrown away */
			for( i=0; i<nZ; ++i )
				RR(i,nZ-1) = 0.0;
		}
	}


	/* IV) UPDATE INDICES */
	THIS->tabularOutput.idxAddC = number;

	if ( Constraints_moveInactiveToActive( &(THIS->constraints),number,C_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDCONSTRAINT_FAILED );


	return SUCCESSFUL_RETURN;
}



/*
 *	a d d C o n s t r a i n t _ c h e c k L I
 */
returnValue QProblem_addConstraint_checkLI( QProblem* THIS, int number )
{
	returnValue returnvalue = RET_LINEARLY_DEPENDENT;

	int i, j, ii;
	int nV  = QProblem_getNV( THIS );
	int nFR = QProblem_getNFR( THIS );
	int nZ  = QProblem_getNZ( THIS );
	int nC  = QProblem_getNC( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nFX = QProblem_getNFX( THIS );
	int dim;
	
	int *FX_idx, *AC_idx, *IAC_idx, *FR_idx;

	myStatic real_t delta_g[NVMAX];
	myStatic real_t delta_xFX[NVMAX];
	myStatic real_t delta_xFR[NVMAX];
	myStatic real_t delta_yAC[NCMAX];
	myStatic real_t delta_yFX[NVMAX];
	
	myStatic real_t nul[NVCMAX];
	myStatic real_t Arow[NVMAX];
	
	real_t weight = 0.0;
	real_t zero = 0.0;
	
	real_t sum, l2;
	
	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );


	if (THIS->options.enableFullLITests)
	{
		/*
		 * expensive LI test. Backsolve with refinement using special right
		 * hand side. This gives an estimate for what should be considered
		 * "zero". We then check linear independence relative to THIS estimate.
		 */
		Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
		Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );
		Indexlist_getNumberArray( Constraints_getInactive( &(THIS->constraints) ),&IAC_idx );

		dim = (nC>nV)?nC:nV;
		for (ii = 0; ii < dim; ++ii)
			nul[ii]=0.0;

		DenseMatrix_getRow(THIS->A,number, 0, 1.0, delta_g);

		returnvalue = QProblem_determineStepDirection(	THIS,delta_g,
														nul, nul, nul, nul,
														BT_FALSE, BT_FALSE,
														delta_xFX, delta_xFR, delta_yAC, delta_yFX);

		/* compute the weight in inf-norm */
		for (ii = 0; ii < nAC; ++ii)
		{
			real_t a = qpOASES_getAbs (delta_yAC[ii]);
			if (weight < a) weight = a;
		}
		for (ii = 0; ii < nFX; ++ii)
		{
			real_t a = qpOASES_getAbs (delta_yFX[ii]);
			if (weight < a) weight = a;
		}

		/* look at the "zero" in a relative inf-norm */
		for (ii = 0; ii < nFX; ++ii)
		{
			real_t a = qpOASES_getAbs (delta_xFX[ii]);
			if (zero < a) zero = a;
		}
		for (ii = 0; ii < nFR; ++ii)
		{
			real_t a = qpOASES_getAbs (delta_xFR[ii]);
			if (zero < a) zero = a;
		}

		/* relative test against zero in inf-norm */
		if (zero > THIS->options.epsLITests * weight)
			returnvalue = RET_LINEARLY_INDEPENDENT;

	}
	else
	{
		/*
		 * cheap LI test for constraint. Check if constraint <number> is
		 * linearly independent from the the active ones (<=> is element of null
		 * space of Afr).
		 */

		DenseMatrix_getRow(THIS->A,number, Bounds_getFree( &(THIS->bounds) ), 1.0, Arow);

		l2  = 0.0;
		for (i = 0; i < nFR; i++)
			l2  += Arow[i]*Arow[i];

		for( j=0; j<nZ; ++j )
		{
			sum = 0.0;
			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				sum += Arow[i] * QQ(ii,j);
			}

			if ( qpOASES_getAbs( sum ) > THIS->options.epsLITests*l2 )
			{
				/*fprintf(stdFile, "LI test: |sum| = %9.2e, l2 = %9.2e, var = %d\n", qpOASES_getAbs(sum), l2, jj+1); */
				returnvalue = RET_LINEARLY_INDEPENDENT;
				break;
			}
		}
	}

	return THROWINFO( returnvalue );
}


/*
 *	a d d C o n s t r a i n t _ e n s u r e L I
 */
returnValue QProblem_addConstraint_ensureLI( QProblem* THIS, int number, SubjectToStatus C_status )
{
	int i, j, ii, jj;
	#ifndef __XPCTARGET__
	myStatic char messageString[QPOASES_MAX_STRING_LENGTH];
	#endif

	int nV  = QProblem_getNV( THIS );
	int nFR = QProblem_getNFR( THIS );
	int nFX = QProblem_getNFX( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nZ  = QProblem_getNZ( THIS );

	int *FR_idx, *FX_idx, *AC_idx;
	
	myStatic real_t xiC[NCMAX];
	myStatic real_t xiC_TMP[NCMAX];
	myStatic real_t xiB[NVMAX];
	myStatic real_t Arow[NVMAX];
	myStatic real_t num[NVMAX];
	
	returnValue returnvalue = SUCCESSFUL_RETURN;

	real_t y_min = THIS->options.maxDualJump;
	int y_min_number = -1;
	int y_min_number_bound = -1;
	BooleanType y_min_isBound = BT_FALSE;


	/* I) Check if new constraint is linearly independent from the active ones. */
	returnValue returnvalueCheckLI = QProblem_addConstraint_checkLI( THIS,number );

	if ( returnvalueCheckLI == RET_INDEXLIST_CORRUPTED )
		return THROWERROR( RET_ENSURELI_FAILED );

	if ( returnvalueCheckLI == RET_LINEARLY_INDEPENDENT )
		return SUCCESSFUL_RETURN;


 	/* II) NEW CONSTRAINT IS LINEARLY DEPENDENT: */
	/* 1) Determine coefficients of linear combination,
	 *    cf. M.J. Best. Applied Mathematics and Parallel Computing, chapter:
	 *    An Algorithm for the Solution of the Parametric Quadratic Programming
	 *    Problem, pages 57-76. Physica-Verlag, Heidelberg, 1996. */
	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );

	DenseMatrix_getRow(THIS->A,number, Bounds_getFree( &(THIS->bounds) ), C_status == ST_LOWER ? 1.0 : -1.0, Arow);

	/* 2) Calculate xiC */
	if ( nAC > 0 )
	{
		for( i=0; i<nAC; ++i )
		{
			xiC_TMP[i] = 0.0;
			for( j=0; j<nFR; ++j )
			{
				jj = FR_idx[j];
				xiC_TMP[i] += QQ(jj,nZ+i) * Arow[j];
			}
		}

		if ( QProblem_backsolveT( THIS,xiC_TMP, BT_TRUE, xiC ) != SUCCESSFUL_RETURN )
		{
			returnvalue = RET_ENSURELI_FAILED_TQ;
			goto farewell;
		}
	}

	/* 3) Calculate xiB. */
	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );

	DenseMatrix_getRow(THIS->A,number, Bounds_getFixed( &(THIS->bounds) ), C_status == ST_LOWER ? 1.0 : -1.0, xiB);
	DenseMatrix_subTransTimes(THIS->A,Constraints_getActive( &(THIS->constraints)), Bounds_getFixed( &(THIS->bounds) ), 1, -1.0, xiC, nAC, 1.0, xiB, nFX);

	/* III) DETERMINE CONSTRAINT/BOUND TO BE REMOVED. */
	
	/* 1) Constraints. */
	for( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];
		num[i] = THIS->y[nV+ii];
	}

	QProblem_performRatioTestC( THIS,nAC, AC_idx, &(THIS->constraints), num, xiC, THIS->options.epsNum, THIS->options.epsDen, &y_min,&y_min_number);	

	/* 2) Bounds. */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = THIS->y[ii];
	}

	QProblem_performRatioTestB( THIS,nFX, FX_idx, &(THIS->bounds),num, xiB, THIS->options.epsNum, THIS->options.epsDen, &y_min,&y_min_number_bound);

	if ( y_min_number_bound >= 0 )
	{
		y_min_number = y_min_number_bound;
		y_min_isBound = BT_TRUE;
	}


	/* IV) REMOVE CONSTRAINT/BOUND FOR RESOLVING LINEAR DEPENDENCE: */
	if ( y_min_number >= 0 )
	{
		/* Update Lagrange multiplier... */
		for( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			THIS->y[nV+ii] -= y_min * xiC[i];
		}
		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];
			THIS->y[ii] -= y_min * xiB[i];
		}

		/* ... also for newly active constraint... */
		if ( C_status == ST_LOWER )
			THIS->y[nV+number] = y_min;
		else
			THIS->y[nV+number] = -y_min;

		/* ... and for constraint to be removed. */
		if ( y_min_isBound == BT_TRUE )
		{
			#ifndef __XPCTARGET__
			snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"bound no. %d.",y_min_number );
			MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( QProblem_removeBound( THIS,y_min_number,BT_TRUE,BT_FALSE,BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				returnvalue = RET_REMOVE_FROM_ACTIVESET_FAILED;
				goto farewell;
			}
			THIS->tabularOutput.excRemB = 1;
			
			THIS->y[y_min_number] = 0.0;
		}
		else
		{
			#ifndef __XPCTARGET__
			snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"constraint no. %d.",y_min_number );
			MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( QProblem_removeConstraint( THIS,y_min_number,BT_TRUE,BT_FALSE,BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				returnvalue = RET_REMOVE_FROM_ACTIVESET_FAILED;
				goto farewell;
			}
			THIS->tabularOutput.excRemC = 1;

			THIS->y[nV+y_min_number] = 0.0;
		}
	}
	else
	{
		if (THIS->options.enableDropInfeasibles == BT_TRUE) {
			/* dropping of infeasible constraints according to drop priorities */
			returnvalue = QProblem_dropInfeasibles( THIS,number, C_status, BT_FALSE, xiB, xiC);
		}
		else
		{
			/* no constraint/bound can be removed => QP is infeasible! */		
			returnvalue = RET_ENSURELI_FAILED_NOINDEX;
			QProblem_setInfeasibilityFlag( THIS,returnvalue,BT_FALSE );
		}
	}

farewell:
	MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_LI_RESOLVED,0,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
	
	return ( (returnvalue != SUCCESSFUL_RETURN) && (returnvalue != RET_ENSURELI_FAILED_NOINDEX ) ) ? THROWERROR (returnvalue) : returnvalue;
}



/*
 *	a d d B o u n d
 */
returnValue QProblem_addBound(	QProblem* THIS, int number, SubjectToStatus B_status,
								BooleanType updateCholesky,
								BooleanType ensureLI
								)
{
	int i, j, ii;
	returnValue ensureLIreturnvalue;

	int nFR, nAC, nZ, tcol, lastfreenumber;
	
	int* FR_idx;
	myStatic real_t w[NVMAX];
	real_t c, s, nu;
	myStatic real_t tmp[NCMAX];


	/* consistency checks */
	if ( Bounds_getStatus( &(THIS->bounds),number ) != ST_INACTIVE )
		return THROWERROR( RET_BOUND_ALREADY_ACTIVE );

	if ( QProblem_getNFR( THIS ) == Bounds_getNUV( &(THIS->bounds) ) )
		return THROWERROR( RET_ALL_BOUNDS_ACTIVE );

	if ( ( QProblem_getStatus( THIS ) == QPS_NOTINITIALISED )    ||
		 ( QProblem_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblem_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
 		 ( QProblem_getStatus( THIS ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}


	/* I) ENSURE LINEAR INDEPENDENCE OF THE WORKING SET,
	 *    i.e. remove a constraint or bound if linear dependence occurs. */
	/* check for LI only if Cholesky decomposition shall be updated! */
	if ( updateCholesky == BT_TRUE && ensureLI == BT_TRUE )
	{
		ensureLIreturnvalue = QProblem_addBound_ensureLI( THIS,number,B_status );

		switch ( ensureLIreturnvalue )
		{
			case SUCCESSFUL_RETURN:
				break;

			case RET_LI_RESOLVED:
				break;

			case RET_ENSURELI_FAILED_NOINDEX:
				return RET_ADDBOUND_FAILED_INFEASIBILITY;

			case RET_ENSURELI_FAILED_CYCLING:
				return RET_ADDBOUND_FAILED_INFEASIBILITY;

			case RET_ENSURELI_DROPPED:
				return SUCCESSFUL_RETURN;
				
			default:
				return THROWERROR( RET_ENSURELI_FAILED );
		}
	}

	/* some definitions */
	nFR = QProblem_getNFR( THIS );
	nAC = QProblem_getNAC( THIS );
	nZ  = QProblem_getNZ( THIS );

	tcol = THIS->sizeT - nAC;


	/* II) SWAP INDEXLIST OF FREE VARIABLES:
	 *     move the variable to be fixed to the end of the list of free variables. */
	lastfreenumber = Indexlist_getLastNumber( Bounds_getFree( &(THIS->bounds) ) );
	if ( lastfreenumber != number )
		if ( Bounds_swapFree( &(THIS->bounds),number,lastfreenumber ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_ADDBOUND_FAILED );

	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );


	/* III) ADD NEW ACTIVE BOUND TO TOP OF MATRIX T: */
	/* 1) add row [wZ wY] = [Z Y](number) at the top of T: assign w */
	for( i=0; i<nFR; ++i )
		w[i] = QQ(FR_idx[nFR-1],i);


	/* 2) Use column-wise Givens rotations to restore reverse triangular form
	 *    of the first row of T, simultanenous change of Q (i.e. Z) and R. */
	for( j=0; j<nZ-1; ++j )
	{
		QProblemB_computeGivens( w[j+1],w[j], &w[j+1],&w[j],&c,&s );
		nu = s/(1.0+c);

		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			QProblemB_applyGivens( c,s,nu,QQ(ii,1+j),QQ(ii,j), &QQ(ii,1+j),&QQ(ii,j) );
		}

		if ( ( updateCholesky == BT_TRUE ) &&
			 ( THIS->hessianType != HST_ZERO ) && ( THIS->hessianType != HST_IDENTITY ) )
		{
			for( i=0; i<=j+1; ++i )
				QProblemB_applyGivens( c,s,nu,RR(i,1+j),RR(i,j), &RR(i,1+j),&RR(i,j) );
		}
	}


	if ( nAC > 0 )	  /* ( nAC == 0 ) <=> ( nZ == nFR ) <=> Y and T are empty => nothing to do */
	{
		/* store new column a in a temporary vector instead of shifting T one column to the left */
		for( i=0; i<nAC; ++i )
			tmp[i] = 0.0;

		{
			j = nZ-1;

			QProblemB_computeGivens( w[j+1],w[j], &w[j+1],&w[j],&c,&s );
			nu = s/(1.0+c);

			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				QProblemB_applyGivens( c,s,nu,QQ(ii,1+j),QQ(ii,j), &QQ(ii,1+j),&QQ(ii,j) );
			}

			QProblemB_applyGivens( c,s,nu,TT(nAC-1,tcol),tmp[nAC-1], &tmp[nAC-1],&TT(nAC-1,tcol) );
		}

		for( j=nZ; j<nFR-1; ++j )
		{
			QProblemB_computeGivens( w[j+1],w[j], &w[j+1],&w[j],&c,&s );
			nu = s/(1.0+c);

			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				QProblemB_applyGivens( c,s,nu,QQ(ii,1+j),QQ(ii,j), &QQ(ii,1+j),&QQ(ii,j) );
			}

			for( i=(nFR-2-j); i<nAC; ++i )
				QProblemB_applyGivens( c,s,nu,TT(i,1+tcol-nZ+j),tmp[i], &tmp[i],&TT(i,1+tcol-nZ+j) );
		}
	}


	if ( ( updateCholesky == BT_TRUE ) &&
		 ( THIS->hessianType != HST_ZERO ) && ( THIS->hessianType != HST_IDENTITY ) )
	{
		/* IV) RESTORE TRIANGULAR FORM OF R:
		 *     use row-wise Givens rotations to restore upper triangular form of R */
		for( i=0; i<nZ-1; ++i )
		{
			QProblemB_computeGivens( RR(i,i),RR(1+i,i), &RR(i,i),&RR(1+i,i),&c,&s );
			nu = s/(1.0+c);

			for( j=(1+i); j<nZ-1; ++j ) /* last column of R is thrown away */
				QProblemB_applyGivens( c,s,nu,RR(i,j),RR(1+i,j), &RR(i,j),&RR(1+i,j) );
		}
		/* last column of R is thrown away */
		for( i=0; i<nZ; ++i )
			RR(i,nZ-1) = 0.0;
	}


	/* V) UPDATE INDICES */
	THIS->tabularOutput.idxAddB = number;
	if ( Bounds_moveFreeToFixed( &(THIS->bounds),number,B_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDBOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	a d d B o u n d _ c h e c k L I
 */
returnValue QProblem_addBound_checkLI( QProblem* THIS, int number )
{
	int i, ii;
	int nV  = QProblem_getNV( THIS );  /* for QQ() macro */
	int nFR = QProblem_getNFR( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nFX = QProblem_getNFX( THIS );
	int nC  = QProblem_getNC( THIS );
	returnValue returnvalue = RET_LINEARLY_DEPENDENT;
	
	myStatic real_t delta_g[NVMAX];
	myStatic real_t delta_xFX[NVMAX];
	myStatic real_t delta_xFR[NVMAX];
	myStatic real_t delta_yAC[NCMAX];
	myStatic real_t delta_yFX[NVMAX];
	
	int dim, nZ;
	myStatic real_t nul[NVCMAX];
	returnValue dsdReturnValue;
	
	real_t weight = 0.0;
	real_t zero = 0.0;

	if (THIS->options.enableFullLITests)
	{
		/*
		 * expensive LI test. Backsolve with refinement using special right
		 * hand side. This gives an estimate for what should be considered
		 * "zero". We then check linear independence relative to THIS estimate.
		 */

		/*
		 * expensive LI test. Backsolve with refinement using special right
		 * hand side. This gives an estimate for what should be considered
		 * "zero". We then check linear independence relative to THIS estimate.
		 */

		for (ii = 0; ii < nV; ++ii)
			delta_g[ii] = 0.0;
		delta_g[number] = 1.0;	/* sign doesn't matter here */

		dim = (nC>nV)?nC:nV;
		for (ii = 0; ii < dim; ++ii)
			nul[ii]=0.0;

		dsdReturnValue = QProblem_determineStepDirection( THIS,
				delta_g, nul, nul, nul, nul, BT_FALSE, BT_FALSE,
				delta_xFX, delta_xFR, delta_yAC, delta_yFX);
		if (dsdReturnValue != SUCCESSFUL_RETURN) 
			returnvalue = dsdReturnValue;

		/* compute the weight in inf-norm */
		for (ii = 0; ii < nAC; ++ii)
		{
			real_t a = qpOASES_getAbs (delta_yAC[ii]);
			if (weight < a) weight = a;
		}
		for (ii = 0; ii < nFX; ++ii)
		{
			real_t a = qpOASES_getAbs (delta_yFX[ii]);
			if (weight < a) weight = a;
		}

		/* look at the "zero" in a relative inf-norm */
		for (ii = 0; ii < nFX; ++ii)
		{
			real_t a = qpOASES_getAbs (delta_xFX[ii]);
			if (zero < a) zero = a;
		}
		for (ii = 0; ii < nFR; ++ii)
		{
			real_t a = qpOASES_getAbs (delta_xFR[ii]);
			if (zero < a) zero = a;
		}

		/* relative test against zero in inf-norm */
		if (zero > THIS->options.epsLITests * weight)
			returnvalue = RET_LINEARLY_INDEPENDENT;
	}
	else
	{
		/*
		 * cheap LI test for simple bound. Check if constraint <number> is
		 * linearly independent from the the active ones (<=> is element of null
		 * space of Afr).
		 */

		/* some definitions */
		nZ  = QProblem_getNZ( THIS );

		for( i=0; i<nZ; ++i )
			if ( qpOASES_getAbs( QQ(number,i) ) > THIS->options.epsLITests )
			{
				returnvalue = RET_LINEARLY_INDEPENDENT;
				break;
			}
	}

	return THROWINFO( returnvalue );
}


/*
 *	a d d B o u n d _ e n s u r e L I
 */
returnValue QProblem_addBound_ensureLI( QProblem* THIS, int number, SubjectToStatus B_status )
{
	int i, ii;
	#ifndef __XPCTARGET__
	myStatic char messageString[QPOASES_MAX_STRING_LENGTH];
	#endif

	int nV  = QProblem_getNV( THIS );
	int nFX = QProblem_getNFX( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nZ  = QProblem_getNZ( THIS );

	int *FR_idx, *FX_idx, *AC_idx;
	
	myStatic real_t xiC[NCMAX];
	myStatic real_t xiC_TMP[NCMAX];
	myStatic real_t xiB[NVMAX];
	myStatic real_t num[NVMAX];

	real_t y_min = THIS->options.maxDualJump;
	int y_min_number = -1;
	int y_min_number_bound = -1;
	BooleanType y_min_isBound = BT_FALSE;
	
	returnValue returnvalue = SUCCESSFUL_RETURN;


	/* I) Check if new constraint is linearly independent from the active ones. */
	returnValue returnvalueCheckLI = QProblem_addBound_checkLI( THIS,number );

	if ( returnvalueCheckLI == RET_INDEXLIST_CORRUPTED )
		return THROWERROR( RET_ENSURELI_FAILED );

	if ( returnvalueCheckLI == RET_LINEARLY_INDEPENDENT )
		return SUCCESSFUL_RETURN;


 	/* II) NEW BOUND IS LINEARLY DEPENDENT: */
	/* 1) Determine coefficients of linear combination,
	 *    cf. M.J. Best. Applied Mathematics and Parallel Computing, chapter:
	 *    An Algorithm for the Solution of the Parametric Quadratic Programming
	 *    Problem, pages 57-76. Physica-Verlag, Heidelberg, 1996. */
	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );

	/* 2) Calculate xiC. */
	if ( nAC > 0 )
	{
		if ( B_status == ST_LOWER )
		{
			for( i=0; i<nAC; ++i )
				xiC_TMP[i] = QQ(number,nZ+i);
		}
		else
		{
			for( i=0; i<nAC; ++i )
				xiC_TMP[i] = -QQ(number,nZ+i);
		}

		if ( QProblem_backsolveT( THIS,xiC_TMP, BT_TRUE, xiC ) != SUCCESSFUL_RETURN )
		{
			returnvalue = RET_ENSURELI_FAILED_TQ;
			goto farewell;
		}
	}

	/* 3) Calculate xiB. */
	DenseMatrix_subTransTimes(THIS->A,Constraints_getActive( &(THIS->constraints)), Bounds_getFixed( &(THIS->bounds) ), 1, -1.0, xiC, nAC, 0.0, xiB, nFX);


	/* III) DETERMINE CONSTRAINT/BOUND TO BE REMOVED. */

	/* 1) Constraints. */
	for( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];
		num[i] = THIS->y[nV+ii];
	}
	
	QProblem_performRatioTestC( THIS,nAC,AC_idx,&(THIS->constraints), num,xiC, THIS->options.epsNum,THIS->options.epsDen, &y_min,&y_min_number );	

	/* 2) Bounds. */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = THIS->y[ii];
	}

	QProblem_performRatioTestB( THIS,nFX,FX_idx,&(THIS->bounds), num,xiB, THIS->options.epsNum,THIS->options.epsDen, &y_min,&y_min_number_bound );

	if ( y_min_number_bound >= 0 )
	{
		y_min_number = y_min_number_bound;
		y_min_isBound = BT_TRUE;
	}

	/* IV) REMOVE CONSTRAINT/BOUND FOR RESOLVING LINEAR DEPENDENCE: */
	if ( y_min_number >= 0 )
	{
		/* Update Lagrange multiplier... */
		for( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			THIS->y[nV+ii] -= y_min * xiC[i];
		}
		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];
			THIS->y[ii] -= y_min * xiB[i];
		}

		/* ... also for newly active bound ... */
		if ( B_status == ST_LOWER )
			THIS->y[number] = y_min;
		else
			THIS->y[number] = -y_min;

		/* ... and for bound to be removed. */
		if ( y_min_isBound == BT_TRUE )
		{
			#ifndef __XPCTARGET__
			snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"bound no. %d.",y_min_number );
			MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( QProblem_removeBound( THIS,y_min_number,BT_TRUE,BT_FALSE,BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				returnvalue = RET_REMOVE_FROM_ACTIVESET_FAILED;
				goto farewell;
			}
			THIS->tabularOutput.excRemB = 1;

			THIS->y[y_min_number] = 0.0;
		}
		else
		{
			#ifndef __XPCTARGET__
			snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"constraint no. %d.",y_min_number );
			MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( QProblem_removeConstraint( THIS,y_min_number,BT_TRUE,BT_FALSE,BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				returnvalue = RET_REMOVE_FROM_ACTIVESET_FAILED;
				goto farewell;
			}
			THIS->tabularOutput.excRemC = 1;

			THIS->y[nV+y_min_number] = 0.0;
		}
	}
	else
	{
		if (THIS->options.enableDropInfeasibles == BT_TRUE) {
			/* dropping of infeasible constraints according to drop priorities */
			returnvalue = QProblem_dropInfeasibles( THIS,number, B_status, BT_TRUE, xiB, xiC);
		}
		else
		{
			/* no constraint/bound can be removed => QP is infeasible! */
			returnvalue = RET_ENSURELI_FAILED_NOINDEX;
			QProblem_setInfeasibilityFlag( THIS,returnvalue,BT_FALSE );
		}
	}

farewell:
	MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_LI_RESOLVED,0,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
	
	return ( (returnvalue != SUCCESSFUL_RETURN) && (returnvalue != RET_ENSURELI_FAILED_NOINDEX ) ) ? THROWERROR (returnvalue) : returnvalue;
}



/*
 *	r e m o v e C o n s t r a i n t
 */
returnValue QProblem_removeConstraint(	QProblem* THIS, int number,
										BooleanType updateCholesky,
										BooleanType allowFlipping,
										BooleanType ensureNZC
										)
{
	int i, j, ii, jj;
	returnValue returnvalue = SUCCESSFUL_RETURN;
	BooleanType hasFlipped = BT_FALSE;

	/* some definitions */
	int nFR = QProblem_getNFR( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nZ  = QProblem_getNZ( THIS );

	int tcol = THIS->sizeT - nAC;
	int number_idx = Indexlist_getIndex( Constraints_getActive( &(THIS->constraints) ),number );

	int addIdx;
	BooleanType addBoundNotConstraint;
	SubjectToStatus addStatus;
	BooleanType exchangeHappened = BT_FALSE;

	int *FR_idx;
	
	myStatic real_t Hz[NVMAX];
	myStatic real_t z[NVMAX];
	real_t rho2 = 0.0;
	myStatic real_t ZHz[NVMAX];
	myStatic real_t r[NVMAX];
	
	real_t c, s, nu;

	
	/* consistency check */
	if ( ( QProblem_getStatus( THIS ) == QPS_NOTINITIALISED )    ||
		 ( QProblem_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblem_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
 		 ( QProblem_getStatus( THIS ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* consistency checks */
	if ( Constraints_getStatus( &(THIS->constraints),number ) == ST_INACTIVE )
		return THROWERROR( RET_CONSTRAINT_NOT_ACTIVE );

	if ( ( number_idx < 0 ) || ( number_idx >= nAC ) )
		return THROWERROR( RET_CONSTRAINT_NOT_ACTIVE );


	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );

	/* N) PERFORM QPOASES_ZERO CURVATURE TEST. */
	if (ensureNZC == BT_TRUE)
	{
		returnvalue = QProblem_ensureNonzeroCurvature( THIS,BT_FALSE, number, &exchangeHappened,&addBoundNotConstraint,&addIdx,&addStatus);

		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
	}

	/* I) REMOVE <number>th ROW FROM T,
	 *    i.e. shift rows number+1 through nAC  upwards (instead of the actual
	 *    constraint number its corresponding index within matrix A is used). */
	if ( number_idx < nAC-1 )
	{
		for( i=(number_idx+1); i<nAC; ++i )
			for( j=(nAC-i-1); j<nAC; ++j )
				TT(i-1,tcol+j) = TT(i,tcol+j);
		/* gimmick: write zeros into the last row of T */
		for( j=0; j<nAC; ++j )
			TT(nAC-1,tcol+j) = 0.0;


		/* II) RESTORE TRIANGULAR FORM OF T,
		 *     use column-wise Givens rotations to restore reverse triangular form
		 *     of T simultanenous change of Q (i.e. Y). */

		for( j=(nAC-2-number_idx); j>=0; --j )
		{
			QProblemB_computeGivens( TT(nAC-2-j,tcol+1+j),TT(nAC-2-j,tcol+j), &TT(nAC-2-j,tcol+1+j),&TT(nAC-2-j,tcol+j),&c,&s );
			nu = s/(1.0+c);

			for( i=(nAC-j-1); i<(nAC-1); ++i )
				QProblemB_applyGivens( c,s,nu,TT(i,tcol+1+j),TT(i,tcol+j), &TT(i,tcol+1+j),&TT(i,tcol+j) );

			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				QProblemB_applyGivens( c,s,nu,QQ(ii,nZ+1+j),QQ(ii,nZ+j), &QQ(ii,nZ+1+j),&QQ(ii,nZ+j) );
			}
		}
	}
	else
	{
		/* gimmick: write zeros into the last row of T */
		for( j=0; j<nAC; ++j )
			TT(nAC-1,tcol+j) = 0.0;
	}


	if ( ( updateCholesky == BT_TRUE ) &&
		 ( THIS->hessianType != HST_ZERO ) && ( THIS->hessianType != HST_IDENTITY ) )
	{
		/* III) UPDATE CHOLESKY DECOMPOSITION,
		 *      calculate new additional column (i.e. [r sqrt(rho2)]')
		 *      of the Cholesky factor R. */

		/* 1) Calculate Hz = H*z, where z is the new rightmost column of Z
		 *    (i.e. the old leftmost column of Y).  */
		for( j=0; j<nFR; ++j )
			z[j] = QQ(FR_idx[j],nZ);
		DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFree( &(THIS->bounds) ), 1, 1.0, z, nFR, 0.0, Hz, nFR, BT_TRUE);

		if ( nZ > 0 )
		{
			for ( i=0; i<nZ; ++i )
				ZHz[i] = 0.0;
			
			/* 2) Calculate ZHz = Z'*Hz (old Z). */
			for( j=0; j<nFR; ++j )
			{
				jj = FR_idx[j];

				for( i=0; i<nZ; ++i )
					ZHz[i] += QQ(jj,i) * Hz[j];
			}

			/* 3) Calculate r = R^-T * ZHz. */
			if ( QProblem_backsolveR( THIS,ZHz,BT_TRUE,r ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_REMOVECONSTRAINT_FAILED );

			/* 4) Calculate rho2 = rho^2 = z'*Hz - r'*r
			 *    and store r into R. */
			for( i=0; i<nZ; ++i )
			{
				rho2 -= r[i]*r[i];
				RR(i,nZ) = r[i];
			}
		}

		/* 5) Store rho into R. */
		for( j=0; j<nFR; ++j )
			rho2 += QQ(FR_idx[j],nZ) * Hz[j];

		if ( exchangeHappened == BT_FALSE )
		{
			if ( rho2 > QPOASES_ZERO )
				RR(nZ,nZ) = qpOASES_getSqrt( rho2 );
			else
			{
				if ( allowFlipping == BT_FALSE )
				{
					RR(nZ,nZ) = 100.0*QPOASES_EPS;
				}
				else
				{
					THIS->hessianType = HST_SEMIDEF;
					return THROWERROR( RET_HESSIAN_NOT_SPD );
				}
			}
		}
	}


	/* IV) UPDATE INDICES */
	THIS->tabularOutput.idxRemC = number;
	if ( hasFlipped == BT_FALSE )
	{
		if ( Constraints_moveActiveToInactive( &(THIS->constraints),number ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
	}

	if (exchangeHappened == BT_TRUE)
	{
		/* add bound or constraint */

		/* THIS->hessianType = HST_SEMIDEF; */
		RR(nZ,nZ) = 0.0;

		if ( addBoundNotConstraint )
		{
			QProblem_addBound( THIS,addIdx, addStatus, BT_TRUE, BT_FALSE );
			THIS->tabularOutput.excAddB = 1;
		}
		else
		{
			QProblem_addConstraint( THIS,addIdx, addStatus, BT_TRUE, BT_FALSE );
			THIS->tabularOutput.excAddC = 1;
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e B o u n d
 */
returnValue QProblem_removeBound(	QProblem* THIS, int number,
									BooleanType updateCholesky,
									BooleanType allowFlipping,
									BooleanType ensureNZC
									)
{
	int i, j, ii, jj;
	returnValue returnvalue = SUCCESSFUL_RETURN;
	int addIdx;
	BooleanType addBoundNotConstraint;
	SubjectToStatus addStatus;
	BooleanType exchangeHappened = BT_FALSE;

	/* some definitions */
	int nFR = QProblem_getNFR( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nZ  = QProblem_getNZ( THIS );

	int tcol = THIS->sizeT - nAC;
	
	int *FR_idx, *AC_idx;
	
	int nnFRp1;
	
	myStatic real_t tmp[NCMAX];
	real_t c, s, nu;
	real_t z2, rho2;

	myStatic real_t Hz[NVMAX+1];
	myStatic real_t z[NVMAX+1];

	myStatic real_t r[NVMAX];
	myStatic real_t rhs[NVMAX];


	/* consistency checks */
	if ( Bounds_getStatus( &(THIS->bounds),number ) == ST_INACTIVE )
		return THROWERROR( RET_BOUND_NOT_ACTIVE );

	if ( ( QProblem_getStatus( THIS ) == QPS_NOTINITIALISED )    ||
		 ( QProblem_getStatus( THIS ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( QProblem_getStatus( THIS ) == QPS_HOMOTOPYQPSOLVED )  ||
 		 ( QProblem_getStatus( THIS ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* 0) PERFORM ZERO CURVATURE TEST. */
	if (ensureNZC == BT_TRUE)
	{
		returnvalue = QProblem_ensureNonzeroCurvature( THIS,BT_TRUE, number, &exchangeHappened,&addBoundNotConstraint,&addIdx,&addStatus);

		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
	}

	/* I) UPDATE INDICES */
	THIS->tabularOutput.idxRemB = number;
	if ( Bounds_moveFixedToFree( &(THIS->bounds),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_REMOVEBOUND_FAILED );

	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );

	/* I) APPEND <nFR+1>th UNITY VECTOR TO Q. */
	nnFRp1 = FR_idx[nFR];
	for( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		QQ(ii,nFR) = 0.0;
		QQ(nnFRp1,i) = 0.0;
	}
	QQ(nnFRp1,nFR) = 1.0;

	if ( nAC > 0 )
	{
		/* store new column a in a temporary vector instead of shifting T one column to the left and appending a */
		Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );

		DenseMatrix_getCol(THIS->A,number, Constraints_getActive( &(THIS->constraints)), 1.0, tmp);


		/* II) RESTORE TRIANGULAR FORM OF T,
		 *     use column-wise Givens rotations to restore reverse triangular form
		 *     of T = [T A(:,number)], simultanenous change of Q (i.e. Y and Z). */
		for( j=(nAC-1); j>=0; --j )
		{
			QProblemB_computeGivens( tmp[nAC-1-j],TT(nAC-1-j,tcol+j), &TT(nAC-1-j,tcol+j),&tmp[nAC-1-j],&c,&s );
			nu = s/(1.0+c);

			for( i=(nAC-j); i<nAC; ++i )
				QProblemB_applyGivens( c,s,nu,tmp[i],TT(i,tcol+j), &TT(i,tcol+j),&tmp[i] );

			for( i=0; i<=nFR; ++i )
			{
				ii = FR_idx[i];
				/* nZ+1+nAC = nFR+1  /  nZ+(1) = nZ+1 */
				QProblemB_applyGivens( c,s,nu,QQ(ii,nZ+1+j),QQ(ii,nZ+j), &QQ(ii,nZ+1+j),&QQ(ii,nZ+j) );
			}
		}
	}


	if ( ( updateCholesky == BT_TRUE ) &&
		 ( THIS->hessianType != HST_ZERO ) && ( THIS->hessianType != HST_IDENTITY ) )
	{
		/* III) UPDATE CHOLESKY DECOMPOSITION,
		 *      calculate new additional column (i.e. [r sqrt(rho2)]')
		 *      of the Cholesky factor R: */
		z2 = QQ(nnFRp1,nZ);
		rho2 = DenseMatrix_diag(THIS->H,nnFRp1)*z2*z2; /* rho2 = h2*z2*z2 */

		if ( nFR > 0 )
		{
			/* Attention: Index list of free variables has already grown by one! */
			/* 1) Calculate R'*r = Zfr'*Hfr*z1 + z2*Zfr'*h1 =: Zfr'*Hz + z2*Zfr'*h1 =: rhs and
			 *    rho2 = z1'*Hfr*z1 + 2*z2*h1'*z1 + h2*z2^2 - r'*r =: z1'*Hz + 2*z2*h1'*z1 + h2*z2^2 - r'r */
			for( j=0; j<nFR; ++j )
				z[j] = QQ(FR_idx[j],nZ);
			z[nFR] = 0.0;
			DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFree( &(THIS->bounds) ), 1, 1.0, z, nFR+1, 0.0, Hz, nFR+1, BT_TRUE);

			DenseMatrix_getCol(THIS->H,nnFRp1, Bounds_getFree( &(THIS->bounds) ), 1.0, z);

			if ( nZ > 0 )
			{
				for( i=0; i<nZ; ++i )
					rhs[i] = 0.0;

				/* 2) Calculate rhs. */
				for( j=0; j<nFR; ++j )
				{
					jj = FR_idx[j];
					for( i=0; i<nZ; ++i )
										/* Zfr' * ( Hz + z2*h1 ) */
						rhs[i] += QQ(jj,i) * ( Hz[j] + z2 * z[j] );
				}

				/* 3) Calculate r = R^-T * rhs. */
				if ( QProblem_backsolveRrem( THIS,rhs,BT_TRUE,BT_TRUE,r ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_REMOVEBOUND_FAILED );


				/* 4) Calculate rho2 = rho^2 = z'*Hz - r'*r
				 *    and store r into R. */
				for( i=0; i<nZ; ++i )
				{
					rho2 -= r[i]*r[i];
					RR(i,nZ) = r[i];
				}
			}

			for( j=0; j<nFR; ++j )
			{
				jj = FR_idx[j];
							/* z1' * ( Hz + 2*z2*h1 ) */
				rho2 += QQ(jj,nZ) * ( Hz[j] + 2.0*z2*z[j] );
			}
		}

		/* 5) Store rho into R. */
		if ( exchangeHappened == BT_FALSE )
		{
			if ( rho2 > QPOASES_ZERO )
				RR(nZ,nZ) = qpOASES_getSqrt( rho2 );
			else
			{
				if ( allowFlipping == BT_FALSE )
					RR(nZ,nZ) = 100.0*QPOASES_EPS;
				else
				{
					THIS->hessianType = HST_SEMIDEF;
					return THROWERROR( RET_HESSIAN_NOT_SPD );
				}
			}
		}
		else
		{
			/* add bound or constraint */

			/* THIS->hessianType = HST_SEMIDEF; */
			RR(nZ,nZ) = 0.0;

			if ( addBoundNotConstraint )
			{
				QProblem_addBound( THIS,addIdx, addStatus, BT_TRUE, BT_FALSE);
				THIS->tabularOutput.excAddB = 1;
			}
			else
			{
				QProblem_addConstraint(THIS,addIdx, addStatus, BT_TRUE, BT_FALSE);
				THIS->tabularOutput.excAddC = 1;
			}
		}
	}

	return SUCCESSFUL_RETURN;
}



returnValue QProblem_performPlainRatioTest(	QProblem* THIS,
											int nIdx,
											const int* const idxList,
											const real_t* const num,
											const real_t* const den,
											real_t epsNum,
											real_t epsDen,
											real_t* t,
											int* BC_idx
											)
{
	int i;
	for (i = 0; i < nIdx; i++)
		if ( (num[i] > epsNum) && (den[i] > epsDen) && ((*t) * den[i] > num[i]) )
		{
			*t = num[i] / den[i];
			*BC_idx = idxList[i];
		}
		
	return SUCCESSFUL_RETURN;
}


returnValue QProblem_ensureNonzeroCurvature(	QProblem* THIS,
												BooleanType removeBoundNotConstraint,
												int remIdx,
												BooleanType* exchangeHappened,
												BooleanType* addBoundNotConstraint,
												int* addIdx,
												SubjectToStatus* addStatus
												)
{
	int i, ii;
	int addLBndIdx = -1, addLCnstrIdx = -1, addUBndIdx = -1, addUCnstrIdx = -1; /* exchange indices */
	int *FX_idx, *AC_idx, *IAC_idx, *FR_idx;
	returnValue returnvalue = SUCCESSFUL_RETURN;
 
	int nV  = QProblem_getNV( THIS );
	int nFR = QProblem_getNFR( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nC  = QProblem_getNC( THIS );
	int nFX = QProblem_getNFX( THIS );
	int nIAC = QProblem_getNIAC( THIS );

	myStatic real_t delta_xFX[NVMAX];
	myStatic real_t delta_xFR[NVMAX];
	myStatic real_t delta_yAC[NCMAX];
	myStatic real_t delta_yFX[NVMAX];

	int dim;
	myStatic real_t nul[NVCMAX];
	myStatic real_t ek[NVMAX]; /* minus e_k (bound k is removed) */

	real_t one = 1.0;

	real_t a;
	real_t normXi = 0.0;
	real_t normS = 0.0;
	real_t sigmaLBnd, sigmaLCnstr, sigmaUBnd, sigmaUCnstr, sigma;
	myStatic real_t x_W[NVMAX];
	
	myStatic real_t As[NCMAX];
	myStatic real_t Ax_W[NCMAX];

	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );
	Indexlist_getNumberArray( Constraints_getInactive( &(THIS->constraints) ),&IAC_idx );

	*addBoundNotConstraint = BT_TRUE;
	*addStatus = ST_INACTIVE;
	*exchangeHappened = BT_FALSE;

	if (removeBoundNotConstraint)
	{
		dim = nV < nC ? nC : nV;
		for (ii = 0; ii < dim; ++ii)
			nul[ii]=0.0;
		for (ii = 0; ii < nV; ++ii)
			ek[ii]=0.0;
		ek[remIdx] = Bounds_getStatus( &(THIS->bounds),remIdx) == ST_LOWER ? 1.0 : -1.0;

		returnvalue = QProblem_determineStepDirection(	THIS,nul, nul, nul, ek, ek,
														BT_FALSE, BT_FALSE,
														delta_xFX, delta_xFR, delta_yAC, delta_yFX
														);
	}
	else
	{
		for (ii = 0; ii < nV; ++ii)
			nul[ii]=0.0;
		for (ii = 0; ii < nC; ++ii)
			ek[ii]=0.0;
		ek[remIdx] = Constraints_getStatus( &(THIS->constraints),remIdx) == ST_LOWER ? 1.0 : -1.0;

		returnvalue = QProblem_determineStepDirection(	THIS,nul,
														ek, ek, nul, nul,
														BT_FALSE, BT_TRUE,
														delta_xFX, delta_xFR, delta_yAC, delta_yFX
														);
	}

	/* compute the weight in inf-norm */
	for (ii = 0; ii < nAC; ++ii)
	{
		a = qpOASES_getAbs (delta_yAC[ii]);
		if (normXi < a) normXi = a;
	}
	for (ii = 0; ii < nFX; ++ii)
	{
		a = qpOASES_getAbs (delta_yFX[ii]);
		if (normXi < a) normXi = a;
	}

	/* look at the "zero" in a relative inf-norm */
	for (ii = 0; ii < nFX; ++ii)
	{
		a = qpOASES_getAbs (delta_xFX[ii]);
		if (normS < a) normS = a;
	}
	for (ii = 0; ii < nFR; ++ii)
	{
		a = qpOASES_getAbs (delta_xFR[ii]);
		if (normS < a) normS = a;
	}

	/* relative test against zero in inf-norm */
	if (normXi < THIS->options.epsNZCTests * normS)
	{
		/* determine jump in x via ratio tests */

		/* bounds */

		/* compress x-u */
		for (i = 0; i < nFR; i++)
		{
			ii = FR_idx[i];
			x_W[i] = THIS->ub[ii] - THIS->x[ii];
		}
		/* performRatioTest( nFR,FR_idx,&(THIS->bounds),x_W,delta_xFR, THIS->options.epsNum,THIS->options.epsDen, sigmaUBnd,addUBndIdx ); */
		sigmaUBnd = THIS->options.maxPrimalJump;
		addUBndIdx = -1;
		QProblem_performPlainRatioTest(THIS,nFR, FR_idx, x_W, delta_xFR, THIS->options.epsNum, THIS->options.epsDen, &sigmaUBnd,&addUBndIdx);
		if (removeBoundNotConstraint == BT_TRUE && Bounds_getStatus( &(THIS->bounds),remIdx) == ST_LOWER)
		{
			/* also consider bound which is to be removed */
			one = 1.0;
			x_W[0] = THIS->ub[remIdx] - THIS->x[remIdx];
			QProblem_performPlainRatioTest(THIS,1, &remIdx, x_W, &one, THIS->options.epsNum, THIS->options.epsDen, &sigmaUBnd,&addUBndIdx);
		}

		/* compress x-l */
		for (i = 0; i < nFR; i++)
		{
			ii = FR_idx[i];
			x_W[i] = THIS->x[ii] - THIS->lb[ii];
		}
		for (i = 0; i < nFR; i++)
			delta_xFR[i] = -delta_xFR[i];
		/* performRatioTest( nFR,FR_idx,&(THIS->bounds),x_W,delta_xFR, THIS->options.epsNum,THIS->options.epsDen, sigmaLBnd,addLBndIdx ); */
		sigmaLBnd = THIS->options.maxPrimalJump;
		addLBndIdx = -1;
		QProblem_performPlainRatioTest(THIS,nFR, FR_idx, x_W, delta_xFR, THIS->options.epsNum, THIS->options.epsDen, &sigmaLBnd,&addLBndIdx);
		if (removeBoundNotConstraint == BT_TRUE && Bounds_getStatus( &(THIS->bounds),remIdx) == ST_UPPER)
		{
			/* also consider bound which is to be removed */
			one = 1.0;
			x_W[0] = THIS->x[remIdx] - THIS->lb[remIdx];
			QProblem_performPlainRatioTest(THIS,1, &remIdx, x_W, &one, THIS->options.epsNum, THIS->options.epsDen, &sigmaLBnd,&addLBndIdx);
		}
		for (i = 0; i < nFR; i++)
			delta_xFR[i] = -delta_xFR[i];

		/* constraints */

		/* compute As (compressed to inactive constraints) */
		DenseMatrix_subTimes(THIS->A,Constraints_getInactive(&(THIS->constraints)), Bounds_getFixed( &(THIS->bounds) ), 1, 1.0, delta_xFX, nFX, 0.0, As, nIAC, BT_TRUE);
		DenseMatrix_subTimes(THIS->A,Constraints_getInactive(&(THIS->constraints)), Bounds_getFree( &(THIS->bounds) ), 1, 1.0, delta_xFR, nFR, 1.0, As, nIAC, BT_TRUE);

		/* compress Ax_u */
		for (i = 0; i < nIAC; i++)
		{
			ii = IAC_idx[i];
			Ax_W[i] = THIS->Ax_u[ii];
		}
		/* performRatioTest( nIAC,IAC_idx,&(THIS->constraints), Ax_W,As, THIS->options.epsNum,THIS->options.epsDen, sigmaUCnstr,addUCnstrIdx ); */
		sigmaUCnstr = THIS->options.maxPrimalJump;
		addUCnstrIdx = -1;
		QProblem_performPlainRatioTest(THIS,nIAC, IAC_idx, Ax_W, As, THIS->options.epsNum, THIS->options.epsDen, &sigmaUCnstr,&addUCnstrIdx);
		if (removeBoundNotConstraint == BT_FALSE && Constraints_getStatus( &(THIS->constraints),remIdx) == ST_LOWER)
		{
			/* also consider constraint which is to be removed */
			one = 1.0;
			QProblem_performPlainRatioTest(THIS,1, &remIdx, &THIS->Ax_u[remIdx], &one, THIS->options.epsNum, THIS->options.epsDen, &sigmaUCnstr,&addUCnstrIdx);
		}

		/* compress Ax_l */
		for (i = 0; i < nIAC; i++)
		{
			ii = IAC_idx[i];
			Ax_W[i] = THIS->Ax_l[ii];
		}
		for (i = 0; i < nIAC; i++)
			As[i] = -As[i];
		/* performRatioTest( nIAC,IAC_idx,&(THIS->constraints), Ax_W,As, THIS->options.epsNum,THIS->options.epsDen, sigmaLCnstr,addLCnstrIdx ); */
		sigmaLCnstr = THIS->options.maxPrimalJump;
		addLCnstrIdx = -1;
		QProblem_performPlainRatioTest(THIS,nIAC, IAC_idx, Ax_W, As, THIS->options.epsNum, THIS->options.epsDen, &sigmaLCnstr,&addLCnstrIdx);
		if (removeBoundNotConstraint == BT_FALSE && Constraints_getStatus( &(THIS->constraints),remIdx) == ST_UPPER)
		{
			/* also consider constraint which is to be removed */
			one = 1.0;
			QProblem_performPlainRatioTest(THIS,1, &remIdx, &THIS->Ax_l[remIdx], &one, THIS->options.epsNum, THIS->options.epsDen, &sigmaLCnstr,&addLCnstrIdx);
		}

		/* perform primal jump */
		sigma = THIS->options.maxPrimalJump;
		if (sigmaUCnstr < sigma) { sigma = sigmaUCnstr; *addStatus = ST_UPPER; *addBoundNotConstraint = BT_FALSE; *addIdx = addUCnstrIdx; }
		if (sigmaLCnstr < sigma) { sigma = sigmaLCnstr; *addStatus = ST_LOWER; *addBoundNotConstraint = BT_FALSE; *addIdx = addLCnstrIdx; }
		if (sigmaUBnd < sigma) { sigma = sigmaUBnd; *addStatus = ST_UPPER; *addBoundNotConstraint = BT_TRUE; *addIdx = addUBndIdx; }
		if (sigmaLBnd < sigma) { sigma = sigmaLBnd; *addStatus = ST_LOWER; *addBoundNotConstraint = BT_TRUE; *addIdx = addLBndIdx; }

		if (sigma >= THIS->options.maxPrimalJump)
		{
			THIS->unbounded = BT_TRUE;
			returnvalue = RET_HOTSTART_STOPPED_UNBOUNDEDNESS;
		}
		else
		{
			for (i = 0; i < nFR; i++)
				THIS->x[FR_idx[i]] += sigma * delta_xFR[i];

			for (i = 0; i < nFX; i++)
				THIS->x[FX_idx[i]] += sigma * delta_xFX[i];

			/* update Ax, Ax_u, and Ax_l */
			DenseMatrix_times(THIS->A,1, 1.0, THIS->x, nV, 0.0, THIS->Ax, nC);
			for (i = 0; i < nC; i++) THIS->Ax_u[i] = THIS->ubA[i] - THIS->Ax[i];
			for (i = 0; i < nC; i++) THIS->Ax_l[i] = THIS->Ax[i] - THIS->lbA[i];

			/* change working set later */
			*exchangeHappened = BT_TRUE;
		}
	}

	return returnvalue;
}



/*
 *	b a c k s o l v e T
 */
returnValue QProblem_backsolveT( QProblem* THIS, const real_t* const b, BooleanType transposed, real_t* const a )
{
	int i, j;
	int nT = QProblem_getNAC( THIS );
	int tcol = THIS->sizeT - nT;

	real_t sum;

	/* nothing to do */
	if ( nT <= 0 )
		return SUCCESSFUL_RETURN;


	/* Solve Ta = b, where T might be transposed. */
	if ( transposed == BT_FALSE )
	{
		/* solve Ta = b */
		for( i=0; i<nT; ++i )
		{
			sum = b[i];
			for( j=0; j<i; ++j )
				sum -= TT(i,THIS->sizeT-1-j) * a[nT-1-j];

			if ( qpOASES_getAbs( TT(i,THIS->sizeT-1-i) ) > QPOASES_EPS )
				a[nT-1-i] = sum / TT(i,THIS->sizeT-1-i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}
	else
	{
		/* solve T^T*a = b */
		for( i=0; i<nT; ++i )
		{
			sum = b[i];
			for( j=0; j<i; ++j )
				sum -= TT(nT-1-j,tcol+i) * a[nT-1-j];

			if ( qpOASES_getAbs( TT(nT-1-i,tcol+i) ) > QPOASES_EPS )
				a[nT-1-i] = sum / TT(nT-1-i,tcol+i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}


	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e D a t a S h i f t
 */
returnValue QProblem_determineDataShift(	QProblem* THIS, const real_t* const g_new, const real_t* const lbA_new, const real_t* const ubA_new,
											const real_t* const lb_new, const real_t* const ub_new,
											real_t* const delta_g, real_t* const delta_lbA, real_t* const delta_ubA,
											real_t* const delta_lb, real_t* const delta_ub,
											BooleanType* Delta_bC_isZero, BooleanType* Delta_bB_isZero
											)
{
	int i, ii;
	int nC  = QProblem_getNC( THIS );
	int nAC = QProblem_getNAC( THIS );
	
	int* FX_idx;
	int* AC_idx;

	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );



	/* I) DETERMINE DATA SHIFT FOR BOUNDS */
	QProblemBCPY_determineDataShift(	THIS,g_new,lb_new,ub_new,
										delta_g,delta_lb,delta_ub,
										Delta_bB_isZero );


	/* II) DETERMINE DATA SHIFT FOR CONSTRAINTS */
	/* 1) Calculate shift directions. */
	for( i=0; i<nC; ++i )
	{
		/* if lower constraints' bounds are to be disabled or do not exist, shift them to -infinity */
		if ( lbA_new != 0 )
			delta_lbA[i] = lbA_new[i] - THIS->lbA[i];
		else
			delta_lbA[i] = -QPOASES_INFTY - THIS->lbA[i];
	}

	for( i=0; i<nC; ++i )
	{
		/* if upper constraints' bounds are to be disabled or do not exist, shift them to infinity */
		if ( ubA_new != 0 )
			delta_ubA[i] = ubA_new[i] - THIS->ubA[i];
		else
			delta_ubA[i] = QPOASES_INFTY - THIS->ubA[i];
	}

	/* 2) Determine if active constraints' bounds are to be shifted. */
	*Delta_bC_isZero = BT_TRUE;

	for ( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];

		if ( ( qpOASES_getAbs( delta_lbA[ii] ) > QPOASES_EPS ) || ( qpOASES_getAbs( delta_ubA[ii] ) > QPOASES_EPS ) )
		{
			*Delta_bC_isZero = BT_FALSE;
			break;
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e S t e p D i r e c t i o n
 */
returnValue QProblem_determineStepDirection(	QProblem* THIS, const real_t* const delta_g, const real_t* const delta_lbA, const real_t* const delta_ubA,
												const real_t* const delta_lb, const real_t* const delta_ub,
												BooleanType Delta_bC_isZero, BooleanType Delta_bB_isZero,
												real_t* const delta_xFX, real_t* const delta_xFR,
												real_t* const delta_yAC, real_t* const delta_yFX
												)
{
	int i, j, ii, jj, r;
	int nFR = QProblem_getNFR( THIS );
	int nFX = QProblem_getNFX( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nZ  = QProblem_getNZ( THIS );
	
	int* FR_idx;
	int* FX_idx;
	int* AC_idx;
	
	real_t rnrm = 0.0;

	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );


	/* I) DETERMINE delta_xFX (THIS is exact, does not need refinement) */
	if ( Delta_bB_isZero == BT_FALSE )
	{
		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];

			if ( Bounds_getStatus( &(THIS->bounds),ii ) == ST_LOWER )
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


	/* tempA and tempB hold the residuals in gFR and bA (= lbA or ubA)
	 * delta_xFR, delta_yAC hold the steps that get refined */
	for ( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		THIS->tempA[i] = delta_g[ii];
		delta_xFR[i] = 0.0;
	}
	for ( i=0; i<nAC; ++i )
		delta_yAC[i] = 0.0;
	if ( Delta_bC_isZero == BT_FALSE )
	{
		for ( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			if ( Constraints_getStatus( &(THIS->constraints),ii ) == ST_LOWER )
				THIS->tempB[i] = delta_lbA[ii];
			else
				THIS->tempB[i] = delta_ubA[ii];
		}
	}
	else
	{
		for ( i=0; i<nAC; ++i )
			THIS->tempB[i] = 0.0;
	}

	/* Iterative refinement loop for delta_xFRz, delta_xFRy, delta_yAC */
	for ( r=0; r<=THIS->options.numRefinementSteps; ++r )
	{
		/* II) DETERMINE delta_xFR */
		if ( nFR > 0 )
		{
			for( i=0; i<nFR; ++i )
				THIS->delta_xFR_TMP[i] = 0.0;

			/* 1) Determine delta_xFRy. */
			if ( nAC > 0 )
			{
				if ( ( Delta_bC_isZero == BT_TRUE ) && ( Delta_bB_isZero == BT_TRUE ) )
				{
					for( i=0; i<nAC; ++i )
						THIS->delta_xFRy[i] = 0.0;
				}
				else
				{
					/* compute bA - A * delta_xFX. tempB already holds bA->
					 * in refinements r>=1, delta_xFX is exactly zero */
					if ( ( Delta_bB_isZero == BT_FALSE ) && ( r == 0 ) )
						DenseMatrix_subTimes(THIS->A,Constraints_getActive( &(THIS->constraints)), Bounds_getFixed( &(THIS->bounds) ), 1, -1.0, delta_xFX, nFX, 1.0, THIS->tempB, nAC, BT_TRUE);

					if ( QProblem_backsolveT( THIS,THIS->tempB, BT_FALSE, THIS->delta_xFRy ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_STEPDIRECTION_FAILED_TQ );

					for( i=0; i<nFR; ++i )
					{
						ii = FR_idx[i];
						for( j=0; j<nAC; ++j )
							THIS->delta_xFR_TMP[i] += QQ(ii,nZ+j) * THIS->delta_xFRy[j];
					}
				}
			}


			/* 2) Determine delta_xFRz. */
			for( i=0; i<nZ; ++i )
				THIS->delta_xFRz[i] = 0.0;

			if ( ( THIS->hessianType == HST_ZERO ) || ( THIS->hessianType == HST_IDENTITY ) )
			{
				/* compute Z*delta_gFR [/eps] (delta_gFR is stored in tempA) */
				for( j=0; j<nFR; ++j )
				{
					jj = FR_idx[j];
					for( i=0; i<nZ; ++i )
						THIS->delta_xFRz[i] -= QQ(jj,i) * THIS->tempA[j];
				}

				if ( THIS->hessianType == HST_ZERO )
				{
					if ( QProblem_usingRegularisation( THIS ) == BT_TRUE )
					{
						for( i=0; i<nZ; ++i )
							THIS->delta_xFRz[i] /= THIS->regVal;
					}
					else
					{
						/* TODO: What if no regularisation?? */
					}
				}
			}
			else
			{
				/* compute HMX*delta_xFX. DESTROY delta_gFR that was in tempA */
				if ( ( Delta_bB_isZero == BT_FALSE ) && ( r == 0 ) )
					DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFixed( &(THIS->bounds) ), 1, 1.0, delta_xFX, nFX, 1.0, THIS->tempA, nFR, BT_TRUE);

				/* compute HFR*delta_xFRy */
				if ( ( nAC > 0 ) && ( ( Delta_bC_isZero == BT_FALSE ) || ( Delta_bB_isZero == BT_FALSE ) ) )
					DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFree( &(THIS->bounds) ), 1, 1.0, THIS->delta_xFR_TMP, nFR, 1.0, THIS->tempA, nFR, BT_TRUE);

				/* compute ZFR_delta_xFRz = (Z'*HFR*Z) \ Z * (HFR*delta_xFR + HMX*delta_xFX + delta_gFR) */
				if ( nZ > 0 )
				{
					for( j=0; j<nFR; ++j )
					{
						jj = FR_idx[j];
						for( i=0; i<nZ; ++i )
							THIS->delta_xFRz[i] -= QQ(jj,i) * THIS->tempA[j];
					}

					if ( QProblem_backsolveR( THIS,THIS->delta_xFRz,BT_TRUE,THIS->delta_xFRz ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );

					if ( QProblem_backsolveR( THIS,THIS->delta_xFRz,BT_FALSE,THIS->delta_xFRz ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );
				}
			}

			/* compute Z * ZFR_delta_xFRz */
			if ( nZ > 0 )
			{
				for( i=0; i<nFR; ++i )
				{
					THIS->ZFR_delta_xFRz[i] = 0.0;

					ii = FR_idx[i];
					for( j=0; j<nZ; ++j )
						THIS->ZFR_delta_xFRz[i] += QQ(ii,j) * THIS->delta_xFRz[j];

					THIS->delta_xFR_TMP[i] += THIS->ZFR_delta_xFRz[i];
				}
			}
		}

		/* III) DETERMINE delta_yAC */
		if ( nAC > 0 ) /* => ( nFR = nZ + nAC > 0 ) */
		{
			if ( ( THIS->hessianType == HST_ZERO ) || ( THIS->hessianType == HST_IDENTITY ) )
			{
				/* if zero:     delta_yAC = (T')^-1 * ( Yfr*delta_gFR + eps*delta_xFRy ),
				 * if identity: delta_yAC = (T')^-1 * ( Yfr*delta_gFR +     delta_xFRy )
				 *
				 * DESTROY residual_bA that was stored in tempB
				 * If we come here, residual_gFR in tempA is STILL VALID
				 */
				if ( THIS->hessianType == HST_IDENTITY )
				{
					for( j=0; j<nAC; ++j )
						THIS->tempB[j] = THIS->delta_xFRy[j];
				}
				else /* hessianType == HST_ZERO */
				{
					if ( QProblem_usingRegularisation( THIS ) == BT_TRUE )
					{
						for( j=0; j<nAC; ++j )
							THIS->tempB[j] = THIS->regVal * THIS->delta_xFRy[j];
					}
					else
					{
						for( j=0; j<nAC; ++j )
							THIS->tempB[j] = 0.0;
					}
				}

				for( j=0; j<nAC; ++j )
				{
					for( i=0; i<nFR; ++i )
					{
						ii = FR_idx[i];
						THIS->tempB[j] += QQ(ii,nZ+j) * THIS->tempA[i];
					}
				}
			}
			else
			{
				/* Compute HFR * delta_xFR + HMX*delta_xFX
				 * Here, tempA holds (HFR*delta_xFRy + HMX*delta_xFX) */
				if ( nZ > 0 )
					DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFree( &(THIS->bounds) ), 1, 1.0, THIS->ZFR_delta_xFRz, nFR, 1.0, THIS->tempA, nFR, BT_TRUE);

				for( i=0; i<nAC; ++i)
				{
					THIS->tempB[i] = 0.0;
					for( j=0; j<nFR; ++j )
					{
						jj = FR_idx[j];
						THIS->tempB[i] += QQ(jj,nZ+i) * THIS->tempA[j];
					}
				}
			}

			if ( QProblem_backsolveT( THIS,THIS->tempB,BT_TRUE,THIS->delta_yAC_TMP ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_STEPDIRECTION_FAILED_TQ );
		}

		/* refine the solution found so far */
		for ( i=0; i<nFR; ++i )
			delta_xFR[i] += THIS->delta_xFR_TMP[i];
		for ( i=0; i<nAC; ++i )
			delta_yAC[i] += THIS->delta_yAC_TMP[i];

		if ( THIS->options.numRefinementSteps > 0 )
		{
			/* compute residuals in tempA and tempB, and max-norm */
			for ( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				THIS->tempA[i] = delta_g[ii];
			}

			switch ( THIS->hessianType )
			{
				case HST_ZERO:
					if ( QProblem_usingRegularisation( THIS ) == BT_TRUE )
						for ( i=0; i<nFR; ++i )
							THIS->tempA[i] += THIS->regVal*delta_xFR[i];
					break;

				case HST_IDENTITY:
					for ( i=0; i<nFR; ++i )
						THIS->tempA[i] += delta_xFR[i];
					break;

				default:
					DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFree( &(THIS->bounds) ),  1, 1.0, delta_xFR, nFR, 1.0, THIS->tempA, nFR, BT_TRUE);
					DenseMatrix_subTimes(THIS->H,Bounds_getFree( &(THIS->bounds) ), Bounds_getFixed( &(THIS->bounds) ), 1, 1.0, delta_xFX, nFX, 1.0, THIS->tempA, nFR, BT_TRUE);
					break;
			}

			DenseMatrix_subTransTimes(THIS->A,Constraints_getActive( &(THIS->constraints)), Bounds_getFree( &(THIS->bounds) ), 1, -1.0, delta_yAC, nAC, 1.0, THIS->tempA, nFR);
			rnrm = 0.0;
			for ( i=0; i<nFR; ++i )
				if (rnrm < qpOASES_getAbs (THIS->tempA[i]))
					rnrm = qpOASES_getAbs (THIS->tempA[i]);

			if (!Delta_bC_isZero)
			{
				for ( i=0; i<nAC; ++i )
				{
					ii = AC_idx[i];
					if ( Constraints_getStatus( &(THIS->constraints),ii ) == ST_LOWER )
						THIS->tempB[i] = delta_lbA[ii];
					else
						THIS->tempB[i] = delta_ubA[ii];
				}
			}
			else
			{
				for ( i=0; i<nAC; ++i )
					THIS->tempB[i] = 0.0;
			}
			DenseMatrix_subTimes(THIS->A,Constraints_getActive( &(THIS->constraints)), Bounds_getFree( &(THIS->bounds) ), 1, -1.0, delta_xFR, nFR, 1.0, THIS->tempB, nAC, BT_TRUE);
			DenseMatrix_subTimes(THIS->A,Constraints_getActive( &(THIS->constraints)), Bounds_getFixed( &(THIS->bounds) ), 1, -1.0, delta_xFX, nFX, 1.0, THIS->tempB, nAC, BT_TRUE);
			for ( i=0; i<nAC; ++i )
				if (rnrm < qpOASES_getAbs (THIS->tempB[i]))
					rnrm = qpOASES_getAbs (THIS->tempB[i]);

			/* early termination of residual norm small enough */
			if ( rnrm < THIS->options.epsIterRef )
				break;
		}
	} /* end of refinement loop for delta_xFRz, delta_xFRy, delta_yAC */


	/* IV) DETERMINE delta_yFX */
	if ( nFX > 0 )
	{
		for( i=0; i<nFX; ++i )
			delta_yFX[i] = delta_g[FX_idx[i]];

		DenseMatrix_subTransTimes(THIS->A,Constraints_getActive( &(THIS->constraints)), Bounds_getFixed( &(THIS->bounds) ), 1, -1.0, delta_yAC, nAC, 1.0, delta_yFX, nFX);

		if ( THIS->hessianType == HST_ZERO )
		{
			if ( QProblem_usingRegularisation( THIS ) == BT_TRUE )
				for( i=0; i<nFX; ++i )
					delta_yFX[i] += THIS->regVal*delta_xFX[i];
		}
		else if ( THIS->hessianType == HST_IDENTITY )
		{
			for( i=0; i<nFX; ++i )
				delta_yFX[i] += delta_xFX[i];
		}
		else
		{
			DenseMatrix_subTimes(THIS->H,Bounds_getFixed( &(THIS->bounds) ), Bounds_getFree( &(THIS->bounds) ), 1, 1.0, delta_xFR, nFR, 1.0, delta_yFX, nFX, BT_TRUE);
			DenseMatrix_subTimes(THIS->H,Bounds_getFixed( &(THIS->bounds) ), Bounds_getFixed( &(THIS->bounds) ), 1, 1.0, delta_xFX, nFX, 1.0, delta_yFX, nFX, BT_TRUE);
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p e r f o r m S t e p
 */
returnValue QProblem_performStep(	QProblem* THIS, const real_t* const delta_g,
									const real_t* const delta_lbA, const real_t* const delta_ubA,
									const real_t* const delta_lb, const real_t* const delta_ub,
									const real_t* const delta_xFX, const real_t* const delta_xFR,
									const real_t* const delta_yAC, const real_t* const delta_yFX,
									int* BC_idx, SubjectToStatus* BC_status, BooleanType* BC_isBound
									)
{
	int i, j, ii, jj;
	#ifndef __XPCTARGET__
	myStatic char messageString[QPOASES_MAX_STRING_LENGTH];
	#endif

	int nV  = QProblem_getNV( THIS );
	int nC  = QProblem_getNC( THIS );
	int nFR = QProblem_getNFR( THIS );
	int nFX = QProblem_getNFX( THIS );
	int nAC = QProblem_getNAC( THIS );
	int nIAC = QProblem_getNIAC( THIS );
	
	int* FR_idx;
	int* FX_idx;
	int* AC_idx;
	int* IAC_idx;

	int BC_idx_tmp = -1;

	myStatic real_t num[NVCMAX];
	myStatic real_t den[NVCMAX];

	myStatic real_t delta_Ax_l[NCMAX];
	myStatic real_t delta_Ax_u[NCMAX];
	myStatic real_t delta_Ax[NCMAX];

	myStatic real_t delta_x[NVMAX];

	/* initialise maximum steplength array */
	THIS->tau = 1.0;
	*BC_idx = -1;
	*BC_status = ST_UNDEFINED;
	*BC_isBound = BT_FALSE;

	
	Indexlist_getNumberArray( Bounds_getFree( &(THIS->bounds) ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );
	Indexlist_getNumberArray( Constraints_getInactive( &(THIS->constraints) ),&IAC_idx );
	
	for( j=0; j<nFR; ++j )
	{
		jj = FR_idx[j];
		delta_x[jj] = delta_xFR[j];
	}
	for( j=0; j<nFX; ++j )
	{
		jj = FX_idx[j];
		delta_x[jj] = delta_xFX[j];
	}


	/* I) DETERMINE MAXIMUM DUAL STEPLENGTH: */
	/* 1) Ensure that active dual constraints' bounds remain valid
	 *    (ignoring inequality constraints).  */
	for( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];

		num[i] = THIS->y[nV+ii];
		den[i] = -delta_yAC[i];
	}

	QProblem_performRatioTestC( THIS,nAC,AC_idx,&(THIS->constraints), num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

	if ( BC_idx_tmp >= 0 )
	{
		*BC_idx = BC_idx_tmp;
		*BC_status = ST_INACTIVE;
		*BC_isBound = BT_FALSE;
	}


	/* 2) Ensure that active dual bounds remain valid
	 *    (ignoring implicitly fixed variables). */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = THIS->y[ii];
		den[i] = -delta_yFX[i];
	}

	QProblem_performRatioTestB( THIS,nFX,FX_idx,&(THIS->bounds),num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

	if ( BC_idx_tmp >= 0 )
	{
		*BC_idx = BC_idx_tmp;
		*BC_status = ST_INACTIVE;
		*BC_isBound = BT_TRUE;
	}


 	/* II) DETERMINE MAXIMUM PRIMAL STEPLENGTH */
	/* 1) Ensure that inactive constraints' bounds remain valid
	 *    (ignoring unbounded constraints). */
	/* calculate product A*x */
	if ( THIS->constraintProduct == 0 )
	{
		DenseMatrix_subTimes(THIS->A,Constraints_getInactive(&(THIS->constraints)), 0, 1, 1.0, delta_x, nV, 0.0, delta_Ax, nC, BT_FALSE);
	}
/*	else
	{
		for( i=0; i<nIAC; ++i )
		{
			ii = IAC_idx[i];

			if ( Constraints_getType( &(THIS->constraints),ii ) != ST_UNBOUNDED )
			{
				if ( (*(THIS->constraintProduct))( ii,delta_x, &(delta_Ax[ii]) ) != 0 )
					return THROWERROR( RET_ERROR_IN_CONSTRAINTPRODUCT );
			}
		}
	}*/

	if ( Constraints_hasNoLower( &(THIS->constraints) ) == BT_FALSE )
	{
		for( i=0; i<nIAC; ++i )
		{
			ii = IAC_idx[i];
			num[i] = qpOASES_getMax( THIS->Ax_l[ii],0.0 );
			den[i] = delta_lbA[ii] - delta_Ax[ii];
		}

		QProblem_performRatioTestC( THIS,nIAC,IAC_idx,&(THIS->constraints), num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			*BC_idx = BC_idx_tmp;
			*BC_status = ST_LOWER;
			*BC_isBound = BT_FALSE;
		}
	}

	if ( Constraints_hasNoUpper( &(THIS->constraints) ) == BT_FALSE )
	{
		for( i=0; i<nIAC; ++i )
		{
			ii = IAC_idx[i];
			num[i] = qpOASES_getMax( THIS->Ax_u[ii],0.0 );
			den[i] = delta_Ax[ii] - delta_ubA[ii];
		}

		QProblem_performRatioTestC( THIS,nIAC,IAC_idx,&(THIS->constraints), num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			*BC_idx = BC_idx_tmp;
			*BC_status = ST_UPPER;
			*BC_isBound = BT_FALSE;
		}
	}


	for( i=0; i<nIAC; ++i )
	{
		ii = IAC_idx[i];

		if ( Constraints_getType( &(THIS->constraints),ii ) != ST_UNBOUNDED )
		{
			delta_Ax_l[ii] = delta_Ax[ii] - delta_lbA[ii];
			delta_Ax_u[ii] = delta_ubA[ii] - delta_Ax[ii];
		}
	}


	/* 2) Ensure that inactive bounds remain valid
	 *    (ignoring unbounded variables). */
	/* inactive lower bounds */
	if ( Bounds_hasNoLower( &(THIS->bounds) ) == BT_FALSE )
	{
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			num[i] = qpOASES_getMax( THIS->x[ii] - THIS->lb[ii],0.0 );
			den[i] = delta_lb[ii] - delta_xFR[i];
		}

		QProblem_performRatioTestB( THIS,nFR,FR_idx,&(THIS->bounds),num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			*BC_idx = BC_idx_tmp;
			*BC_status = ST_LOWER;
			*BC_isBound = BT_TRUE;
		}
	}

	/* inactive upper bounds */
	if ( Bounds_hasNoUpper( &(THIS->bounds) ) == BT_FALSE )
	{
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			num[i] = qpOASES_getMax( THIS->ub[ii] - THIS->x[ii],0.0 );
			den[i] = delta_xFR[i] - delta_ub[ii];
		}

		QProblem_performRatioTestB( THIS,nFR,FR_idx,&(THIS->bounds),num,den, THIS->options.epsNum,THIS->options.epsDen, &(THIS->tau),&BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			*BC_idx = BC_idx_tmp;
			*BC_status = ST_UPPER;
			*BC_isBound = BT_TRUE;
		}
	}


	#ifndef __XPCTARGET__
	if ( *BC_status == ST_UNDEFINED )
		snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"Stepsize is %.15e!",THIS->tau );
	else
		snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"Stepsize is %.15e! (idx = %d, isBound = %d, status = %d)",THIS->tau,*BC_idx,*BC_isBound,*BC_status );

	MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_STEPSIZE_NONPOSITIVE,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
	#endif


	/* III) PERFORM STEP ALONG HOMOTOPY PATH */
	if ( THIS->tau > QPOASES_ZERO )
	{
		/* 1) Perform step in primal und dual space... */
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

		for( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			THIS->y[nV+ii] += THIS->tau * delta_yAC[i];
		}

		/* 2) Shift QP data. */
		for( i=0; i<nV; ++i )
		{
			THIS->g[i]  += THIS->tau * delta_g[i];
			THIS->lb[i] += THIS->tau * delta_lb[i];
			THIS->ub[i] += THIS->tau * delta_ub[i];
		}

		for( i=0; i<nC; ++i )
		{
			THIS->lbA[i] += THIS->tau * delta_lbA[i];
			THIS->ubA[i] += THIS->tau * delta_ubA[i];
		}

		DenseMatrix_subTimes(THIS->A,Constraints_getActive( &(THIS->constraints)),0, 1, 1.0, THIS->x, nV, 0.0, THIS->Ax, nC, BT_FALSE );
		for( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			THIS->Ax_u[ii] = THIS->ubA[ii] - THIS->Ax[ii];
			THIS->Ax_l[ii] = THIS->Ax[ii] - THIS->lbA[ii];
		}
		for( i=0; i<nIAC; ++i )
		{
			ii = IAC_idx[i];
			if ( Constraints_getType( &(THIS->constraints),ii ) != ST_UNBOUNDED )
			{
				THIS->Ax[ii]   += THIS->tau * delta_Ax[ii];
				THIS->Ax_l[ii] += THIS->tau * delta_Ax_l[ii];
				THIS->Ax_u[ii] += THIS->tau * delta_Ax_u[ii];
			}
		}
	}
	else
	{
		/* print a stepsize warning if stepsize is zero */
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
returnValue QProblem_changeActiveSet(	QProblem* THIS, 
										int BC_idx, SubjectToStatus BC_status, BooleanType BC_isBound
										)
{
	int nV = QProblem_getNV( THIS );
	myStatic char messageString[QPOASES_MAX_STRING_LENGTH];
	returnValue returnvalue;

	switch ( BC_status )
	{
		/* Optimal solution found as no working set change detected. */
		case ST_UNDEFINED:
			return SUCCESSFUL_RETURN;

		/* Remove one variable from active set. */
		case ST_INACTIVE:
			if ( BC_isBound == BT_TRUE )
			{
				#ifndef __XPCTARGET__
				snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"bound no. %d.", BC_idx );
				MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
				#endif

				if ( QProblem_removeBound( THIS,BC_idx,BT_TRUE,BT_TRUE,THIS->options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_REMOVE_FROM_ACTIVESET_FAILED );

				THIS->y[BC_idx] = 0.0;
			}
			else
			{
				#ifndef __XPCTARGET__
				snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"constraint no. %d.", BC_idx );
				MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
				#endif

				if ( QProblem_removeConstraint( THIS,BC_idx,BT_TRUE,BT_TRUE,THIS->options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_REMOVE_FROM_ACTIVESET_FAILED );

				THIS->y[nV+BC_idx] = 0.0;
			}
			break;


		/* Add one variable to active set. */
		default:
			if ( BC_isBound == BT_TRUE )
			{
				#ifndef __XPCTARGET__
				if ( BC_status == ST_LOWER )
					snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"lower bound no. %d.", BC_idx );
				else
					snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"upper bound no. %d.", BC_idx );
				MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_ADD_TO_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
				#endif

				returnvalue = QProblem_addBound( THIS,BC_idx,BC_status,BT_TRUE,BT_TRUE );
				if ( returnvalue == RET_ADDBOUND_FAILED_INFEASIBILITY )
					return returnvalue;
				if ( returnvalue != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ADD_TO_ACTIVESET_FAILED );
			}
			else
			{
				#ifndef __XPCTARGET__
				if ( BC_status == ST_LOWER )
					snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"lower constraint's bound no. %d.", BC_idx );
				else
					snprintf( messageString,QPOASES_MAX_STRING_LENGTH,"upper constraint's bound no. %d.", BC_idx );
				MessageHandling_throwInfo( qpOASES_getGlobalMessageHandler(),RET_ADD_TO_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
				#endif

				returnvalue = QProblem_addConstraint( THIS,BC_idx,BC_status,BT_TRUE,BT_TRUE );
				if ( returnvalue == RET_ADDCONSTRAINT_FAILED_INFEASIBILITY )
					return returnvalue;
				if ( returnvalue != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ADD_TO_ACTIVESET_FAILED );
			}
	}

	return SUCCESSFUL_RETURN;
}



/*
 * g e t R e l a t i v e H o m o t o p y L e n g t h
 */
real_t QProblem_getRelativeHomotopyLength(	QProblem* THIS,
											const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new,
											const real_t* const lbA_new, const real_t* const ubA_new
											)
{
	int i;
	int nC = QProblem_getNC( THIS );
	real_t len = QProblemBCPY_getRelativeHomotopyLength( THIS, g_new,lb_new,ub_new );
	real_t d, s;

	/* fprintf( stdFile, "len in homotopyLength = %.3e\n",len ); */

	/* lower constraint bounds */
	for (i = 0; i < nC && lbA_new; i++)
	{
		s = qpOASES_getAbs(lbA_new[i]);
		if (s < 1.0) s = 1.0;
		d = qpOASES_getAbs(lbA_new[i] - THIS->lbA[i]) / s;
		if (d > len) len = d;
	}
	/* fprintf( stdFile, "len in homotopyLength = %.3e\n",len ); */

	/* upper constraint bounds */
	for (i = 0; i < nC && ubA_new; i++)
	{
		s = qpOASES_getAbs(ubA_new[i]);
		if (s < 1.0) s = 1.0;
		d = qpOASES_getAbs(ubA_new[i] - THIS->ubA[i]) / s;
		if (d > len) len = d;
	}
	/* fprintf( stdFile, "len in homotopyLength = %.3e\n",len ); */

	return len;
}


/*
 * u p d a t e F a r B o u n d s
 */
returnValue QProblem_updateFarBounds(	QProblem* THIS,
										real_t curFarBound, int nRamp,
                                        const real_t* const lb_new, real_t* const lb_new_far,
                                        const real_t* const ub_new, real_t* const ub_new_far,
                                        const real_t* const lbA_new, real_t* const lbA_new_far,
                                        const real_t* const ubA_new, real_t* const ubA_new_far
                                        )
{
	int i;
	real_t rampVal, t;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

    returnValue returnvalue = QProblemBCPY_updateFarBounds(	THIS,curFarBound,nRamp,
                                                        	lb_new,lb_new_far, ub_new,ub_new_far
                                                            );
    if ( returnvalue != SUCCESSFUL_RETURN )
        return returnvalue;

	if ( THIS->options.enableRamping == BT_TRUE )
	{
		for ( i=0; i<nC; ++i )
		{
			t = (real_t)((nV+i + THIS->rampOffset) % nRamp) / (real_t)(nRamp-1);
			rampVal = curFarBound * (1.0 + (1.0-t)*THIS->ramp0 + t*THIS->ramp1);

			if ( lbA_new == 0 )
				lbA_new_far[i] = -rampVal;
			else
				lbA_new_far[i] = qpOASES_getMax( -rampVal,lbA_new[i] );

			if ( ubA_new == 0 )
				ubA_new_far[i] = rampVal;
			else
				ubA_new_far[i] = qpOASES_getMin( rampVal,ubA_new[i] );
		}
	}
	else
	{
		for ( i=0; i<nC; ++i )
		{
			if ( lbA_new == 0 )
				lbA_new_far[i] = -curFarBound;
			else
				lbA_new_far[i] = qpOASES_getMax( -curFarBound,lbA_new[i] );

			if ( ubA_new == 0 )
				ubA_new_far[i] = curFarBound;
			else
				ubA_new_far[i] = qpOASES_getMin( curFarBound,ubA_new[i] );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 * u p d a t e F a r B o u n d s
 */
returnValue QProblemBCPY_updateFarBounds(	QProblem* THIS,
											real_t curFarBound, int nRamp,
											const real_t* const lb_new, real_t* const lb_new_far,
											const real_t* const ub_new, real_t* const ub_new_far
											)
{
	int i;
	real_t rampVal, t;
	int nV = QProblem_getNV( THIS );

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
 *  p e r f o r m R a m p i n g
 */
returnValue QProblem_performRamping( QProblem* THIS )
{
	int nV = QProblem_getNV( THIS ), nC = QProblem_getNC( THIS ), bstat, cstat, i, nRamp;
	real_t tP, rampValP, tD, rampValD, sca;

	/* compute number of values in ramp */
	nRamp = nV + nC + nC + nV;
	
	/* ramp inactive variable bounds and active dual bound variables */
	for (i = 0; i < nV; i++)
	{
		switch (Bounds_getType(&(THIS->bounds),i))
		{
			case ST_EQUALITY: 
				THIS->lb[i] = THIS->x[i]; THIS->ub[i] = THIS->x[i];  /* reestablish exact feasibility */
				continue; 

			case ST_BOUNDED:
				tP = (real_t)((i+THIS->rampOffset) % nRamp) / (real_t)(nRamp-1);
				rampValP = (1.0-tP) * THIS->ramp0 + tP * THIS->ramp1;
				tD = (real_t)((nV+nC+nC+i+THIS->rampOffset) % nRamp) / (real_t)(nRamp-1);
				rampValD = (1.0-tD) * THIS->ramp0 + tD * THIS->ramp1;
				bstat = Bounds_getStatus( &(THIS->bounds),i);
				if (bstat != ST_LOWER) { sca = qpOASES_getMax(qpOASES_getAbs(THIS->x[i]), 1.0); THIS->lb[i] = THIS->x[i] - sca * rampValP; }
				if (bstat != ST_UPPER) { sca = qpOASES_getMax(qpOASES_getAbs(THIS->x[i]), 1.0); THIS->ub[i] = THIS->x[i] + sca * rampValP; }
				if (bstat == ST_LOWER) { THIS->lb[i] = THIS->x[i]; THIS->y[i] = +rampValD; }
				if (bstat == ST_UPPER) { THIS->ub[i] = THIS->x[i]; THIS->y[i] = -rampValD; }
				if (bstat == ST_INACTIVE) THIS->y[i] = 0.0; /* reestablish exact complementarity */
				break;
				
			case ST_UNBOUNDED: 
			case ST_DISABLED:
			default: 
				 continue;
		}
	}

	/* ramp inactive constraints and active dual constraint variables */
	for (i = 0; i < nC; i++)
	{
		switch (Constraints_getType( &(THIS->constraints),i))
		{
			case ST_EQUALITY: 
				THIS->lbA[i] = THIS->Ax[i]; THIS->ubA[i] = THIS->Ax[i];  /* reestablish exact feasibility */
				continue; 
							
			case ST_BOUNDED:
				tP = (real_t)((nV+i+THIS->rampOffset) % nRamp) / (real_t)(nRamp-1);
				rampValP = (1.0-tP) * THIS->ramp0 + tP * THIS->ramp1;
				tD = (real_t)((nV+nC+i+THIS->rampOffset) % nRamp) / (real_t)(nRamp-1);
				rampValD = (1.0-tD) * THIS->ramp0 + tD * THIS->ramp1;
				cstat = Constraints_getStatus( &(THIS->constraints),i);
				if (cstat != ST_LOWER) { sca = qpOASES_getMax(qpOASES_getAbs(THIS->Ax[i]), 1.0); THIS->lbA[i] = THIS->Ax[i] - sca * rampValP; }
				if (cstat != ST_UPPER) { sca = qpOASES_getMax(qpOASES_getAbs(THIS->Ax[i]), 1.0); THIS->ubA[i] = THIS->Ax[i] + sca * rampValP; }
				if (cstat == ST_LOWER) { THIS->lbA[i] = THIS->Ax[i]; THIS->y[nV+i] = +rampValD; }
				if (cstat == ST_UPPER) { THIS->ubA[i] = THIS->Ax[i]; THIS->y[nV+i] = -rampValD; }
				if (cstat == ST_INACTIVE) THIS->y[nV+i] = 0.0; /* reestablish exact complementarity */
					
				THIS->Ax_l[i] = THIS->Ax[i] - THIS->lbA[i];
				THIS->Ax_u[i] = THIS->ubA[i] - THIS->Ax[i];
				break;

			case ST_UNBOUNDED: 
			case ST_DISABLED: 
			default: 
				continue;
		}
	}

	/* reestablish exact stationarity */
	QProblem_setupAuxiliaryQPgradient( THIS );

	/* advance ramp offset to avoid Ramping cycles */
	(THIS->rampOffset)++;

	return SUCCESSFUL_RETURN;
}



/*
 *	p e r f o r m R a t i o T e s t
 */
returnValue QProblem_performRatioTestC(	QProblem* THIS, 
										const int nIdx,
										const int* const idxList,
										Constraints* const subjectTo,
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

		if ( Constraints_getType( subjectTo,ii ) != ST_EQUALITY )
		{
			if ( ( Constraints_getStatus( subjectTo,ii ) == ST_LOWER ) || ( Constraints_getStatus( subjectTo,ii ) == ST_INACTIVE ) )
			{
				if ( QProblem_isBlocking( THIS,num[i],den[i],epsNum,epsDen,t ) == BT_TRUE )
				{
					*t = num[i] / den[i];
					*BC_idx = ii;
				}
			}
			else
			if ( Constraints_getStatus( subjectTo,ii ) == ST_UPPER )
			{
				if ( QProblem_isBlocking( THIS,-num[i],-den[i],epsNum,epsDen,t ) == BT_TRUE )
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
 * p e r f o r m D r i f t C o r r e c t i o n
 */
returnValue QProblem_performDriftCorrection( QProblem* THIS )
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	for ( i=0; i<nV; ++i )
	{
		switch ( Bounds_getType( &(THIS->bounds),i ) )
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

	for ( i=0; i<nC; ++i )
	{
		switch ( Constraints_getType( &(THIS->constraints),i ) )
		{
			case ST_BOUNDED:
				switch ( Constraints_getStatus( &(THIS->constraints),i ) )
				{
					case ST_LOWER:
						THIS->lbA[i] = THIS->Ax[i];
						THIS->Ax_l[i] = 0.0;
						THIS->ubA[i] = qpOASES_getMax (THIS->ubA[i], THIS->Ax[i]);
						THIS->Ax_u[i] = THIS->ubA[i] - THIS->Ax[i];
						THIS->y[i+nV] = qpOASES_getMax (THIS->y[i+nV], 0.0);
						break;
					case ST_UPPER:
						THIS->lbA[i] = qpOASES_getMin (THIS->lbA[i], THIS->Ax[i]);
						THIS->Ax_l[i] = THIS->Ax[i] - THIS->lbA[i];
						THIS->ubA[i] = THIS->Ax[i];
						THIS->Ax_u[i] = 0.0;
						THIS->y[i+nV] = qpOASES_getMin (THIS->y[i+nV], 0.0);
						break;
					case ST_INACTIVE:
						THIS->lbA[i] = qpOASES_getMin (THIS->lbA[i], THIS->Ax[i]);
						THIS->Ax_l[i] = THIS->Ax[i] - THIS->lbA[i];
						THIS->ubA[i] = qpOASES_getMax (THIS->ubA[i], THIS->Ax[i]);
						THIS->Ax_u[i] = THIS->ubA[i] - THIS->Ax[i];
						THIS->y[i+nV] = 0.0;
						break;
					case ST_UNDEFINED:
					case ST_INFEASIBLE_LOWER:
					case ST_INFEASIBLE_UPPER:
						break;
				}
				break;
			case ST_EQUALITY:
				THIS->lbA[i] = THIS->Ax[i];
				THIS->Ax_l[i] = 0.0;
				THIS->ubA[i] = THIS->Ax[i];
				THIS->Ax_u[i] = 0.0;
				break;
			case ST_UNBOUNDED:
			case ST_UNKNOWN:
            case ST_DISABLED:
				break;
		}
	}

	return QProblem_setupAuxiliaryQPgradient( THIS );
}

/*
 *	s e t u p A u x i l i a r y Q P
 */
returnValue QProblem_setupAuxiliaryQP( QProblem* THIS, Bounds* const guessedBounds, Constraints* const guessedConstraints )
{
	int i, j;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	/* consistency check */
	if ( ( guessedBounds == 0 ) || ( guessedConstraints == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* nothing to do */
	if ( ( guessedBounds == &(THIS->bounds) ) && ( guessedConstraints == &(THIS->constraints) ) )
		return SUCCESSFUL_RETURN;

	THIS->status = QPS_PREPARINGAUXILIARYQP;


	/* I) SETUP WORKING SET ... */
	if ( QProblem_shallRefactorise( THIS,guessedBounds,guessedConstraints ) == BT_TRUE )
	{
		/* ... WITH REFACTORISATION: */
		/* 1) Reset bounds/constraints ... */
		Bounds_init( &(THIS->bounds),nV );
		Constraints_init( &(THIS->constraints),nC );

		/*    ... and set them up afresh. */
		if ( QProblem_setupSubjectToType( THIS ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( Bounds_setupAllFree( &(THIS->bounds) ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( Constraints_setupAllInactive( &(THIS->constraints) ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 2) Setup TQ factorisation. */
		if ( QProblem_setupTQfactorisation( THIS ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 3) Setup guessed working sets afresh. */
		if ( QProblem_setupAuxiliaryWorkingSet( THIS,guessedBounds,guessedConstraints,BT_TRUE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 4) Computes Cholesky decomposition of projected Hessian
		 *    This now handles all special cases (no active bounds/constraints, no nullspace) */
		if ( QProblem_computeProjectedCholesky( THIS ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}
	else
	{
		/* ... WITHOUT REFACTORISATION: */
		if ( QProblem_setupAuxiliaryWorkingSet( THIS,guessedBounds,guessedConstraints,BT_FALSE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}


	/* II) SETUP AUXILIARY QP DATA: */
	/* 1) Ensure that dual variable is zero for free bounds and inactive constraints. */
	for ( i=0; i<nV; ++i )
		if ( Bounds_getStatus( &(THIS->bounds),i ) == ST_INACTIVE )
			THIS->y[i] = 0.0;

	for ( i=0; i<nC; ++i )
		if ( Constraints_getStatus( &(THIS->constraints),i ) == ST_INACTIVE )
			THIS->y[nV+i] = 0.0;

	/* 2) Setup gradient and (constraints') bound vectors. */
	if ( QProblem_setupAuxiliaryQPgradient( THIS ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	DenseMatrix_times(THIS->A,1, 1.0, THIS->x, nV, 0.0, THIS->Ax, nC);
	for ( j=0; j<nC; ++j )
	{
		THIS->Ax_l[j] = THIS->Ax[j];
		THIS->Ax_u[j] = THIS->Ax[j];
	}

	/* (also sets Ax_l and Ax_u) */
	if ( QProblem_setupAuxiliaryQPbounds( THIS,0,0,BT_FALSE ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	s h a l l R e f a c t o r i s e
 */

BooleanType QProblem_shallRefactorise(	QProblem* THIS,
										Bounds* const guessedBounds,
										Constraints* const guessedConstraints
										)
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	int differenceNumberBounds = 0;
	int differenceNumberConstraints = 0;

	/* always refactorise if Hessian is not known to be positive definite */
	if ( ( THIS->hessianType == HST_SEMIDEF ) || ( THIS->hessianType == HST_INDEF ) )
		return BT_TRUE;

	/* 1) Determine number of bounds that have same status
	 *    in guessed AND current bounds.*/
	for( i=0; i<nV; ++i )
		if ( Bounds_getStatus( guessedBounds,i ) != Bounds_getStatus( &(THIS->bounds),i ) )
			++differenceNumberBounds;

	/* 2) Determine number of constraints that have same status
	 *    in guessed AND current constraints.*/
	for( i=0; i<nC; ++i )
		if ( Constraints_getStatus( guessedConstraints,i ) != Constraints_getStatus( &(THIS->constraints),i ) )
			++differenceNumberConstraints;

	/* 3) Decide wheter to refactorise or not. */
	if ( 2*(differenceNumberBounds+differenceNumberConstraints) > Constraints_getNAC( guessedConstraints )+Bounds_getNFX( guessedBounds ) )
		return BT_TRUE;
	else
		return BT_FALSE;
}


/*
 *	s e t u p Q P d a t a
 */
returnValue QProblem_setupQPdataM(	QProblem* THIS,
									DenseMatrix *_H, const real_t* const _g, DenseMatrix *_A,
									const real_t* const _lb, const real_t* const _ub,
									const real_t* const _lbA, const real_t* const _ubA
									)
{
	if ( _H == 0 )
	{
		if ( _A == 0 )
			return QProblem_setupQPdata( THIS,(real_t*)0,_g,(real_t*)0,_lb,_ub,_lbA,_ubA );
		else
			return QProblem_setupQPdata( THIS,(real_t*)0,_g,DenseMatrix_getVal(_A),_lb,_ub,_lbA,_ubA );
	}
	else
	{
		if ( _A == 0 )
			return QProblem_setupQPdata( THIS,DenseMatrix_getVal(_H),_g,(real_t*)0,_lb,_ub,_lbA,_ubA );
		else
			return QProblem_setupQPdata( THIS,DenseMatrix_getVal(_H),_g,DenseMatrix_getVal(_A),_lb,_ub,_lbA,_ubA );
	}
}


/*
 *	s e t u p Q P d a t a
 */
returnValue QProblem_setupQPdata(	QProblem* THIS,
									real_t* const _H, const real_t* const _g, real_t* const _A,
									const real_t* const _lb, const real_t* const _ub,
									const real_t* const _lbA, const real_t* const _ubA
									)
{
	int i;
	int nC = QProblem_getNC( THIS );


	/* 1) Load Hessian matrix as well as lower and upper bounds vectors. */
	if ( QProblemBCPY_setupQPdata( THIS,_H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Load constraint matrix. */
	if ( ( nC > 0 ) && ( _A == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( nC > 0 )
	{
		QProblem_setA( THIS,_A );

		/* 3) Setup lower constraints' bounds vector. */
		if ( _lbA != 0 )
		{
			QProblem_setLBA( THIS,_lbA );
		}
		else
		{
			/* if no lower constraints' bounds are specified, set them to -infinity */
			for( i=0; i<nC; ++i )
				THIS->lbA[i] = -QPOASES_INFTY;
		}

		/* 4) Setup upper constraints' bounds vector. */
		if ( _ubA != 0 )
		{
			QProblem_setUBA( THIS,_ubA );
		}
		else
		{
			/* if no upper constraints' bounds are specified, set them to infinity */
			for( i=0; i<nC; ++i )
				THIS->ubA[i] = QPOASES_INFTY;
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p Q P d a t a F r o m F i l e
 */
returnValue QProblem_setupQPdataFromFile(	QProblem* THIS,
											const char* const H_file, const char* const g_file, const char* const A_file,
											const char* const lb_file, const char* const ub_file,
											const char* const lbA_file, const char* const ubA_file
											)
{
	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );

	returnValue returnvalue;

	myStatic real_t _A[NCMAX*NVMAX];


	/* 1) Load Hessian matrix as well as lower and upper bounds vectors from files. */
	returnvalue = QProblemBCPY_setupQPdataFromFile( THIS,H_file,g_file,lb_file,ub_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
		return THROWERROR( returnvalue );

	/* 2) Load constraint matrix from file. */
	if ( ( nC > 0 ) && ( A_file == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( nC > 0 )
	{
		returnvalue = qpOASES_readFromFileM( _A, nC,nV, A_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
		QProblem_setA( THIS,_A );

		/* 3) Load lower constraints' bounds vector from file. */
		if ( lbA_file != 0 )
		{
			returnvalue = qpOASES_readFromFileV( THIS->lbA, nC, lbA_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* if no lower constraints' bounds are specified, set them to -infinity */
			for( i=0; i<nC; ++i )
				THIS->lbA[i] = -QPOASES_INFTY;
		}

		/* 4) Load upper constraints' bounds vector from file. */
		if ( ubA_file != 0 )
		{
			returnvalue = qpOASES_readFromFileV( THIS->ubA, nC, ubA_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* if no upper constraints' bounds are specified, set them to infinity */
			for( i=0; i<nC; ++i )
				THIS->ubA[i] = QPOASES_INFTY;
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	l o a d Q P v e c t o r s F r o m F i l e
 */
returnValue QProblem_loadQPvectorsFromFile(	QProblem* THIS,
											const char* const g_file, const char* const lb_file, const char* const ub_file,
											const char* const lbA_file, const char* const ubA_file,
											real_t* const g_new, real_t* const lb_new, real_t* const ub_new,
											real_t* const lbA_new, real_t* const ubA_new
											)
{
	int nC = QProblem_getNC( THIS );

	returnValue returnvalue;


	/* 1) Load gradient vector as well as lower and upper bounds vectors from files. */
	returnvalue = QProblemBCPY_loadQPvectorsFromFile( THIS,g_file,lb_file,ub_file, g_new,lb_new,ub_new );
	if ( returnvalue != SUCCESSFUL_RETURN )
		return THROWERROR( returnvalue );

	if ( nC > 0 )
	{
		/* 2) Load lower constraints' bounds vector from file. */
		if ( lbA_file != 0 )
		{
			if ( lbA_new != 0 )
			{
				returnvalue = qpOASES_readFromFileV( lbA_new, nC, lbA_file );
				if ( returnvalue != SUCCESSFUL_RETURN )
					return THROWERROR( returnvalue );
			}
			else
			{
				/* If filename is given, storage must be provided! */
				return THROWERROR( RET_INVALID_ARGUMENTS );
			}
		}

		/* 3) Load upper constraints' bounds vector from file. */
		if ( ubA_file != 0 )
		{
			if ( ubA_new != 0 )
			{
				returnvalue = qpOASES_readFromFileV( ubA_new, nC, ubA_file );
				if ( returnvalue != SUCCESSFUL_RETURN )
					return THROWERROR( returnvalue );
			}
			else
			{
				/* If filename is given, storage must be provided! */
				return THROWERROR( RET_INVALID_ARGUMENTS );
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t I t e r a t i o n
 */
returnValue QProblem_printIteration( 	QProblem* THIS, 
										int iter,
										int BC_idx,	SubjectToStatus BC_status, BooleanType BC_isBound,
										real_t homotopyLength
		  								)
{
	#ifndef __XPCTARGET__

	int i;
	int nV = QProblem_getNV( THIS );
	int nC = QProblem_getNC( THIS );
	int nAC = QProblem_getNAC( THIS );

	real_t stat, bfeas, cfeas, bcmpl, ccmpl, Tmaxomin;
	real_t Tmin, Tmax;

	myStatic real_t grad[NVMAX];
	myStatic real_t AX[NCMAX];

	myStatic char myPrintfString[QPOASES_MAX_STRING_LENGTH];
	myStatic char info[QPOASES_MAX_STRING_LENGTH];
	const char excStr[] = " ef";

	/* consistency check */
	if ( iter < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	switch ( THIS->options.printLevel )
	{
		case PL_DEBUG_ITER:
			stat = bfeas = cfeas = bcmpl = ccmpl = Tmaxomin = 0.0;

			/* stationarity */
			for (i = 0; i < nV; i++) grad[i] = THIS->g[i] - THIS->y[i];
			DenseMatrix_times(THIS->H,1, 1.0, THIS->x, nV, 1.0, grad, nV);
			DenseMatrix_transTimes(THIS->A,1, -1.0, THIS->y+nV, nC, 1.0, grad, nV);
			for (i = 0; i < nV; i++) if (qpOASES_getAbs(grad[i]) > stat) stat = qpOASES_getAbs(grad[i]);

			/* feasibility */
			for (i = 0; i < nV; i++) if (THIS->lb[i] - THIS->x[i] > bfeas) bfeas = THIS->lb[i] - THIS->x[i];
			for (i = 0; i < nV; i++) if (THIS->x[i] - THIS->ub[i] > bfeas) bfeas = THIS->x[i] - THIS->ub[i];
			DenseMatrix_times(THIS->A,1, 1.0, THIS->x, nV, 0.0, AX, nC);
			for (i = 0; i < nC; i++) if (THIS->lbA[i] - AX[i] > cfeas) cfeas = THIS->lbA[i] - AX[i];
			for (i = 0; i < nC; i++) if (AX[i] - THIS->ubA[i] > cfeas) cfeas = AX[i] - THIS->ubA[i];

			/* complementarity */
			for (i = 0; i < nV; i++) if (THIS->y[i] > +QPOASES_EPS && qpOASES_getAbs((THIS->lb[i] - THIS->x[i])*THIS->y[i]) > bcmpl) bcmpl = qpOASES_getAbs((THIS->lb[i] - THIS->x[i])*THIS->y[i]);
			for (i = 0; i < nV; i++) if (THIS->y[i] < -QPOASES_EPS && qpOASES_getAbs((THIS->ub[i] - THIS->x[i])*THIS->y[i]) > bcmpl) bcmpl = qpOASES_getAbs((THIS->ub[i] - THIS->x[i])*THIS->y[i]);
			for (i = 0; i < nC; i++) if (THIS->y[nV+i] > +QPOASES_EPS && qpOASES_getAbs((THIS->lbA[i]-AX[i])*THIS->y[nV+i]) > ccmpl) ccmpl = qpOASES_getAbs((THIS->lbA[i]-AX[i])*THIS->y[nV+i]);
			for (i = 0; i < nC; i++) if (THIS->y[nV+i] < -QPOASES_EPS && qpOASES_getAbs((THIS->ubA[i]-AX[i])*THIS->y[nV+i]) > ccmpl) ccmpl = qpOASES_getAbs((THIS->ubA[i]-AX[i])*THIS->y[nV+i]);

			Tmin = 1.0e16; Tmax = 0.0;
			for (i = 0; i < nAC; i++)
				if (qpOASES_getAbs(TT(i,THIS->sizeT-i-1)) < Tmin)
					Tmin = qpOASES_getAbs(TT(i,THIS->sizeT-i-1));
				else if (qpOASES_getAbs(TT(i,THIS->sizeT-i-1)) > Tmax)
					Tmax = qpOASES_getAbs(TT(i,THIS->sizeT-i-1));
			Tmaxomin = Tmax/Tmin;

			if (iter % 10 == 0)
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "\n%5s %4s %4s %4s %4s %9s %9s %9s %9s %9s %9s %9s %9s\n",
						"iter", "addB", "remB", "addC", "remC", "hom len", "tau", "stat",
						"bfeas", "cfeas", "bcmpl", "ccmpl", "Tmin");
				qpOASES_myPrintf( myPrintfString );
			}

			snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%5d ", iter);
			qpOASES_myPrintf( myPrintfString );

			if (THIS->tabularOutput.idxAddB >= 0)
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%4d ", THIS->tabularOutput.idxAddB);
				qpOASES_myPrintf( myPrintfString );
			}
			else
			{
				qpOASES_myPrintf( "     " );
			}

			if (THIS->tabularOutput.idxRemB >= 0)
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%4d ", THIS->tabularOutput.idxRemB);
				qpOASES_myPrintf( myPrintfString );
			}
			else
			{
				qpOASES_myPrintf( "     " );
			}

			if (THIS->tabularOutput.idxAddC >= 0)
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%4d ", THIS->tabularOutput.idxAddC);
				qpOASES_myPrintf( myPrintfString );
			}
			else
			{
				qpOASES_myPrintf( "     " );
			}

			if (THIS->tabularOutput.idxRemC >= 0)
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%4d ", THIS->tabularOutput.idxRemC);
				qpOASES_myPrintf( myPrintfString );
			}
			else
			{
				qpOASES_myPrintf( "     " );
			}

			snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
					homotopyLength, THIS->tau, stat, bfeas, cfeas, bcmpl, ccmpl, Tmin);
			qpOASES_myPrintf( myPrintfString );
			break;

		case PL_TABULAR:
			if (iter % 10 == 0)
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "\n%5s %6s %6s %6s %6s %9s %9s\n",
						"iter", "addB", "remB", "addC", "remC", "hom len", "tau" );
				qpOASES_myPrintf( myPrintfString );
			}
			snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%5d ",iter);
			qpOASES_myPrintf( myPrintfString );

			if (THIS->tabularOutput.idxAddB >= 0)
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%5d%c ", THIS->tabularOutput.idxAddB, excStr[THIS->tabularOutput.excAddB]);
				qpOASES_myPrintf( myPrintfString );
			}
			else
			{
				qpOASES_myPrintf( "       " );
			}

			if (THIS->tabularOutput.idxRemB >= 0)
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%5d%c ", THIS->tabularOutput.idxRemB, excStr[THIS->tabularOutput.excRemB]);
				qpOASES_myPrintf( myPrintfString );
			}
			else 
			{
				qpOASES_myPrintf( "       " );
			}

			if (THIS->tabularOutput.idxAddC >= 0) 
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%5d%c ", THIS->tabularOutput.idxAddC, excStr[THIS->tabularOutput.excAddC]);
				qpOASES_myPrintf( myPrintfString );
			}
			else 
			{
				qpOASES_myPrintf( "       " );
			}

			if (THIS->tabularOutput.idxRemC >= 0) 
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%5d%c ", THIS->tabularOutput.idxRemC, excStr[THIS->tabularOutput.excRemC]);
				qpOASES_myPrintf( myPrintfString );
			}
			else 
			{
				qpOASES_myPrintf( "       " );
			}

			snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH, "%9.2e %9.2e\n", homotopyLength, THIS->tau);
			qpOASES_myPrintf( myPrintfString );
			break;

		case PL_MEDIUM:
			/* 1) Print header at first iteration. */
 			if ( iter == 0 )
			{
				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"\n\n####################   qpOASES  --  QP NO. %3.0d   #####################\n\n", THIS->count );
				qpOASES_myPrintf( myPrintfString );

				qpOASES_myPrintf( "    Iter   |    StepLength    |       Info       |   nFX   |   nAC    \n" );
				qpOASES_myPrintf( " ----------+------------------+------------------+---------+--------- \n" );
			}

			/* 2) Print iteration line. */
			if ( BC_status == ST_UNDEFINED )
			{
				if ( THIS->hessianType == HST_ZERO )
					snprintf( info,3,"LP" );
				else
					snprintf( info,3,"QP" );

				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |    %s SOLVED     |  %4.1d   |  %4.1d   \n", iter,THIS->tau,info,QProblem_getNFX( THIS ),QProblem_getNAC( THIS ) );
				qpOASES_myPrintf( myPrintfString );
			}
			else
			{
				if ( BC_status == ST_INACTIVE )
					snprintf( info,5,"REM " );
				else
					snprintf( info,5,"ADD " );

				if ( BC_isBound == BT_TRUE )
					snprintf( &(info[4]),4,"BND" );
				else
					snprintf( &(info[4]),4,"CON" );

				snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |   %s %4.1d   |  %4.1d   |  %4.1d   \n", iter,THIS->tau,info,BC_idx,QProblem_getNFX( THIS ),QProblem_getNAC( THIS ) );
				qpOASES_myPrintf( myPrintfString );
			}
			break;

		default:
			/* nothing to display */
			break;
	}
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	a r e B o u n d s C o n s i s t e n t
 */
returnValue QProblem_areBoundsConsistent(	QProblem* THIS,
											const real_t* const lb_new, const real_t* const ub_new,
											const real_t* const lbA_new, const real_t* const ubA_new
											)
{
	int i;

	if (lb_new && ub_new) {
		for (i = 0; i < QProblem_getNV(THIS); ++i) {
			if (lb_new[i] > ub_new[i]+QPOASES_EPS) {
				return RET_QP_INFEASIBLE;
			}
		}
	}

	if (lbA_new && ubA_new) {
		for (i = 0; i < QProblem_getNC(THIS); ++i) {
			if (lbA_new[i] > ubA_new[i]+QPOASES_EPS) {
				return RET_QP_INFEASIBLE;
			}
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 * d r o p I n f e a s i b l e s
 */
returnValue QProblem_dropInfeasibles(	QProblem* THIS,
										int BC_number, SubjectToStatus BC_status, BooleanType BC_isBound,
										real_t *xiB, real_t *xiC 
										)
{
	int i;
	
	int nAC                   = QProblem_getNAC( THIS );
	int nFX                   = QProblem_getNFX( THIS );	
	int blockingPriority      = (BC_isBound) ? THIS->options.dropBoundPriority : THIS->options.dropIneqConPriority;
	int y_min_number          = -1;
	BooleanType y_min_isBound = BC_isBound;
	int y_min_priority        = blockingPriority;

	int *AC_idx, *FX_idx;

	Indexlist_getNumberArray( Constraints_getActive( &(THIS->constraints) ),&AC_idx );
	Indexlist_getNumberArray( Bounds_getFixed( &(THIS->bounds) ),&FX_idx );
	
	if (THIS->options.dropEqConPriority <= y_min_priority) 
	{
		/* look for an equality constraint we can drop according to priorities */
		for ( i = 0; i < nAC; ++i )
			if ( (Constraints_getType( &(THIS->constraints),i) == ST_EQUALITY)
				&& (qpOASES_getAbs (xiC[i]) > THIS->options.epsDen) ) 
			{
				y_min_number = AC_idx[i];
				y_min_isBound = BT_FALSE;
				y_min_priority = THIS->options.dropEqConPriority;
				break;
			}
	}
	
	if (THIS->options.dropIneqConPriority <= y_min_priority) 
	{
		/* look for an inequality constraint we can drop according to priorities */
		for ( i = 0; i < nAC; ++i )
			if ( (Constraints_getType( &(THIS->constraints),i) == ST_BOUNDED)
				&& (qpOASES_getAbs (xiC[i]) > THIS->options.epsDen) ) 
			{
				y_min_number = AC_idx[i];
				y_min_isBound = BT_FALSE;
				y_min_priority = THIS->options.dropIneqConPriority;
				break;
			}
	}
	
	if (THIS->options.dropBoundPriority <= y_min_priority) 
	{
		/* look for a simple bound we can drop according to priorities */
		for ( i = 0; i < nFX; ++i )
			if (qpOASES_getAbs (xiB[i]) > THIS->options.epsDen) 
			{
				y_min_number = FX_idx[i];
				y_min_isBound = BT_TRUE;
				y_min_priority = THIS->options.dropBoundPriority;
				break;
			}
	}
	
	if (y_min_number >= 0) {
		
		/* drop active equality or active bound we have found */
		if (y_min_isBound) {
			SubjectToStatus status_ = Bounds_getStatus( &(THIS->bounds),y_min_number);
			QProblem_removeBound( THIS,y_min_number, BT_TRUE, BT_FALSE, BT_FALSE);
			Bounds_setStatus(&(THIS->bounds),y_min_number, (status_ == ST_LOWER) ? ST_INFEASIBLE_LOWER : ST_INFEASIBLE_UPPER);
			/* TODO: fix duals THIS->y[] */
			/*fprintf (stdFile, "Dropping bounds %d for %s %d\n", y_min_number, BC_isBound?"bound":"constraint", BC_number);*/
		} else {
			SubjectToStatus status_ = Constraints_getStatus( &(THIS->constraints),y_min_number);
			QProblem_removeConstraint( THIS,y_min_number, BT_TRUE, BT_FALSE, BT_FALSE);
			Constraints_setStatus( &(THIS->constraints),y_min_number, (status_ == ST_LOWER) ? ST_INFEASIBLE_LOWER : ST_INFEASIBLE_UPPER);
			/* TODO: fix duals THIS->y[] */
			/*fprintf (stdFile, "Dropping constraint %d for %s %d\n", y_min_number, BC_isBound?"bound":"constraint", BC_number);*/
		}
		
		/*/ ... now return, add the blocking constraint, and continue solving QP with dropped bound/constraint */
		return SUCCESSFUL_RETURN;
		
	} else {
		
		/* nothing found, then drop the blocking (still inactive) constraint */
		if (BC_isBound)
			Bounds_setStatus(&(THIS->bounds),BC_number, (BC_status == ST_LOWER) ? ST_INFEASIBLE_LOWER : ST_INFEASIBLE_UPPER);
		else
			Constraints_setStatus( &(THIS->constraints),BC_number, (BC_status == ST_LOWER) ? ST_INFEASIBLE_LOWER : ST_INFEASIBLE_UPPER);
		
		/*fprintf (stdFile, "Dropping %s %d itself\n", BC_isBound?"bound":"constraint", BC_number);*/
		
		/* ... now return, and continue solving QP with dropped bound/constraint */
		return RET_ENSURELI_DROPPED;
	}	
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
