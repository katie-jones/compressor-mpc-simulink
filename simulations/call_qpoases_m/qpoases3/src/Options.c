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
 *	\file src/Options.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0embedded
 *	\date 2007-2014
 *
 *	Implementation of the Options class designed to manage working sets of
 *	constraints and bounds within a QProblem.
 */


#include <qpOASES/Options.h>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	O p t i o n s
 */
void OptionsCON(	Options* THIS
					)
{
	#ifdef __CODE_GENERATION__
	Options_setToFast( THIS );
	#else
	Options_setToDefault( THIS );
	#endif /* __CODE_GENERATION__ */
}


/*
 *	c o p y
 */
void OptionsCPY(	Options* FROM,
					Options* TO
					)
{
	TO->printLevel                    =  FROM->printLevel;

	TO->enableRamping                 =  FROM->enableRamping;
	TO->enableFarBounds               =  FROM->enableFarBounds;
	TO->enableFlippingBounds          =  FROM->enableFlippingBounds;
	TO->enableRegularisation          =  FROM->enableRegularisation;
	TO->enableFullLITests             =  FROM->enableFullLITests;
	TO->enableNZCTests                =  FROM->enableNZCTests;
	TO->enableDriftCorrection         =  FROM->enableDriftCorrection;
	TO->enableCholeskyRefactorisation =  FROM->enableCholeskyRefactorisation;
	TO->enableEqualities              =  FROM->enableEqualities;

	TO->terminationTolerance          =  FROM->terminationTolerance;
	TO->boundTolerance                =  FROM->boundTolerance;
	TO->boundRelaxation               =  FROM->boundRelaxation;
	TO->epsNum                        =  FROM->epsNum;
	TO->epsDen                        =  FROM->epsDen;
	TO->maxPrimalJump                 =  FROM->maxPrimalJump;
	TO->maxDualJump                   =  FROM->maxDualJump;

	TO->initialRamping                =  FROM->initialRamping;
	TO->finalRamping                  =  FROM->finalRamping;
	TO->initialFarBounds              =  FROM->initialFarBounds;
	TO->growFarBounds                 =  FROM->growFarBounds;
 	TO->initialStatusBounds           =  FROM->initialStatusBounds;
	TO->epsFlipping                   =  FROM->epsFlipping;
	TO->numRegularisationSteps        =  FROM->numRegularisationSteps;
	TO->epsRegularisation             =  FROM->epsRegularisation;
	TO->numRefinementSteps            =  FROM->numRefinementSteps;
	TO->epsIterRef                    =  FROM->epsIterRef;
	TO->epsLITests                    =  FROM->epsLITests;
	TO->epsNZCTests                   =  FROM->epsNZCTests;

	TO->enableDropInfeasibles         =  FROM->enableDropInfeasibles;
    TO->dropBoundPriority             =  FROM->dropBoundPriority;
    TO->dropEqConPriority             =  FROM->dropEqConPriority;
    TO->dropIneqConPriority           =  FROM->dropIneqConPriority;
}



/*
 *	s e t T o D e f a u l t
 */
returnValue Options_setToDefault(	Options* THIS
									)
{
	THIS->printLevel = PL_MEDIUM;
	#ifdef __DEBUG__
	THIS->printLevel = PL_HIGH;
	#endif
	#ifdef __XPCTARGET__
	THIS->printLevel = PL_NONE;
	#endif
	#ifdef __SUPPRESSANYOUTPUT__
	THIS->printLevel = PL_NONE;
	#endif
	#ifdef __CODE_GENERATION__
	THIS->printLevel = QPOASES_PRINTLEVEL;
	#endif

	THIS->enableRamping                 =  BT_TRUE;
	THIS->enableFarBounds               =  BT_TRUE;
	THIS->enableFlippingBounds          =  BT_TRUE;
	THIS->enableRegularisation          =  BT_FALSE;
	THIS->enableFullLITests             =  BT_FALSE;
	THIS->enableNZCTests                =  BT_TRUE;
	THIS->enableDriftCorrection         =  1;
	THIS->enableCholeskyRefactorisation =  0;
	THIS->enableEqualities              =  BT_FALSE;

	#ifdef __USE_SINGLE_PRECISION__
	THIS->terminationTolerance          =  1.0e2 * QPOASES_EPS;
	THIS->boundTolerance                =  1.0e2 * QPOASES_EPS;
	#else
	THIS->terminationTolerance          =  5.0e6 * QPOASES_EPS;
	THIS->boundTolerance                =  1.0e6 * QPOASES_EPS;
	#endif
	THIS->boundRelaxation               =  1.0e4;
	#ifdef __USE_SINGLE_PRECISION__
	THIS->epsNum                        = -1.0e2 * QPOASES_EPS;
	THIS->epsDen                        =  1.0e2 * QPOASES_EPS;
	#else
	THIS->epsNum                        = -1.0e3 * QPOASES_EPS;
	THIS->epsDen                        =  1.0e3 * QPOASES_EPS;
	#endif
	THIS->maxPrimalJump                 =  1.0e8;
	THIS->maxDualJump                   =  1.0e8;

	THIS->initialRamping                =  0.5;
	THIS->finalRamping                  =  1.0;
	THIS->initialFarBounds              =  1.0e6;
	THIS->growFarBounds                 =  1.0e3;
 	THIS->initialStatusBounds           =  ST_LOWER;
	#ifdef __USE_SINGLE_PRECISION__
	THIS->epsFlipping                   =  5.0e1 * QPOASES_EPS;
	#else
	THIS->epsFlipping                   =  1.0e3 * QPOASES_EPS;
	#endif
	THIS->numRegularisationSteps        =  0;
	#ifdef __USE_SINGLE_PRECISION__
	THIS->epsRegularisation             =  2.0e1 * QPOASES_EPS;
	THIS->numRefinementSteps            =  2;
	#else
	THIS->epsRegularisation             =  1.0e3 * QPOASES_EPS;
	THIS->numRefinementSteps            =  1;
	#endif
	THIS->epsIterRef                    =  1.0e2 * QPOASES_EPS;
	#ifdef __USE_SINGLE_PRECISION__
	THIS->epsLITests                    =  5.0e1 * QPOASES_EPS;
	THIS->epsNZCTests                   =  1.0e2 * QPOASES_EPS;
	#else
	THIS->epsLITests                    =  1.0e5 * QPOASES_EPS;
	THIS->epsNZCTests                   =  3.0e3 * QPOASES_EPS;
	#endif

	THIS->enableDropInfeasibles         =  BT_FALSE;
    THIS->dropBoundPriority             =  1;
    THIS->dropEqConPriority             =  1;
    THIS->dropIneqConPriority           =  1;
    
	return SUCCESSFUL_RETURN;
}


/*
 *	s e t T o R e l i a b l e
 */
returnValue Options_setToReliable(	Options* THIS
									)
{
	Options_setToDefault( THIS );

	THIS->enableFullLITests             =  BT_TRUE;
	THIS->enableCholeskyRefactorisation =  1;

	#ifdef __USE_SINGLE_PRECISION__
	THIS->numRefinementSteps            =  3;
	#else
 	THIS->numRefinementSteps            =  2;
	#endif /*__USE_SINGLE_PRECISION__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t T o M P C
 */
returnValue Options_setToMPC(	Options* THIS
								)
{
	Options_setToDefault( THIS );

	THIS->enableRamping                 =  BT_FALSE;
	THIS->enableFarBounds               =  BT_TRUE;
	THIS->enableFlippingBounds          =  BT_FALSE;
	THIS->enableRegularisation          =  BT_TRUE;
	THIS->enableNZCTests                =  BT_FALSE;
	THIS->enableDriftCorrection         =  0;
	THIS->enableEqualities              =  BT_TRUE;

	#ifdef __USE_SINGLE_PRECISION__
	THIS->terminationTolerance          =  1.0e3 * QPOASES_EPS;
	#else
	THIS->terminationTolerance          =  1.0e9 * QPOASES_EPS;
	#endif

	THIS->initialStatusBounds           =  ST_INACTIVE;
	THIS->numRegularisationSteps        =  1;
	#ifdef __USE_SINGLE_PRECISION__
	THIS->numRefinementSteps            =  2;
	#else
	THIS->numRefinementSteps            =  0;
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t T o F a s t
 */
returnValue Options_setToFast(	Options* THIS
								)
{
	return Options_setToMPC( THIS );
}



/*
 *	e n s u r e C o n s i s t e n c y
 */
returnValue Options_ensureConsistency(	Options* THIS
										)
{
	BooleanType needToAdjust = BT_FALSE;

	/* flipping bounds require far bounds */
    /* (ckirches) Removed this as per filter's trust region
	if( enableFlippingBounds == BT_TRUE )
		enableFarBounds = BT_TRUE;
    */
	
	if( THIS->enableDriftCorrection < 0 )
	{
		THIS->enableDriftCorrection = 0;
		needToAdjust = BT_TRUE;
	}
	
	if( THIS->enableCholeskyRefactorisation < 0 )
	{
		THIS->enableCholeskyRefactorisation = 0;
		needToAdjust = BT_TRUE;
	}


	if ( THIS->terminationTolerance <= 0.0 )
	{
		THIS->terminationTolerance = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->epsIterRef <= 0.0 )
	{
		THIS->epsIterRef = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->epsRegularisation <= 0.0 )
	{
		THIS->epsRegularisation = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->boundTolerance <= 0.0 )
	{
		THIS->boundTolerance = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->boundRelaxation <= 0.0 )
	{
		THIS->boundRelaxation = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}
	
	if ( THIS->maxPrimalJump <= 0.0 )
	{
		THIS->maxPrimalJump = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->maxDualJump <= 0.0 )
	{
		THIS->maxDualJump = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}


	if ( THIS->initialRamping < 0.0 )
	{
		THIS->initialRamping = 0.0;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->finalRamping < 0.0 )
	{
		THIS->finalRamping = 0.0;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->initialFarBounds <= THIS->boundRelaxation )
	{
		THIS->initialFarBounds = THIS->boundRelaxation+QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}
	
	if ( THIS->growFarBounds < 1.1 )
	{
		THIS->growFarBounds = 1.1;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->epsFlipping <= 0.0 )
	{
		THIS->epsFlipping = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->numRegularisationSteps < 0 )
	{
		THIS->numRegularisationSteps = 0;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->epsRegularisation < 0.0 )
	{
		THIS->epsRegularisation = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->numRefinementSteps < 0 )
	{
		THIS->numRefinementSteps = 0;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->epsIterRef < 0.0 )
	{
		THIS->epsIterRef = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->epsLITests < 0.0 )
	{
		THIS->epsLITests = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( THIS->epsNZCTests < 0.0 )
	{
		THIS->epsNZCTests = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( needToAdjust == BT_TRUE)
		return THROWWARNING( RET_OPTIONS_ADJUSTED );

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue Options_print(	Options* THIS
							)
{
	#ifndef __XPCTARGET__
	#ifndef __DSPACE__
	static char myPrintfString[QPOASES_MAX_STRING_LENGTH];
	static char info[QPOASES_MAX_STRING_LENGTH];

	qpOASES_myPrintf( "\n###################   qpOASES  --  QP OPTIONS   ##################\n" );
	qpOASES_myPrintf( "\n" );

	qpOASES_convertPrintLevelToString( THIS->printLevel,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"printLevel                     =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_myPrintf( "\n" );

	qpOASES_convertBooleanTypeToString( THIS->enableRamping,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableRamping                  =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( THIS->enableFarBounds,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableFarBounds                =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( THIS->enableFlippingBounds,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableFlippingBounds           =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( THIS->enableRegularisation,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableRegularisation           =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( THIS->enableFullLITests,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableFullLITests              =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );
	
	qpOASES_convertBooleanTypeToString( THIS->enableNZCTests,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableNZCTests                 =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableDriftCorrection          =  %d\n",THIS->enableDriftCorrection );
	qpOASES_myPrintf( myPrintfString );
	
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableCholeskyRefactorisation  =  %d\n",THIS->enableCholeskyRefactorisation );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( THIS->enableEqualities,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableEqualities               =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_myPrintf( "\n" );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"terminationTolerance           =  %e\n",THIS->terminationTolerance );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"boundTolerance                 =  %e\n",THIS->boundTolerance );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"boundRelaxation                =  %e\n",THIS->boundRelaxation );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsNum                         =  %e\n",THIS->epsNum );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsDen                         =  %e\n",THIS->epsDen );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"maxPrimalJump                  =  %e\n",THIS->maxPrimalJump );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"maxDualJump                    =  %e\n",THIS->maxDualJump );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_myPrintf( "\n" );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"initialRamping                 =  %e\n",THIS->initialRamping );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"finalRamping                   =  %e\n",THIS->finalRamping );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"initialFarBounds               =  %e\n",THIS->initialFarBounds );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"growFarBounds                  =  %e\n",THIS->growFarBounds );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertSubjectToStatusToString( THIS->initialStatusBounds,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"initialStatusBounds            =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsFlipping                    =  %e\n",THIS->epsFlipping );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"numRegularisationSteps         =  %d\n",THIS->numRegularisationSteps );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsRegularisation              =  %e\n",THIS->epsRegularisation );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"numRefinementSteps             =  %d\n",THIS->numRefinementSteps );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsIterRef                     =  %e\n",THIS->epsIterRef );
	qpOASES_myPrintf( myPrintfString );
	
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsLITests                     =  %e\n",THIS->epsLITests );
	qpOASES_myPrintf( myPrintfString );
	
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsNZCTests                    =  %e\n",THIS->epsNZCTests );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_myPrintf( "\n\n" );
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/



END_NAMESPACE_QPOASES


/*
 *	end of file
 */
