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
 *	\file interfaces/matlab/qpOASES.c
 *	\author Hans Joachim Ferreau, Aude Perrin
 *	\version 3.0embedded
 *	\date 2007-2014
 *
 *	Interface for Matlab(R) that enables to call qpOASES as a MEX function.
 *
 */


#include <qpOASES.h>

USING_NAMESPACE_QPOASES

#include "qpOASES_matlab_utils.c"


#define NPMAX 200

/*
 *	q p O A S E S m e x _ c o n s t r a i n t s
 */
int qpOASESmex_constraints(		int nV, int nC, int nP,
								DenseMatrix* H, real_t* g, DenseMatrix* A,
								real_t* lb, real_t* ub, real_t* lbA, real_t* ubA,
								int nWSRin, real_t* x0, Options* options,
								int nOutputs, mxArray* plhs[],
								double* guessedBounds, double* guessedConstraints
								)
{
	int i,k;
	returnValue returnvalue;

	QProblem QP;
	Bounds bounds;
	Constraints constraints;
    
	real_t *g_current, *lb_current, *ub_current, *lbA_current, *ubA_current;

    int nWSRout;
        
	/* 1) Setup initial QP. */
	QProblemCON( &QP,nV,nC,HST_SEMIDEF ); //ensure Hessian regularisation if necessary
	QProblem_setOptions( &QP,*options );

	/* 2) Solve initial QP. */
	BoundsCON( &bounds,nV );
	ConstraintsCON( &constraints,nC );

	if (guessedBounds != 0) {
		mexPrintf( "GUESSED BOUNDS!!\n" );
		for (i = 0; i < nV; i++) {
			if ( qpOASES_isEqual(guessedBounds[i],-1.0,QPOASES_TOL) == BT_TRUE ) {
				Bounds_setupBound( &bounds,i, ST_LOWER);
			} else if ( qpOASES_isEqual(guessedBounds[i],1.0,QPOASES_TOL) == BT_TRUE ) {
				Bounds_setupBound( &bounds,i, ST_UPPER);
			} else if ( qpOASES_isEqual(guessedBounds[i],0.0,QPOASES_TOL) == BT_TRUE ) {
				Bounds_setupBound( &bounds,i, ST_INACTIVE);
			} else {
				char msg[200];
				snprintf(msg, 199,
						"ERROR (qpOASES): Only {-1, 0, 1} allowed for status of bounds!");
				myMexErrMsgTxt(msg);
				return -1;
			}
		}
	}

	if (guessedConstraints != 0) {
		mexPrintf( "GUESSED CONSTRAINTS!!\n" );
		for (i = 0; i < nC; i++) {
			if ( qpOASES_isEqual(guessedConstraints[i],-1.0,QPOASES_TOL) == BT_TRUE ) {
				Constraints_setupConstraint( &constraints,i, ST_LOWER);
			} else if ( qpOASES_isEqual(guessedConstraints[i],1.0,QPOASES_TOL) == BT_TRUE ) {
				Constraints_setupConstraint( &constraints,i, ST_UPPER);
			} else if ( qpOASES_isEqual(guessedConstraints[i],0.0,QPOASES_TOL) == BT_TRUE ) {
				Constraints_setupConstraint( &constraints,i, ST_INACTIVE);
			} else {
				char msg[200];
				snprintf(msg, 199,
						"ERROR (qpOASES): Only {-1, 0, 1} allowed for status of constraints!");
				myMexErrMsgTxt(msg);
				return -1;
			}
		}
	}
	
    nWSRout = nWSRin;
	if (x0 == 0 && guessedBounds == 0 && guessedConstraints == 0)
	{
		returnvalue = QProblem_initM( &QP,H, g, A, lb, ub, lbA, ubA, &nWSRout, 0);
	}
	else
	{
		returnvalue = QProblem_initMW( &QP,H, g, A, lb, ub, lbA, ubA, &nWSRout, 0, x0, 0,
				guessedBounds != 0 ? &bounds : 0,
				guessedConstraints != 0 ? &constraints : 0);
	}

	/* 3) Solve remaining QPs and assign lhs arguments. */
	/*    Set up pointers to the current QP vectors */
	g_current   = g;
	lb_current  = lb;
	ub_current  = ub;
	lbA_current = lbA;
	ubA_current = ubA;

	/* Loop through QP sequence. */
	for ( k=0; k<nP; ++k )
	{
		if ( k != 0 )
		{
			/* update pointers to the current QP vectors */
			g_current = &(g[k*nV]);
			if ( lb != 0 )
				lb_current = &(lb[k*nV]);
			if ( ub != 0 )
				ub_current = &(ub[k*nV]);
			if ( lbA != 0 )
				lbA_current = &(lbA[k*nC]);
			if ( ubA != 0 )
				ubA_current = &(ubA[k*nC]);

			nWSRout = nWSRin;
			returnvalue = QProblem_hotstart( &QP,g_current,lb_current,ub_current,lbA_current,ubA_current, &nWSRout,0 );
		}

		/* write results into output vectors */
		obtainOutputs(	k,&QP,returnvalue,nWSRout,
						nOutputs,plhs,nV,nC,-1 );
	}

	//QP.writeQpDataIntoMatFile( "qpDataMat0.mat" );

	return 0;
}


/*
 *	q p O A S E S m e x _ b o u n d s
 */
int qpOASESmex_bounds(	int nV, int nP,
						DenseMatrix* H, real_t* g,
						real_t* lb, real_t* ub,
						int nWSRin, real_t* x0, Options* options,
						int nOutputs, mxArray* plhs[],
						double* guessedBounds
						)
{
	int i,k;
	returnValue returnvalue;

	QProblemB QP;
	Bounds bounds;

	real_t *g_current, *lb_current, *ub_current;
    
    int nWSRout;

	/* 1) Setup initial QP. */
	QProblemBCON( &QP,nV,HST_SEMIDEF ); //ensure Hessian regularisation if necessary
	QProblemB_setOptions( &QP,*options );

	/* 2) Solve initial QP. */
	BoundsCON( &bounds,nV );

	/* 2) Solve initial QP. */
	if (guessedBounds != 0) {
		for (i = 0; i < nV; i++) {
			if ( qpOASES_isEqual(guessedBounds[i],-1.0,QPOASES_TOL) == BT_TRUE ) {
				Bounds_setupBound( &bounds,i, ST_LOWER);
			} else if ( qpOASES_isEqual(guessedBounds[i],1.0,QPOASES_TOL) == BT_TRUE ) {
				Bounds_setupBound( &bounds,i, ST_UPPER);
			} else if ( qpOASES_isEqual(guessedBounds[i],0.0,QPOASES_TOL) == BT_TRUE ) {
				Bounds_setupBound( &bounds,i, ST_INACTIVE);
			} else {
				char msg[200];
				snprintf(msg, 199,
						"ERROR (qpOASES): Only {-1, 0, 1} allowed for status of bounds!");
				myMexErrMsgTxt(msg);
				return -1;
			}
		}
	}

    nWSRout = nWSRin;
	if (x0 == 0 && guessedBounds == 0)
		returnvalue = QProblemB_initM( &QP,H, g, lb, ub, &nWSRout, 0);
	else
		returnvalue = QProblemB_initMW( &QP,H, g, lb, ub, &nWSRout, 0, x0, 0,
				guessedBounds != 0 ? &bounds : 0);

	/* 3) Solve remaining QPs and assign lhs arguments. */
	/*    Set up pointers to the current QP vectors */
	g_current  = g;
	lb_current = lb;
	ub_current = ub;

	/* Loop through QP sequence. */
	for ( k=0; k<nP; ++k )
	{
		if ( k != 0 )
		{
			/* update pointers to the current QP vectors */
			g_current = &(g[k*nV]);
			if ( lb != 0 )
				lb_current = &(lb[k*nV]);
			if ( ub != 0 )
				ub_current = &(ub[k*nV]);

            nWSRout = nWSRin;
			returnvalue = QProblemB_hotstart( &QP,g_current,lb_current,ub_current, &nWSRout,0 );
		}

		/* write results into output vectors */
		obtainOutputsB(	k,&QP,returnvalue,nWSRout,
						nOutputs,plhs,nV,-1 );
	}

	return 0;
}



/*
 *	m e x F u n c t i o n
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
	int i;
	//unsigned int i;

	/* inputs */
	DenseMatrix H;
	DenseMatrix A;
	real_t g[NVMAX*NPMAX];
	real_t lb[NVMAX*NPMAX];
	real_t ub[NVMAX*NPMAX];
	real_t lbA[NCMAX*NPMAX];
	real_t ubA[NCMAX*NPMAX];
	real_t x0[NVMAX];

	int H_idx, g_idx, A_idx, lb_idx, ub_idx, lbA_idx, ubA_idx, x0_idx=-1, options_idx=-1;

	real_t* A_for=0;
	real_t A_mem[NCMAX*NVMAX];

	double guessedBoundsAndConstraints[NVMAX+NCMAX];
	double guessedBounds[NVMAX];
	double guessedConstraints[NCMAX];
	int guessedBoundsAndConstraintsIndex = -1;

	/* dimensions */
	unsigned int nV=0, nC=0, nP=0;
	int nWSRin;
	int numberOfColumns;

    /* Setup default options */
	Options options;
	Options_setToDefault( &options );
	options.printLevel = PL_LOW;
	#ifdef __DEBUG__
	options.printLevel = PL_HIGH;
	#endif
	#ifdef __SUPPRESSANYOUTPUT__
	options.printLevel = PL_NONE;
	#endif


	/* I) CONSISTENCY CHECKS: */
	/* 1) Ensure that qpOASES is called with a feasible number of input arguments. */
	if ( ( nrhs < 4 ) || ( nrhs > 10 ) )
	{
		myMexErrMsgTxt( "ERROR (qpOASES): Invalid number of input arguments!\nType 'help qpOASES' for further information." );
		return;
	}

	/* 2) Check for proper number of output arguments. */
	if ( nlhs > 6 )
	{
		myMexErrMsgTxt( "ERROR (qpOASES): At most six output arguments are allowed: \n    [x,fval,exitflag,iter,lambda,activeSet]!" );
		return;
	}
	if ( nlhs < 1 )
	{
		myMexErrMsgTxt( "ERROR (qpOASES): At least one output argument is required: [x,...]!" );
		return;
	}


	/* II) PREPARE RESPECTIVE QPOASES FUNCTION CALL: */
	/*     Choose between QProblem and QProblemB object and assign the corresponding
	 *     indices of the input pointer array in to order to access QP data correctly. */
	H_idx = 0;
	g_idx = 1;
	nV = mxGetM( prhs[ H_idx ] ); /* row number of Hessian matrix */
	nP = mxGetN( prhs[ g_idx ] ); /* number of columns of the gradient matrix (vectors series have to be stored columnwise!) */

	/* 0) Check whether options are specified .*/
	for (i = 4; i < nrhs; i++) {
		if ((!mxIsEmpty(prhs[i])) && (mxIsStruct(prhs[i]))) {
			options_idx = i;
			break;
		}
	}

	// Is the third argument constraint Matrix A?
	numberOfColumns = mxGetN(prhs[2]);

	/* 1) Simply bounded QP. */
	if ( ( nrhs <= 6 ) || 
		 ( ( nrhs == 7 ) && ( options_idx == 6 ) ) || 
		 ( ( numberOfColumns == 1 ) && ( nV != 1 ) ) )
	{
		lb_idx   = 2;
		ub_idx   = 3;

		if ( nrhs >= 6 ){ /* x0 specified */
			x0_idx = 5;
			if (nrhs >= 7) {
				guessedBoundsAndConstraintsIndex = 6;
			} else {
				guessedBoundsAndConstraintsIndex = -1;
			}
		}
		else {
			x0_idx = -1;
			guessedBoundsAndConstraintsIndex = -1;
		}
	}
	else
	{
		A_idx = 2;

		/* If constraint matrix is empty, use a QProblemB object! */
		if ( mxIsEmpty( prhs[ A_idx ] ) )
		{
			lb_idx   = 3;
			ub_idx   = 4;

			if ( nrhs >= 9 ) /* x0 specified */
				x0_idx = 8;
			else
				x0_idx = -1;

			/* guesses for bounds and constraints specified? */
			if (nrhs >= 10) {
				guessedBoundsAndConstraintsIndex = 9;
			} else {
				guessedBoundsAndConstraintsIndex = -1;
			}
		}
		else
		{
			lb_idx   = 3;
			ub_idx   = 4;
			lbA_idx  = 5;
			ubA_idx  = 6;

			if ( nrhs >= 9 ) /* x0 specified */
				x0_idx = 8;
			else
				x0_idx = -1;

			/* guesses for bounds and constraints specified? */
			if (nrhs >= 10) {
				guessedBoundsAndConstraintsIndex = 9;
			} else {
				guessedBoundsAndConstraintsIndex = -1;
			}

			nC = mxGetM( prhs[ A_idx ] ); /* row number of constraint matrix */
		}
	}

	/*mexPrintf( "H_idx = %d\ng_idx = %d\nA_idx = %d\nlb_idx = %d\nub_idx = %d\nlbA_idx = %d\nubA_idx = %d\nx0_idx = %d\nopt_idx = %d\nguessedBoundsAndConstraintsIndex = %d\n",
				H_idx,g_idx,A_idx,lb_idx,ub_idx,lbA_idx,ubA_idx,x0_idx,options_idx,guessedBoundsAndConstraintsIndex );*/

	/* III) ACTUALLY PERFORM QPOASES FUNCTION CALL: */
	nWSRin = 5*(nV+nC);
	if ( options_idx > 0 )
		setupOptions( &options,prhs[options_idx],&nWSRin );

	/* ensure that data is given in real_t precision */
	if ( ( mxIsDouble( prhs[ H_idx ] ) == 0 ) ||
		 ( mxIsDouble( prhs[ g_idx ] ) == 0 ) )
	{
		myMexErrMsgTxt( "ERROR (qpOASES): All data has to be provided in double precision!" );
		return;
	}

	/* Check inputs dimensions and assign pointers to inputs. */
	if ( mxGetN( prhs[ H_idx ] ) != nV )
	{
		char msg[200]; 
		snprintf(msg, 199, "ERROR (qpOASES): Hessian matrix input dimension mismatch (%ld != %d)!", 
				(long int)mxGetN(prhs[H_idx]), nV);
		myMexErrMsgTxt(msg);
		return;
	}

	/* check for sparsity */
	if ( mxIsSparse( prhs[ H_idx ] ) != 0 )
	{
		myMexErrMsgTxt( "ERROR (qpOASES): Cannot handle sparse matrices in embedded variant!" );
		return;
	}
	else
	{
		DenseMatrixCON( &H,nV, nV, nV, (real_t*) mxGetPr(prhs[H_idx]));
	}

	/*DenseMatrix_print( &H );*/

	if ( smartDimensionCheck( g,nV,nP, BT_FALSE,prhs,&g_idx ) != SUCCESSFUL_RETURN )
		return;
	/*if ( g_idx >= 0 ) qpOASES_printV( g,nV );*/
	
	if ( smartDimensionCheck( lb,nV,nP, BT_TRUE,prhs,&lb_idx ) != SUCCESSFUL_RETURN )
		return;
	/*if ( lb_idx >= 0 ) qpOASES_printV( lb,nV );*/

	if ( smartDimensionCheck( ub,nV,nP, BT_TRUE,prhs,&ub_idx ) != SUCCESSFUL_RETURN )
		return;
	/*if ( ub_idx >= 0 ) qpOASES_printV( ub,nV );*/
	
	if ( smartDimensionCheck( x0,nV,1, BT_TRUE,prhs,&x0_idx ) != SUCCESSFUL_RETURN )
		return;

	if (guessedBoundsAndConstraintsIndex >= 0) {
		if (smartDimensionCheck(guessedBoundsAndConstraints, nV + nC, 1,
				BT_TRUE, prhs, &guessedBoundsAndConstraintsIndex)
				!= SUCCESSFUL_RETURN) {
			return;
		}
	}

	if (guessedBoundsAndConstraintsIndex >= 0) {
		mexPrintf( "guessedBoundsAndConstraintsIndex >=0\n" );
		for (i = 0; i < nV; i++) {
			guessedBounds[i] = guessedBoundsAndConstraints[i];
		}

		if ( nC > 0 )
		{
			for (i = 0; i < nC; i++) {
				guessedConstraints[i] = guessedBoundsAndConstraints[i + nV];
			}
		}
	}

	if ( nC > 0 )
	{
		/* ensure that data is given in real_t precision */
		if ( mxIsDouble( prhs[ A_idx ] ) == 0 )
		{
			myMexErrMsgTxt( "ERROR (qpOASES): All data has to be provided in real_t precision!" );
			return;
		}

		/* Check inputs dimensions and assign pointers to inputs. */
		if ( mxGetN( prhs[ A_idx ] ) != nV )
		{
			char msg[200]; 
			snprintf(msg, 199, "ERROR (qpOASES): Constraint matrix input dimension mismatch (%ld != %d)!", 
					(long int)mxGetN(prhs[A_idx]), nV);
			//fprintf(stderr, "%s\n", msg);
			myMexErrMsgTxt(msg);
			return;
		}

		A_for = (real_t*) mxGetPr( prhs[ A_idx ] );
		/*mexPrintf( "nV = %d\nnC = %d\nnP= %d\n", nV,nC,nP );*/
		/*mexPrintf( "nV = %d\nnC = %d\nnP= %d\nlbA_idx = %d\n", nV,nC,nP,lbA_idx );*/

		if ( smartDimensionCheck( lbA,nC,nP, BT_TRUE,prhs,&lbA_idx ) != SUCCESSFUL_RETURN )
			return;
		/* if ( lbA_idx >= 0 ) qpOASES_printV( lbA,nC ); */

		/*mexPrintf( "nV = %d\nnC = %d\nnP= %d\nubA_idx = %d\n", nV,nC,nP,ubA_idx );*/

		if ( smartDimensionCheck( ubA,nC,nP, BT_TRUE,prhs,&ubA_idx ) != SUCCESSFUL_RETURN )
			return;
		/* if ( ubA_idx >= 0 ) qpOASES_printV( ubA,nC ); */

		/* Check for sparsity. */
		if ( mxIsSparse( prhs[ A_idx ] ) != 0 )
		{
			myMexErrMsgTxt( "ERROR (qpOASES): All data has to be provided in real_t precision!" );
			return;
		}
		else
		{
			/* Convert constraint matrix A from FORTRAN to C style
			 * (not necessary for H as it should be symmetric!). */
			convertFortranToC( A_for,nV,nC, A_mem );
			DenseMatrixCON( &A,nC, nV, nV, A_mem);
		}
	}


	/*if ( A_idx >= 0 ) DenseMatrix_print( &A );*/

	allocateOutputs( nlhs,plhs,nV,nC,nP,-1 );

	if ( nC == 0 )
	{
		/* call qpOASES */
		qpOASESmex_bounds(	nV,nP,
							&H,
							g_idx >= 0  ? g  : 0,
							lb_idx >= 0 ? lb : 0,
							ub_idx >= 0 ? ub : 0,
							nWSRin,
							x0_idx >= 0 ? x0 : 0,
							&options,
							nlhs,plhs,
							guessedBoundsAndConstraintsIndex >= 0 ? guessedBounds : 0
							);

		return;
		/* 2) Call usual version including constraints (using QProblem class) */
	}
	else
	{
		/* Call qpOASES. */
		qpOASESmex_constraints(	nV,nC,nP,
								&H,
								g_idx >= 0   ? g   : 0,
								&A,
								lb_idx >= 0  ? lb  : 0,
								ub_idx >= 0  ? ub  : 0,
								lbA_idx >= 0 ? lbA : 0,
								ubA_idx >= 0 ? ubA : 0,
								nWSRin,
								x0_idx >= 0  ? x0  : 0,
								&options,
								nlhs,plhs,
								guessedBoundsAndConstraintsIndex >= 0 ? guessedBounds      : 0,
								guessedBoundsAndConstraintsIndex >= 0 ? guessedConstraints : 0
								);
		return;
	}
}

/*
 *	end of file
 */
