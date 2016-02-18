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
 *	\file src/Matrices.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0embedded
 *	\date 2007-2014
 *
 *	Implementation of the matrix classes.
 */



#include <qpOASES/Matrices.h>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


void DenseMatrixCON(	DenseMatrix* THIS,
						int m,
						int n,
						int lD,
						real_t *v
						)
{
	DenseMatrix_init( THIS,m,n,lD,v );
}


void DenseMatrixCPY(	DenseMatrix* FROM,
						DenseMatrix* TO
						)
{
	TO->nRows  = FROM->nRows;
	TO->nCols  = FROM->nCols;
	TO->leaDim = FROM->leaDim;
	TO->val    = FROM->val; 
}



void DenseMatrix_free( DenseMatrix* THIS )
{
	THIS->val = 0;
}


returnValue DenseMatrix_init(	DenseMatrix* THIS,
								int m,
								int n,
								int lD,
								real_t *v
								)
{
	THIS->nRows  = m; 
	THIS->nCols  = n;
	THIS->leaDim = lD;
	THIS->val    = v;
	return SUCCESSFUL_RETURN;
}



real_t DenseMatrix_diag(	DenseMatrix* THIS,
							int i
							)
{
	return THIS->val[i*((THIS->leaDim)+1)];
}


BooleanType DenseMatrix_isDiag( DenseMatrix* THIS )
{
	int i, j;

	if (THIS->nRows != THIS->nCols)
		return BT_FALSE;

	for ( i=0; i<THIS->nRows; ++i )
		for ( j=0; j<i; ++j )
			if ( ( fabs( THIS->val[i*(THIS->leaDim)+j] ) > QPOASES_EPS ) || ( fabs( THIS->val[j*(THIS->leaDim)+i] ) > QPOASES_EPS ) )
				return BT_FALSE;

	return BT_TRUE;
}


real_t DenseMatrix_getNorm( DenseMatrix* THIS, int type )
{
    return REFER_NAMESPACE_QPOASES qpOASES_getNorm( THIS->val,(THIS->nCols)*(THIS->nRows),type ); 
}


real_t DenseMatrix_getRowNorm( DenseMatrix* THIS, int rNum, int type )
{
    return REFER_NAMESPACE_QPOASES qpOASES_getNorm( &(THIS->val[rNum*(THIS->leaDim)]),THIS->nCols,type ); 
}


returnValue DenseMatrix_getRow( DenseMatrix* THIS,int rNum, const Indexlist* const icols, real_t alpha, real_t *row)
{
	int i;
    if (icols != 0)
    {
	    if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		    for (i = 0; i < icols->length; i++)
			    row[i] = THIS->val[rNum*(THIS->leaDim)+icols->number[i]];
	    else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		    for (i = 0; i < icols->length; i++)
			    row[i] = -THIS->val[rNum*(THIS->leaDim)+icols->number[i]];
	    else
		    for (i = 0; i < icols->length; i++)
			    row[i] = alpha*THIS->val[rNum*(THIS->leaDim)+icols->number[i]];
    }
    else
    {
	    if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		    for (i = 0; i < THIS->nCols; i++)
			    row[i] = THIS->val[rNum*(THIS->leaDim)+i];
	    else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		    for (i = 0; i < THIS->nCols; i++)
			    row[i] = -THIS->val[rNum*(THIS->leaDim)+i];
	    else
		    for (i = 0; i < THIS->nCols; i++)
			    row[i] = alpha*THIS->val[rNum*(THIS->leaDim)+i];
    }
	return SUCCESSFUL_RETURN;
}


returnValue DenseMatrix_getCol( DenseMatrix* THIS,int cNum, const Indexlist* const irows, real_t alpha, real_t *col)
{
	int i;

	if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (i = 0; i < irows->length; i++)
			col[i] = THIS->val[irows->number[i]*(THIS->leaDim)+cNum];
	else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (i = 0; i < irows->length; i++)
			col[i] = -THIS->val[irows->number[i]*(THIS->leaDim)+cNum];
	else
		for (i = 0; i < irows->length; i++)
			col[i] = alpha*THIS->val[irows->number[i]*(THIS->leaDim)+cNum];

	return SUCCESSFUL_RETURN;
}



returnValue DenseMatrix_times(	DenseMatrix* THIS,
								int xN, real_t alpha, const real_t *x, int xLD, real_t beta, real_t *y, int yLD
								)
{
	unsigned long _xN     = (unsigned long)xN;
	unsigned long _nRows  = (unsigned long)(THIS->nRows);
	unsigned long _nCols  = (unsigned long)(THIS->nCols);
	unsigned long _leaDim = (unsigned long)qpOASES_getMax(1,THIS->nCols);
	unsigned long _xLD    = (unsigned long)qpOASES_getMax(1,xLD);
	unsigned long _yLD    = (unsigned long)qpOASES_getMax(1,yLD);

	/* Call BLAS. Mind row major format! */
	GEMM("TRANS", "NOTRANS", &_nRows, &_xN, &_nCols, &alpha, THIS->val, &_leaDim, x, &_xLD, &beta, y, &_yLD);
	return SUCCESSFUL_RETURN;
}


returnValue DenseMatrix_transTimes(	DenseMatrix* THIS,
									int xN, real_t alpha, const real_t *x, int xLD, real_t beta, real_t *y, int yLD
									)
{
	unsigned long _xN     = (unsigned long)xN;
	unsigned long _nRows  = (unsigned long)(THIS->nRows);
	unsigned long _nCols  = (unsigned long)(THIS->nCols);
	unsigned long _leaDim = (unsigned long)qpOASES_getMax(1,THIS->nCols);
	unsigned long _xLD    = (unsigned long)qpOASES_getMax(1,xLD);
	unsigned long _yLD    = (unsigned long)qpOASES_getMax(1,yLD);

	/* Call BLAS. Mind row major format! */
	GEMM("NOTRANS", "NOTRANS", &_nCols, &_xN, &_nRows, &alpha, THIS->val, &_leaDim, x, &_xLD, &beta, y, &_yLD);
	return SUCCESSFUL_RETURN;
}


returnValue DenseMatrix_subTimes(	DenseMatrix* THIS,
									const Indexlist* const irows, const Indexlist* const icols,
									int xN, real_t alpha, const real_t *x, int xLD, real_t beta, real_t *y, int yLD,
									BooleanType yCompr
									)
{
	int i, j, k, row, col, iy, irA;

	if (yCompr == BT_TRUE)
	{
		if ( REFER_NAMESPACE_QPOASES qpOASES_isZero(beta,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = 0.0;
		else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(beta,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = -y[j+k*yLD];
		else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(beta,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_FALSE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] *= beta;

		if (icols == 0)
			if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * (THIS->leaDim);
						for (i = 0; i < THIS->nCols; i++)
							y[iy] += THIS->val[irA+i] * x[k*xLD+i];
					}
			else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * (THIS->leaDim);
						for (i = 0; i < THIS->nCols; i++)
							y[iy] -= THIS->val[irA+i] * x[k*xLD+i];
					}
			else
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * (THIS->leaDim);
						for (i = 0; i < THIS->nCols; i++)
							y[iy] += alpha * THIS->val[irA+i] * x[k*xLD+i];
					}
		else /* icols != 0 */
			if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * (THIS->leaDim);
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] += THIS->val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
			else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * (THIS->leaDim);
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] -= THIS->val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
			else
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * (THIS->leaDim);
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] += alpha * THIS->val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
	}
	else /* y not compressed */
	{
		if ( REFER_NAMESPACE_QPOASES qpOASES_isZero(beta,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = 0.0;
		else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(beta,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = -y[j+k*yLD];
		else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(beta,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_FALSE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] *= beta;

		if (icols == 0)
			if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * (THIS->leaDim);
						for (i = 0; i < THIS->nCols; i++)
							y[iy] += THIS->val[irA+i] * x[k*xLD+i];
					}
			else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * (THIS->leaDim);
						for (i = 0; i < THIS->nCols; i++)
							y[iy] -= THIS->val[irA+i] * x[k*xLD+i];
					}
			else
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * (THIS->leaDim);
						for (i = 0; i < THIS->nCols; i++)
							y[iy] += alpha * THIS->val[irA+i] * x[k*xLD+i];
					}
		else /* icols != 0 */
			if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * (THIS->leaDim);
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] += THIS->val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
			else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * (THIS->leaDim);
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] -= THIS->val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
			else
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * (THIS->leaDim);
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] += alpha * THIS->val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
	}

	return SUCCESSFUL_RETURN;
}

returnValue DenseMatrix_subTransTimes(	DenseMatrix* THIS,
										const Indexlist* const irows, const Indexlist* const icols,
										int xN, real_t alpha, const real_t *x, int xLD, real_t beta, real_t *y, int yLD
										)
{
	int i, j, k, row, col;

	if ( REFER_NAMESPACE_QPOASES qpOASES_isZero(beta,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] = 0.0;
	else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(beta,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] = -y[j+k*yLD];
	else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(beta,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_FALSE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] *= beta;

	if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < irows->length; j++)
			{
				row = irows->iSort[j];
				for (i = 0; i < icols->length; i++)
				{
					col = icols->iSort[i];
					y[col+k*yLD] += THIS->val[irows->number[row]*(THIS->leaDim)+icols->number[col]] * x[row+k*xLD];
				}
			}
	else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(alpha,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < irows->length; j++)
			{
				row = irows->iSort[j];
				for (i = 0; i < icols->length; i++)
				{
					col = icols->iSort[i];
					y[col+k*yLD] -= THIS->val[irows->number[row]*(THIS->leaDim)+icols->number[col]] * x[row+k*xLD];
				}
			}
	else
		for (k = 0; k < xN; k++)
			for (j = 0; j < irows->length; j++)
			{
				row = irows->iSort[j];
				for (i = 0; i < icols->length; i++)
				{
					col = icols->iSort[i];
					y[col+k*yLD] += alpha * THIS->val[irows->number[row]*(THIS->leaDim)+icols->number[col]] * x[row+k*xLD];
				}
			}

	return SUCCESSFUL_RETURN;
}


returnValue DenseMatrix_addToDiag( DenseMatrix* THIS,real_t alpha )
{
	int i;
	for (i = 0; i < THIS->nRows && i < THIS->nCols; i++)
		THIS->val[i*((THIS->leaDim)+1)] += alpha;

	return SUCCESSFUL_RETURN;
}


returnValue DenseMatrix_print( DenseMatrix* THIS )
{
	return qpOASES_printM( THIS->val,THIS->nRows,THIS->nCols );
}


returnValue DenseMatrix_bilinear(	DenseMatrix* THIS,
									const Indexlist* const icols, int xN, const real_t *x, int xLD, real_t *y, int yLD
									)
{
	int ii, jj, kk, col;
	int i,j,k,irA;

	myStatic real_t Ax[NVCMAX*NVMAX];
	real_t h;
	
	for (ii = 0; ii < xN; ii++)
		for (jj = 0; jj < xN; jj++)
			y[ii*yLD+jj] = 0.0;

	for (i=0;i<icols->length * xN;++i)
		Ax[i]=0.0;

	/* exploit symmetry of A ! */
	for (j = 0; j < icols->length; j++) {
		irA = icols->number[j] * (THIS->leaDim);
		for (i = 0; i < icols->length; i++)
		{
			h = THIS->val[irA+icols->number[i]];
			for (k = 0; k < xN; k++)
				Ax[j + k * icols->length] += h * x[k*xLD+icols->number[i]];
		}
	}

	for (ii = 0; ii < icols->length; ++ii) {
		col = icols->number[ii];
		for (jj = 0; jj < xN; ++jj) {
			for (kk = 0; kk < xN; ++kk) {
				y[kk + jj*yLD] += x[col + jj*xLD] * Ax[ii + kk*icols->length];
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


void dgemm_(	const char *TRANSA, const char *TRANSB,
				const unsigned long *M, const unsigned long *N, const unsigned long *K,
				const double *ALPHA, const double *A, const unsigned long *LDA, const double *B, const unsigned long *LDB,
				const double *BETA, double *C, const unsigned long *LDC
				)
{
	unsigned int i, j, k;

	if ( REFER_NAMESPACE_QPOASES qpOASES_isZero(*BETA,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] = 0.0;
	else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*BETA,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] = -C[j+(*LDC)*k];
	else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*BETA,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_FALSE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] *= *BETA;

	if (TRANSA[0] == 'N')
		if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*ALPHA,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += A[j+(*LDA)*i] * B[i+(*LDB)*k];
		else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*ALPHA,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] -= A[j+(*LDA)*i] * B[i+(*LDB)*k];
		else
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += *ALPHA * A[j+(*LDA)*i] * B[i+(*LDB)*k];
	else
		if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*ALPHA,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += A[i+(*LDA)*j] * B[i+(*LDB)*k];
		else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*ALPHA,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] -= A[i+(*LDA)*j] * B[i+(*LDB)*k];
		else
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += *ALPHA * A[i+(*LDA)*j] * B[i+(*LDB)*k];
}

void sgemm_(	const char *TRANSA, const char *TRANSB,
				const unsigned long *M, const unsigned long *N, const unsigned long *K,
				const float *ALPHA, const float *A, const unsigned long *LDA, const float *B, const unsigned long *LDB,
				const float *BETA, float *C, const unsigned long *LDC
				)
{
	unsigned int i, j, k;

	if ( REFER_NAMESPACE_QPOASES qpOASES_isZero(*BETA,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] = 0.0;
	else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*BETA,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] = -C[j+(*LDC)*k];
	else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*BETA,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_FALSE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] *= *BETA;

	if (TRANSA[0] == 'N')
		if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*ALPHA,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += A[j+(*LDA)*i] * B[i+(*LDB)*k];
		else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*ALPHA,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] -= A[j+(*LDA)*i] * B[i+(*LDB)*k];
		else
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += *ALPHA * A[j+(*LDA)*i] * B[i+(*LDB)*k];
	else
		if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*ALPHA,1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += A[i+(*LDA)*j] * B[i+(*LDB)*k];
		else if ( REFER_NAMESPACE_QPOASES qpOASES_isEqual(*ALPHA,-1.0,QPOASES_TOL) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] -= A[i+(*LDA)*j] * B[i+(*LDB)*k];
		else
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += *ALPHA * A[i+(*LDA)*j] * B[i+(*LDB)*k];
}



void dpotrf_(	const char *uplo, const unsigned long *_n, double *a,
				const unsigned long *_lda, long *info
				)
{
	double sum;
	long i, j, k;
	long n = (long)(*_n);
	long lda = (long)(*_lda);

	for( i=0; i<n; ++i )
	{
		/* j == i */
		sum = a[i + lda*i];

		for( k=(i-1); k>=0; --k )
			sum -= a[k+lda*i] * a[k+lda*i];

		if ( sum > 0.0 )
			a[i+lda*i] = qpOASES_getSqrt( sum );
		else
		{
			a[0] = sum; /* tunnel negative diagonal element to caller */
			if (info != 0)
				*info = (long)i+1;
			return;
		}

		for( j=(i+1); j<n; ++j )
		{
			sum = a[j*lda + i];

			for( k=(i-1); k>=0; --k )
				sum -= a[k+lda*i] * a[k+lda*j];

			a[i+lda*j] = sum / a[i+lda*i];
		}
	}
	if (info != 0)
		*info = 0;
}


void spotrf_(	const char *uplo, const unsigned long *_n, float *a,
				const unsigned long *_lda, long *info
				)
{
	float sum;
	long i, j, k;
	long n = (long)(*_n);
	long lda = (long)(*_lda);

	for( i=0; i<n; ++i )
	{
		/* j == i */
		sum = a[i + lda*i];

		for( k=(i-1); k>=0; --k )
			sum -= a[k+lda*i] * a[k+lda*i];

		if ( sum > 0.0 )
			a[i+lda*i] = (float)(REFER_NAMESPACE_QPOASES qpOASES_getSqrt( sum ));
		else
		{
			a[0] = sum; /* tunnel negative diagonal element to caller */
			if (info != 0)
				*info = (long)i+1;
			return;
		}

		for( j=(i+1); j<n; ++j )
		{
			sum = a[j*lda + i];

			for( k=(i-1); k>=0; --k )
				sum -= a[k+lda*i] * a[k+lda*j];

			a[i+lda*j] = sum / a[i+lda*i];
		}
	}
	if (info != 0)
		*info = 0;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
