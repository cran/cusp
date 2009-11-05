/*
 *  cusp.nc.c: calculate normalizing constant for cusp- (or Cobb-) density
 *  
 *
 *  Created by Raoul Grasman at the University of Amsterdam on 11/29/07.
 *  Copyright 2007 Raoul Grasman. 
 *
 *  This file is part of the cusp package for R.
 *
 *  The cusp package is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 */

#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>

#define fequal(a,b) (fabs((a)-(b))<1.0e-10)

/* Only this function is called via .External(.) :*/
SEXP cuspnc(SEXP args);

typedef struct par_struct
{
	double alpha; /* alpha (or 'normal' parameter-) values */
	double beta;  /* beta  (or 'splitting' parameter-) values */
} par_struct, *ParStruct;

static void cuspfun(double *x, int n, void *ex)
{
	double alpha, beta, y, y2;
	int i;

	ParStruct PS = (ParStruct) ex;
	
	alpha = PS->alpha;
	beta  = PS->beta;
	for(i=0; i<n; i++) {
		y = x[i];
		y2 = y*y;
		x[i] = exp(alpha*y+beta*0.5*y2-0.25*y2*y2);
		if(!R_FINITE(x[i])) {
			error("non-finite function value (alpha and/or beta are probably to extreme)");
		}
	}
	return;
}

SEXP cuspnc(SEXP args)
{
    par_struct par;
    SEXP ans, ansnames;
	SEXP alpha, beta;
    double bound, epsabs, epsrel, result, abserr, *work;
	double akeep, bkeep;
    int n, inf, neval, ier, limit, lenw, last, *iwork;

    args = CDR(args);
//	n = asInteger(CAR(args)); args = CDR(args);
	alpha = CAR(args); args = CDR(args);
	beta  = CAR(args); args = CDR(args);
    bound = asReal(CAR(args)); args = CDR(args);
    inf = asInteger(CAR(args)); args = CDR(args);
    epsabs = asReal(CAR(args)); args = CDR(args);
    epsrel = asReal(CAR(args)); args = CDR(args);
    limit = asInteger(CAR(args)); args = CDR(args);
    lenw = 4 * limit;
    iwork = (int *) R_alloc(limit, sizeof(int));
    work = (double *) R_alloc(lenw, sizeof(double));

	if(length(alpha) != length(beta)) {
		error("alpha and beta should have the same length");
	}
	n = length(alpha);


    PROTECT(ans = allocVector(VECSXP, 4));
    PROTECT(ansnames = allocVector(STRSXP, 4));
    SET_STRING_ELT(ansnames, 0, mkChar("value"));
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, n));
    SET_STRING_ELT(ansnames, 1, mkChar("abs.error"));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, n));
    SET_STRING_ELT(ansnames, 2, mkChar("subdivisions"));
    SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, n));
    SET_STRING_ELT(ansnames, 3, mkChar("ierr"));
    SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, n));

	akeep = bkeep = -999.0; // impossible value
	for(int i=0;i<n;i++){
		par.alpha = REAL(alpha)[i];
		par.beta  = REAL(beta)[i];
		
		// if alpha or beta is different from previous values, calculate normalizing constant
		if(!fequal(par.alpha, akeep) || !fequal(par.beta, bkeep)) {
			Rdqagi(cuspfun, (void*)&par,&bound,&inf,&epsabs,&epsrel,&result, // integrate
				&abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);
		}

		REAL(VECTOR_ELT(ans, 0))[i]    = result;
		REAL(VECTOR_ELT(ans, 1))[i]    = abserr;
		INTEGER(VECTOR_ELT(ans, 2))[i] = last;
		INTEGER(VECTOR_ELT(ans, 3))[i] = ier;
		akeep = par.alpha;
		bkeep = par.beta;
	}

    setAttrib(ans, R_NamesSymbol, ansnames);
    UNPROTECT(2);
    return ans;
}


