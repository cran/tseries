/* Copyright (C) 1997-1999  Adrian Trapletti
  
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.
  
   You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   Support functions for time series utilities */


#include <math.h>
#include "S.h"


extern void muin2ser_ ();


void R_intgrt_vec (double* x, double* y, int* lag, int* n)
{
  int i;

  for (i=*lag; i<*lag+*n; i++)
    y[i] = x[i-*lag]+y[i-*lag];
}

void R_embed_vec (double* x, double* y, int* dim, int* n)
{
  int i, j, nr, nc;

  nr = *n-*dim+1;
  nc = *dim;
  for (i=0; i<nc; i++)
    for (j=0; j<nr; j++)  
      y[j+nr*(nc-i-1)] = x[i+j];
}

void R_pp_sum (double* u, int* n, int* l, double* sum)
{
  int i, j;
  double tmp1, tmp2;
  
  tmp1 = 0.0;
  for (i=1; i<=(*l); i++)
  {
    tmp2 = 0.0;
    for (j=i; j<(*n); j++)  
    {
      tmp2 += u[j]*u[j-i];
    }
    tmp2 *= 1.0-((double)i/((double)(*l)+1.0));
    tmp1 += tmp2;
  }
  tmp1 /= (double)(*n);
  tmp1 *= 2.0;
  (*sum) += tmp1;
}

void R_amif (double *x, int *lx, double *inf, int *k, int *maxbit, 
	     double *confidence, int *norm, int *trace)
{
  double information, freq;
  double *x0, *y0;
  int *intx, *inty, *s, *q, *q_unsort, *indices_x, *indices_y, *position;
  int ixmax, i, j;

  freq = 1.0;
  (*lx) -= (*k);
  x0 = Calloc ((*lx), double); 
  y0 = Calloc ((*lx), double);
  intx = Calloc ((*lx), int);
  inty = Calloc ((*lx), int);
  s = Calloc ((*lx), int);
  q = Calloc ((*lx), int);
  q_unsort = Calloc ((*lx), int);
  indices_x = Calloc ((*lx), int);
  indices_y = Calloc ((*lx), int);
  ixmax = (int) pow (2.0, (double) (*maxbit));
  position = Calloc (ixmax+1, int);
  for (i=0; i<=(*k); i++)
  {
    for (j=0; j<(*lx); j++)
    {
      x0[j] = x[j]; y0[j] = x[j+i];
    }
    muin2ser_ (&information, lx, &freq, x0, y0, intx, inty, s, q, q_unsort, indices_x, indices_y, 
	       position, maxbit, confidence, trace);
    if (*norm) inf[i] = sqrt(1.0-exp((-2.0)*information));
    else inf[i] = information;
  }
  if (*norm) inf[0] = 1.0;
  Free (x0); Free (y0);
  Free (intx); Free (inty);
  Free (s); Free (q);
  Free (q_unsort);
  Free (indices_x); Free (indices_y);
  Free (position);
}

void R_quad_map (double *x, double *xi, double *a, int *n)
{
  int i;

  x[0] = (*xi);
  for (i=1; i<(*n); i++)
    x[i] = (*a)*(1-x[i-1])*x[i-1];
}

void R_Durbin (double* cor, double* pacf, int* lag)
{
  int l, j;
  double sum1, sum2;

  pacf[0] = cor[0];
  for (l=2; l<=(*lag); l++)
  {
    sum1 = sum2 = 0.0;
    for (j=1; j<=l-1; j++)
    {
      sum1 += pacf[l-2+(*lag)*(j-1)]*cor[l-j-1];
      sum2 += pacf[l-2+(*lag)*(j-1)]*cor[j-1];
    }
    pacf[l-1+(*lag)*(l-1)] = (cor[l-1]-sum1)/(1-sum2);
    for (j=1; j<=l-1; j++)
      pacf[l-1+(*lag)*(j-1)] = pacf[l-2+(*lag)*(j-1)]
	-pacf[l-1+(*lag)*(l-1)]*pacf[l-2+(*lag)*(l-j-1)];
  }
}
