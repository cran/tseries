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
#include "R.h"


extern void F77_SYMBOL(muin2ser) ();

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
    F77_SYMBOL(muin2ser) (&information, lx, &freq, x0, y0, intx, inty, 
			  s, q, q_unsort, indices_x, indices_y, 
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

