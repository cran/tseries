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

   time series bootstrap functions */


#include <math.h>
#include <R.h>


static int geodev (double p)
     /* Return geometric distributed random deviate with p(i) = p*(1-p)^i,
	where 0<p<1 and i=0,1,..
	See Fishman, G. S., 1996, Monte Carlo (Springer, NY, Berlin, Heidelberg), 221. */
{
  double b, y, z;
  
  b = (-1.0)/log(1.0-p);
  y = exp_rand();
  z = b*y;
  return (int)z;
}

static int disuni (int n)
     /* Return discrete uniform distributed random deviate on {1,..,n} */
{
  double temp;
  
  temp = unif_rand()*(double)n+1.0;
  return (int)temp;
}

static int WRAP (int i, int n)
     /* WRAP (i=..,-n+1,..,-1,0,1,..,n,n+1,.., n) = ..,1,..,n-1,n,1,..,n,1,.. */
{
  if (i < 1)
    return i%n+n;
  else if (i > n)
    return (i-1)%n+1;
  else
    return i;
}

static void StatBoot (double x[], double xBoot[], int n, double p)
     /* Generates a bootstrap sample xBoot[1..n] from x[1..n] 
	according to the stationary bootstrap resampling scheme
	(Politis, D. N. and Romano, J. P., 1994, The Stationary Bootstrap, 
	J. Amer. Statist. Assoc. 89, 1303-1313).
	Input is n, x[1..n], xBoot[1..n] and p, the parameter for 
	the geometric distribution, output xBoot[1..n]. */
{
  int i, j, I, L;
  
  i = 1;
  while (i <= n)
  {
    I = disuni(n);
    L = geodev(p);
    j = 0;
    while ((j < L) && (i <= n))
    {
      xBoot[i] = x[WRAP(I+j,n)];
      i++;
      j++;
    }
  }
}

static void BlockBoot (double x[], double xBoot[], int n, int L)
     /* Generates a bootstrap sample xBoot[1..n] from x[1..n] 
	according to the blockwise bootstrap resampling scheme
	(Kuensch, H. R., 1989, The Jackknife and the Bootstrap 
	for General Stationary Observations, Ann. Stat. 17, 1217-1241).
	Input is n, x[1..n], xBoot[1..n] and L, the blocklength, 
	output xBoot[1..n]. */
{
  int i, j, I;
  
  i = 1;
  while (i <= n)
  {
    I = disuni(n-L+1);
    j = 0;
    while ((j < L) && (i <= n))
    {
      xBoot[i] = x[I+j];
      i++;
      j++;
    }
  }
}

void tseries_boot (double *x, double *xb, int *n, double *b, int *type)
{
  GetRNGstate();
  if (*type == 0) StatBoot (x-1, xb-1, *n, *b);
  else if (*type == 1) BlockBoot (x-1, xb-1, *n, (int)*b);
  else error ("this type of bootstrap is not yet implemented\n");
  PutRNGstate();
}



