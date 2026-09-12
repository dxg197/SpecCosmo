#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

CCTK_REAL determinant(CCTK_REAL a[4][4], CCTK_REAL k);
void cofactor(CCTK_REAL num[4][4], CCTK_REAL inverse[4][4], CCTK_REAL f);
void transpose(CCTK_REAL num[4][4], CCTK_REAL fac[4][4], CCTK_REAL inverse[4][4], CCTK_REAL r);
 
/*For calculating Determinant of the Matrix */
CCTK_REAL determinant(CCTK_REAL a[4][4], CCTK_REAL k)
{
  CCTK_REAL s = 1, det = 0, b[4][4];
  int i, j, m, n, c;
  if (k == 1)
    {
     return (a[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] = a[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (a[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }
 
    return (det);
}
 
void cofactor(CCTK_REAL num[4][4], CCTK_REAL inverse[4][4], CCTK_REAL f)
{
 CCTK_REAL b[4][4], fac[4][4];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  transpose(num, fac, inverse, f);
}
/*Finding transpose of matrix*/ 
void transpose(CCTK_REAL num[4][4], CCTK_REAL fac[4][4], CCTK_REAL inverse[4][4], CCTK_REAL r)
{
  int i, j;
  CCTK_REAL b[4][4], d;
 
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num,r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        inverse[i][j] = b[i][j] / d;
        }
    }
}