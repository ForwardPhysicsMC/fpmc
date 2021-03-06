/************************************************
*    CompHEP version 4.5.1   *
*------------------------------------------------
* Copyright (C) 2001-2009, CompHEP Collaboration*
************************************************/
#include<math.h>
#define real double
#include"out_int.h"
#include"out_ext.h"

namespace anom_aaaa {

extern real *Q0,*Q1,*Q2;
extern double va[3];
DNN  d_1;
int d_1(real * momenta)
{int I,err=0;
real s0max=0;
 for(I=0;I<nin_;I++) s0max+=momenta[4*I];
s0max=computer_eps*s0max*s0max;
return err;
}

} //namespace anom_aaaa
