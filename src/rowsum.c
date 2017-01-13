
void rowsum_kahan(int *len,  double *x, double *group)
{
     register int i, k;
     int  nrow,  newrow;
     double  tgrp,	dummy;
     long double  sum, c,t,y;
     
     nrow = *len;
     
     dummy =0;
     for (i=0; i<nrow; i++) if (group[i] < dummy) dummy = group[i];
     dummy = (dummy/2) -1;    /*no group uses this number */
     
     newrow =0;
     for (i=0; i<nrow; i++) {
	  if (group[i] > dummy) {
	       tgrp = group[i];
	       sum =0;
	       c=0;
	       for (k=i; k<nrow; k++) {
		    if (group[k] == tgrp) {
			 y = x[k]-c;
			 t = sum+y;
			 c = (t-sum) -y;
			 sum = t;
		    }
	       }    
	       x[newrow] = sum;
	       for (k=i; k<nrow; k++){
		    if (group[k] == tgrp) group[k] = dummy;
	       }
	       newrow++;
	  }
     }
     *len = newrow;
     return;
}
