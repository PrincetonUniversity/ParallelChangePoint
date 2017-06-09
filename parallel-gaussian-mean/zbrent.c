/***************************************************************************
	zbrent.c  -  Recursively find critical value
	-------------------
	email : hawyang@princeton.edu

	See:
	L. Watkins and H. Yang, J Phys Chem B (2005), 109 1):617-28.
		doi:10.1021/jp0467548
	Marc Noe, Ann. Math. Statist (1972), 43(1):58-64.
 ***************************************************************************/
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define ITMAX 100
#define EPS 3.0e-8

double zbrent(func,x1,x2,tol,p1,p2)
double x1,x2,tol,p1,p2;
double (*func)();	/* ANSI: double (*func)(double); */
{
	int iter;
	double a=x1,b=x2,c,d,e,min1,min2;
	double fa=(*func)(a,p1,p2),fb=(*func)(b,p1,p2),fc,p,q,r,s,tol1,xm;

	if (fb*fa > 0.0) {
		printf("Root must be bracketed in ZBRENT\n");
		exit(1);
	}
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {

		if (fb*fc > 0.0) {
			c=a;
			fc=fa;
			e=d=b-a;
		}

		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) 
			return b;

		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			}
			else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  
				q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			}
			else {
				d=xm;
				e=d;
			}
		}
		else {
			d=xm;
			e=d;
		}

		a=b;
		fa=fb;

		if (fabs(d) > tol1)
			b += d;
		else
			b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));

		fb=(*func)(b,p1,p2);
	} 

	printf("Maximum number of iterations exceeded in ZBRENT\n");
	exit(1);
}

#undef ITMAX
#undef EPS
