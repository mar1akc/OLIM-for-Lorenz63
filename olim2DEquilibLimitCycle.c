// Ordered Line Integral Method - MIDPOINT (OLIM-MID) in deformed polar coordinates
// for computing the quasi-potential with respect to the stable equilibrium  
// up to unstable limit cycle 
// Copyright: Maria Cameron, January 14, 2018

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.141592653589793
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-12
// Parameters for the Lorenz system
#define BETA 2.6666666666666666
#define RHO 20.0
#define SIGMA 10.0
#define SSF 4


struct myindex {
	int il; // longitude index
	int ir; // radial index
};	

struct myvector {
  double x;  
  double y;
  double z;		
};

struct mysol {
  double g;
  char c;
};  



int main(void);
struct myvector myfield(struct myvector x); /* B */
double angle(double x,double y);
double length(double x,double y);
void param(void);
void init_strange_attractor(void);
void olim(void);
struct mysol triangle_update(int ind,int ind0,int ind1);
double one_pt_update(int ind,int ind0);			   
void addtree(int ind); /* adds a node to the binary tree
                                 of the "considered" points */ 
void updatetree(int ind); /* updates the binary tree */
void deltree(void); /* deletes the root of the binary tree */
double dot_product(struct myvector a,struct myvector b); 
struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b);
struct myvector vec_difference(struct myvector v1,struct myvector v2);
struct myvector vec_sum(struct myvector v1,struct myvector v2);
struct myvector getpoint(int ind);
double length_vec(struct myvector x);
struct mysol hybrid_nonlin_solver(double a,double b,double u0,double u1,
			struct myvector x0,struct myvector x1,
			struct myvector bm0,struct myvector bm1,struct myvector x);
double myfun(double s,double u1mu0,struct myvector xmx0,struct myvector x1mx0,
			struct myvector bm0,struct myvector bm1,
			struct myvector b1mb0m,struct myvector x);
double geometric_action_line(struct myvector x0, struct myvector x1);
int mod(int x,int m);
int getneighbors(int ind,int Kr,int kl,int *nei);
int getneighbors_simple(int ind,struct myindex sar,int Kr,int kl,int *nei,int nneib);
struct myindex get_index(int ind);

/***************************************/

int Nrad,Nloop,Kr,Kl,nr1,nl1,Nmesh; 
int NCURVE;
int count=0; /* # of considered points */
char *ms,*stype; /* 0 = 'Unknown', 1 = 'Considered', 2 = "in Accepted Front", 3 = "Accepted" */
double *g; /* function to be computed */
int *pos; /* pos(index of mesh pt) = position in binary tree */
int *tree; /* tree(position in the tree) = index of mesh pt */

// variables for the potential 
double *xmesh,*ymesh,*zmesh;
int PRINT = 0;

/**************************************/
struct myvector myfield(struct myvector x) {
  struct myvector b;
  double aux;

	// Lorenz'63
	b.x = SIGMA*(x.y - x.x);
	b.y = x.x*(RHO - x.z) - x.y;
	b.z = x.x*x.y - BETA*x.z;

  return b;
}

/*************************************/

void param() {
  int i,j,ind,n; 
  double t; 
  FILE *fid;

  printf("in param()\n");
  
  ms = (char *)malloc(Nmesh*sizeof(char)); 
  stype = (char *)malloc(Nmesh*sizeof(char)); 
  pos = (int *)malloc((Nmesh + 1)*sizeof(int)); 
  tree = (int *)malloc((Nmesh + 1)*sizeof(int)); 
  g = (double *)malloc(Nmesh*sizeof(double));
  
  nl1 = Nloop - 1;
  nr1 = Nrad - 1;
  printf("Nrad = %i, Nloop = %i, Kr = %i, Kl = %i, Nmesh = %i\n",Nrad,Nloop,Kr,Kl,Nmesh);
  
  for( i = 0; i < Nmesh; i++ ) {
  		ms[i] = 0;
  		g[i] = INFTY;
  		stype[i] = 'n';	  
  }
}

/************************************/

void init_strange_attractor() {
  int ind;
  int Ninit = 0;

  for( ind = 0; ind < Nloop; ind++ ) {
	 g[ind] = 0.0;
	 ms[ind] = 1;
	 stype[ind] = 0;
	 addtree(ind);
	 Ninit++;
 }
 printf("Ninit: %i\n",Ninit);
printf("count = %i\n",count); 
} 
 


/**********************************************/
/*** ordered line integral method ***/

void olim(void) {
  int k,m,ind,ind0,ind1,indupdate,imin,n;
  int mycount = 0; /* counter for accepted points */
  double gmin,gtemp,gold;
  int Naf, AFneib[8], Nc, NCneib[8]; /* accepted front neighbors and new considered of the newly accepted point */
  struct mysol sol; 
  struct myvector vec,b,bm0,bm1,b0,b1,v0,v1,vnew; 
  double s, s1;
  char pr = 'n'; // a logic variable: if pr == 'y', solution will be printed out in the terminal window
  int nneib, nnaux, fneib, *fnei, *nei, *naux;
//   char do_flag = 'y';	
  FILE *fid,*fs;
  int fneibmax = (2*Kr+1)*(2*Kl+1)-1, nneibmax = 8;
 struct myindex ilw;

  printf("in olim()\n");
  
//   fs = fopen("fstink1.txt","w");


  fnei = (int *)malloc(fneibmax*sizeof(int));	// far neighbors
  nei = (int *)malloc(nneibmax*sizeof(int));	// nearest neighbors
  naux = (int *)malloc(nneibmax*sizeof(int));	// nearest neighbors

  while( count > 0 ) {
  
//   					for( k = 1; k <= count; k++ )  if( pos[tree[k]] != k ) {
// 						printf("pos[%i] = %i, tree[%i] = %i\n",tree[k],pos[tree[k]],k,tree[k]);
// 						exit(1);
// 					}		

  
    ind = tree[1];
    ilw = get_index(ind);
    vnew = getpoint(ind);
    /* x and y of the newly accepted point */
    ms[ind] = 2;
    deltree();
	mycount++;
	
	// FINISH if the top boundary is reached
// 	if( ilw.iw == nw1 ) {
// 		printf("Boundary Nrad is reached: ind = %i, iw = %i, il = %i, g = %.4e\n",ind,ilw.iw,ilw.il,g[ind]);
// 		break;
// 	}	   

 	if( mycount%10000 == 0 ) 
	printf("mycount = %i, Accept %i, g = %.4e\n",mycount,ind,g[ind]);

// 	if( ind < Nloop ) {
// 		printf("mycount = %i: Accept: %d, u = %.4e, stype = %i\n",mycount,ind,g[ind],stype[ind]);
// 	}
	// CONTINUE
    /* Inspect the neighbors of the new Accepted point */
    nneib = getneighbors(ind,1,1,nei);
    if( nneib > nneibmax ) {
    	printf("ind = %i, nneib = %i\n",ind,nneib); 
    	exit(1);
    }	
    Naf = 0;
    Nc = 0;
    for( k=0; k<nneib; k++ ) {
		  ind1 = *(nei + k); // neighbor of the new Accepted point
// 		  printf("ind1 = %i, ms = %i\n",ind1,ms[ind1]);
		  // update AcceptedFront
		  if( ms[ind1] == 2 ) {
				m = 0;
				nnaux = getneighbors(ind1,1,1,naux);
				if( nnaux > nneibmax ) {
					printf("ind1 = %i, nnaux = %i\n",ind1,nnaux); 
					exit(1);
				}	
				for( n = 0; n < nnaux; n++ ) {
				   ind0 = *(naux + n);
				   if( ms[ind0] < 2 ) m++;
				}   
				if( m == 0 ) { /* ind1 has no considered neighbors */
				  ms[ind1] = 3;
				}
				else {
				  AFneib[Naf] = ind1;
				  Naf++;
				}  
		  }
		  else if( ms[ind1] == 0 ) { // the neighbor ind1 will be a new Considered point
		    	vec = getpoint(ind);
				NCneib[Nc] = ind1;
				Nc++;
		  }  
	}
	
	
	/* update considered points */
	fneib = getneighbors(ind,Kr,Kl,fnei);
	if( fneib > fneibmax ) {
		printf("ind = %i, fneib = %i\n",ind,fneib);
		exit(1);
	}
		
	for( k = 0; k < fneib; k++ ) {
		indupdate = *(fnei + k);
// 		if( pr == 'y' ) printf("k = %i, indupdate = %i\n",k,indupdate);
		if( ms[indupdate] == 1 && stype[indupdate] != 0 ) { 
			// existing Considered neighbors not lying on the strange attractor
			gold = g[indupdate];
			gtemp  = one_pt_update(indupdate,ind);
			vec = getpoint(indupdate);
			if( gtemp < g[indupdate] ) {
				g[indupdate] = gtemp;
				stype[indupdate] = 1;
// 		   		   if ( indupdate == SUSPECT) {
// 		        printf("1-pt-update: ind = %i, indupdate = %i: ind0 = %i, g0 = %.4e, g = %.4e\n",
// 		        	ind, indupdate, ind0,g[ind0],g[indupdate]);
// 		  }							
			}
			bm0 = myfield(vec_lin_comb(vnew,vec,0.5,0.5)); 
			for( m = 0; m < Naf; m++ ) {
			  ind1 = AFneib[m];
			  v1 = getpoint(ind1);
			  bm1 = myfield(vec_lin_comb(v1,vec,0.5,0.5)); 
			  sol = hybrid_nonlin_solver(0.0,1.0,g[ind],g[ind1],vnew,v1,bm0,bm1,vec);
			  
			  PRINT = 0;
			  
			  if( sol.c == 'y' ) {
				  s = sol.g;
				  s1 = 1.0 - s;	
				  gtemp = s1*g[ind] + s*g[ind1] + geometric_action_line(vec_lin_comb(vnew,v1,s1,s),vec);	
				  if( gtemp < g[indupdate] ) {
					g[indupdate] = gtemp;
					stype[indupdate] = 2;
// 		   		   if ( indupdate == SUSPECT) {
// 		        printf("triangle-update: ind = %i, indupdate = %i:  g = %.4e\n",
// 		        	ind, indupdate,g[indupdate]);
// 		  }							
				  }	
			  }	
			} // end for( m = 0; m < Naf; m++ ) 
			if( gold > g[indupdate] ) {
			  updatetree(indupdate);
			}   
		} 
	}
		
     /* Shift Unknown neighbors of the new Accepted point to Considered and compute values at them */ 			  
	 for( m = 0; m < Nc; m++ ) { /* for all points that switched from unknown to considered */
		   indupdate = NCneib[m];
		   vec = getpoint(indupdate);
		   gmin = INFTY;
		   imin = ind;
		   fneib = getneighbors(indupdate,Kr,Kl,fnei);
		   	if( fneib > fneibmax ) {
				printf("indupdate = %i, fneib = %i\n",indupdate,fneib);
				exit(1);
			}

		   for( n = 0; n < fneib; n++ )  {
			   ind0 = *(fnei + n);
				 /* compute tentative values using poins of the accepted front or close accepted poonts */
				 if( ms[ind0] == 2 ) {//|| (ms[ind0]==3 && labs(i-i0)<1.5 && labs(j-j0)<1.5) ) {
					 gtemp = one_pt_update(indupdate,ind0);
					 if( gtemp < gmin ) {
						 gmin = gtemp;
						 imin = ind0;
					 }
				 }
		   }
		   ind0 = imin;	 
		   g[indupdate] = gmin;
		   stype[indupdate] = 1;
		   v0 = getpoint(ind0);
		   bm0 = myfield(vec_lin_comb(v0,vec,0.5,0.5)); 
// 		   		   if ( indupdate == SUSPECT) {
// 		        printf("New: 1-pt-update: ind = %i, indupdate = %i: ind0 = %i, g0 = %.4e, gmin = %.4e\n",
// 		        	ind, indupdate, ind0,g[ind0],g[indupdate]);
// 		  }							

		   nneib = getneighbors(ind0,1,1,nei);
    		if( nneib > nneibmax ) {
    			printf("ind0 = %i, nneib = %i\n",ind0,nneib); 
    			exit(1);
    		}	
		   for( k = 0; k < nneib; k++ ) {
			 ind1 = *(nei + k);
			 if( ms[ind1] == 2 ) {
				 v1 = getpoint(ind1);
				 bm1 = myfield(vec_lin_comb(v1,vec,0.5,0.5));				     
				 sol = hybrid_nonlin_solver(0.0,1.0,g[ind0],g[ind1],v0,v1,bm0,bm1,vec);
				  if( sol.c == 'y' ) {
					  s = sol.g;
					  s1 = 1.0 - s;	
					  gtemp = s1*g[ind0] + s*g[ind1] + geometric_action_line(vec_lin_comb(v0,v1,s1,s),vec);
					  if( gtemp < g[indupdate] ) {
						g[indupdate] = gtemp;
						stype[indupdate] = 2;
// 		   		   if ( indupdate == SUSPECT) {
// 		        printf("New: triangle-update: ind = %i, indupdate = %i: g = %.4e\n",
// 		        	ind, indupdate, g[indupdate]);
// 		  }							
					  }	
				  }				  
			  } /* end if( ms[ind1] == 2 ) */
		  }	      
		  addtree(indupdate);
		  ms[indupdate] = 1;
	} /* end for( m = 0; m < Nc; m++ ) */
  } /* end while ( count > 0 ) */
}

  
/*********************************************/

struct myindex get_index(int ind) {
	struct myindex mdx;
		mdx.ir = ind/Nloop;
		mdx.il = ind%Nloop;
	return mdx;
}

/*********************************************/


double one_pt_update(int ind,int ind0) {
  struct myvector x,x0,l,b,bb;
  double gtemp,ab,len;
  
  x = getpoint(ind);
  x0 = getpoint(ind0);
//    printf("ind0 = %i, ind = %i:",ind0,ind);
  gtemp = g[ind0] + geometric_action_line(x0,x);

  return gtemp;
}
/*-------------*/

double geometric_action_line(struct myvector x0, struct myvector x1) {
  struct myvector dx,v0,vm,v1,b0,bm,b1;
  
  vm = vec_difference(x1,x0);
  bm = myfield(vec_lin_comb(x0,x1,0.5,0.5));

  return  length_vec(bm)*length_vec(vm) - dot_product(bm,vm);
}


/*-------------*/  
  
struct myvector getpoint(int ind) {
  struct myvector x;
  
	x.x = xmesh[ind];
	x.y = ymesh[ind];
	x.z = zmesh[ind];
  return x;
}

/*----------------*/

int getneighbors(int ind,int kr,int kl,int *nei) {
	int mr,ml;
	int i,m,nneib = 0;
    struct myindex ilr; 
  	
    ilr = get_index(ind);
	for( i = -kr; i <= kr; i++ ) {
		mr = ilr.ir + i;
		if( mr >= 0 && mr <= nr1 ) {			
			for( m = -kl; m <= kl; m++ ) {
			    if( m != 0 || i != 0 ) {
			    	ml = mod(ilr.il + m,Nloop);
					*(nei + nneib) = mr*Nloop + ml;		
					nneib++;
				}
			}
		}
	}
	
	return nneib;
		
}

			
/***** N o n l i n e a r   1D    s o l v e r *****/


struct mysol hybrid_nonlin_solver(double a,double b,double u0,double u1,
			struct myvector x0,struct myvector x1,
			struct myvector bm0,struct myvector bm1,struct myvector x) {
	double c, fa, fb, fc, d, fd, dd, df, dm, ds, t, u1mu0;
	struct myvector x1mx0,b1mb0m,xmx0;
	struct mysol sol;
	double NONLIN_TOL = 1.0e-6;
	int iter = 0, itermax = 100;

	x1mx0 = vec_difference(x1,x0);
	b1mb0m = vec_difference(bm1,bm0);
	xmx0 = vec_difference(x,x0);
	u1mu0 = u1 - u0;
	
	c = a;
	fa = myfun(a,u1mu0,xmx0,x1mx0,bm0,bm1,b1mb0m,x); 
	fb = myfun(b,u1mu0,xmx0,x1mx0,bm0,bm1,b1mb0m,x); 
	
// 	if( PRINT == 1 ) printf("(a,fa) = (%.3e,%.3e), (b,fb) = (%.3e,%.3e)\n",a,fa,b,fb);
	
	
	fc = fa;
	if( (fa > 0 && fb > 0 ) || (fa < 0 && fb < 0) ) {
//	 root is not bracketed 
		sol.c = 'n';
		sol.g = INFTY;
		return sol;
	}
	while( iter < itermax ) {
		if( fabs(fc) < fabs(fb) ) { 
			t = c; c = b; b = t;
			t = fc; fc = fb; fb = t;
			a = c; fa = fc;
		}		
		if( fabs(b - c) < NONLIN_TOL ) break;		
		dm = 0.5*(c - b);
		df = fa - fb;

		if( fabs(df) < NONLIN_TOL ) ds = dm;
		else ds = -fb*(a - b)/df;
		if( (ds > 0 && dm < 0) || (ds < 0 && dm > 0) || (fabs(ds) > fabs(dm)) ) dd = dm;
		else dd = ds;
		
		if( fabs(dd) < NONLIN_TOL ) dd = 0.5*sgn(dm)*NONLIN_TOL;

		d = b + dd;
		fd = myfun(d,u1mu0,xmx0,x1mx0,bm0,bm1,b1mb0m,x); 
		
// 		if( PRINT == 1 ) printf("(d,fd) = (%.3e,%.3e)\n",d,fd);

		
		if( fabs(fd) < NONLIN_TOL ) {
			b = d;
			break;
		}
		a = b; b = d; fa = fb; fb = fd;
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
			c = a; fc = fa;
		}
		iter++;
	}
	sol.c = 'y';
	sol.g = b;
	
	return sol;
}


/*--------------*/

double myfun(double s,double u1mu0,struct myvector xmx0,struct myvector x1mx0,
			struct myvector bm0,struct myvector bm1,
			struct myvector b1mb0m,struct myvector x) {
	double s1 = 1.0 - s,lbms,dlppds,lpp,dlbms,dp1,dp2;
	struct myvector pp,bms;
	
	pp = vec_lin_comb(xmx0,x1mx0,1,-s); // psi'		
	bms = vec_lin_comb(bm0,bm1,s1,s);  // bms
	lbms = length_vec(bms); // ||b_{ms}||
	lpp = length_vec(pp); // ||psi'||
	dlppds = -dot_product(pp,x1mx0)/lpp; // d ||psi'|| / ds 
	dlbms = dot_product(bms,b1mb0m)/lbms; // d ||bms|| / ds
	dp1 = -dot_product(bms,x1mx0);
	dp2 = dot_product(b1mb0m,pp);
	
	return u1mu0 + lbms*dlppds + dlbms*lpp - dp1 - dp2;
	
}
			



/**********************************************/
/*** linear algebra ***/

double length_vec(struct myvector x) {
  return sqrt(x.x*x.x + x.y*x.y + x.z*x.z);
}


double dot_product(struct myvector a,struct myvector b) {
   return a.x*b.x + a.y*b.y + a.z*b.z;
}


struct myvector vec_lin_comb(struct myvector x1,struct myvector x2,double a,double b) {
	struct myvector v;
	
	v.x = a*x1.x + b*x2.x;
	v.y = a*x1.y + b*x2.y;
	v.z = a*x1.z + b*x2.z;
	return v;
}
	
struct myvector vec_difference(struct myvector x1,struct myvector x2) {
	struct myvector v;
	
	v.x = x1.x - x2.x;
	v.y = x1.y - x2.y;
	v.z = x1.z - x2.z;
	return v;
}
	
struct myvector vec_sum(struct myvector x1,struct myvector x2) {
	struct myvector v;
	
	v.x = x1.x + x2.x;
	v.y = x1.y + x2.y;
	v.z = x1.z + x2.z;
	return v;
}

int mod(int x,int m) {

	while( x < 0 ) x += m;
	while( x >= m ) x -= m;

    return x;
}



/****************************************/
/**************************************************************/
/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(int ind) {
  int loc, ptemp;
  int indp, indc;
  char ch;
//   int SUSPECT = 10371;

//  if( ind == SUSPECT ) printf("addtree(%i)\n",ind);
  count++;
  if( count > Nmesh ) {
  	printf("Nmesh = %i, count = %i\n",Nmesh,count);
  	exit(1);
  }	
  tree[count]=ind;
  pos[ind]=count;
  if( count > 1 ) {
    loc=count;
    indc=tree[loc];
    indp=tree[loc/2];
    ch=( g[indc] < g[indp] ) ? 'y' : 'n';
    while( ch == 'y' ) {
      ptemp=pos[indc];
      pos[indc]=pos[indp];
      tree[loc/2]=indc;
      pos[indp]=ptemp;
      tree[loc]=indp;
      loc=loc/2;
      if( loc > 1 ) {
        indc=tree[loc];
        indp=tree[loc/2];
        ch=( g[indc] < g[indp] ) ? 'y' : 'n';
      }
      else ch='n';
    }
  }
}

/*------------------------------------------------------------------*/

void updatetree(int ind) {
  int loc, lcc;
  double g0;
//   int SUSPECT = 4353;
//   char pr = 'n';

//  if( ind == SUSPECT ) printf("updatetree(%i)\n",ind);

  g0 = g[ind];
  loc = pos[ind];
//   if( ind == SUSPECT ) {
//   		printf("updatetree(%i): g0 = %.4e, pos[%i] = %i\n",ind,g0,ind,pos[ind]);
//   		printf("loc/2 = %i\n",loc/2);
//   		printf("tree[loc/2] = %i\n",tree[loc/2]);
//   }		
  while( loc > 1 && g0 < g[tree[loc/2]] ) {
//     if( ind == SUSPECT ) printf("go up: loc = %i, loc/2 = %i\n",loc,loc/2);
    tree[loc]=tree[loc/2];
    pos[tree[loc]]=loc;
    loc=loc/2;
    tree[loc]=ind;
    pos[tree[loc]]=loc;
//      if( ind == SUSPECT ) printf("updatetree(%i): loc = %i\n",ind,loc);
  } 
  
  lcc = count;
  while( (loc*2 <= count && g0 > g[tree[loc*2]]) || (loc*2+1 <= count && g0 > g[tree[loc*2+1]]) )  {
    lcc=( loc*2+1 <=count && g[tree[loc*2+1]] < g[tree[loc*2]] ) ? loc*2+1 : loc*2;
    tree[loc]=tree[lcc];
    pos[tree[loc]]=loc;
    loc=lcc;
    tree[loc]=ind; 
    pos[tree[loc]]=loc;

  }
}

/*---------------------------------------------------------------------*/


/* deletes root of the binary tree */
void deltree() {
  int loc, ptemp, ind, lcc, ic, ic1, ic2, mind;
  char chd, ch='n';
  char pr = 'n';

  mind=tree[1];
  
  pos[tree[1]] = 0;
  tree[1] = tree[count];
  pos[tree[1]] = 1;
  count--;
  loc=1;
  ind=tree[1];
  lcc=2*loc;
  if( lcc < count )  {
	  ic1 = tree[lcc];
	  ic2 = tree[lcc+1];
	  if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
		  if( (g[ic1]) <= (g[ic2]) )  {
			chd = 'l';
			ic = ic1;
		  }
		  else {
			chd = 'r';
			ic = ic2;
			lcc++;
		  }
	  }
	  else chd = 'n';
  }
  else if( lcc == count ) {
	  ic = tree[lcc];
	  if( (g[ind]) > (g[ic]) ) {
		chd='l'; 
	  	if(ch=='y') 
	  	printf("left\n");
	  }
	  else chd = 'n';
  }
  else chd='n';
  

  while( chd != 'n' ) {    
  // swap parent and child
    ptemp = pos[ind];
    pos[ind] = pos[ic];
    tree[loc] = ic;
    pos[ic] = ptemp;
    tree[lcc] = ind;
    loc = lcc;        
    lcc=2*loc;
    if( lcc < count )  {
      ic1=tree[lcc];
      ic2=tree[lcc+1];
      if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
        if( (g[ic1]) <= (g[ic2]) )  {
          chd='l';
	      ic=ic1;
        }
        else {
          chd='r';
	      ic=ic2;
	      lcc++;
        }
      }
      else chd='n';
    }
    else if( lcc == count ) {
      ic=tree[lcc];
      if(ch=='y') printf("child: loc(%i)=%i, t1=%.12e\n",ic1,lcc,g[ic1]);
      if( (g[ind]) > (g[ic]) ) { chd='l';if(ch=='y') printf("left\n");}
      else chd='n';
    }
    else chd='n';
  } /* end while( chd != 'n' ) */
  
}


/********************************************************/		    
/*** main ***/

 int main() {
  int jr,jl,ind,k,nprint,i,j; 
  double tmp;
  clock_t CPUbegin;
  double cpu;
  double umin = INFTY,umax = 0.0 ,umean = 0.0;
  struct myvector vec;
  struct myindex sar;
  FILE *fid,*fid1,*fid2,*fid3,*fid4;
  char fname[100];
  
  // read file with parameters
  sprintf(fname,"Mesh2Ddata_rho%.2f.txt",RHO);
  fid = fopen(fname,"r");
  if( fid == NULL ) {
  	  printf("Cannot find the file %s\n",fname);
  	  exit(1);
  }
  fscanf(fid,"%i\n%i\n",&Nloop,&Nrad);
  Nloop--; // take into account that the mesh is periodic
  fclose(fid);
  // set up update radii
  Kl = Nloop/40; // angular
  Kr = Nrad/40;  // radial
  Nmesh = Nrad*Nloop;
  xmesh = (double *)malloc(Nmesh*sizeof(double));
  ymesh = (double *)malloc(Nmesh*sizeof(double));
  zmesh = (double *)malloc(Nmesh*sizeof(double));
  
  // read mesh data
  sprintf(fname,"Mesh2Drho%.2f.txt",RHO);
  // format: x Nrad by Nloop  
  // y Nrad by Nloop
  // z Nrad by Nloop
  fid = fopen(fname,"r");
  if( fid == NULL ) {
  	  printf("Cannot find the file %s\n",fname);
  	  exit(1);
  }
  for( jr = 0; jr < Nrad; jr++ ) {
  	  for( jl = 0; jl < Nloop; jl++ ) {
  	  	  fscanf(fid,"%le\t",&tmp);
  	  	  ind = jl + jr*Nloop;
  	  	  xmesh[ind] = tmp;
  	  }
  } 
  for( jr = 0; jr < Nrad; jr++ ) {
  	  for( jl = 0; jl < Nloop; jl++ ) {
  	  	  fscanf(fid,"%le\t",&tmp);
  	  	  ind = jl + jr*Nloop;
  	  	  ymesh[ind] = tmp;
  	  }
  } 
  for( jr = 0; jr < Nrad; jr++ ) {
  	  for( jl = 0; jl < Nloop; jl++ ) {
  	  	  fscanf(fid,"%le\t",&tmp);
  	  	  ind = jl + jr*Nloop;
  	  	  zmesh[ind] = tmp;
  	  }
  } 
  fclose(fid);	  	  
  	  	  
  	  	  	  
  param();
  init_strange_attractor();
  
  CPUbegin=clock();
  
  olim();
  
  printf("Returned from olim()\n");
  
  cpu=(clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
  printf("cputime of olim() = %g\n",cpu);
  
  // OUTPUT FILES
  sprintf(fname,"Qpot2Drho%.2f.txt",RHO); // output file with the quasipotential

  fid = fopen(fname,"w");

  for( jr = 0; jr < Nrad; jr+=SSF ) {
  	 for( jl = 0; jl < Nloop; jl+=SSF ) {	
  	 	ind = jr*Nloop + jl;
	 	if( ms[ind] < 2 ) g[ind] = INFTY;			
		fprintf(fid,"%.6e\t",g[ind]);
	 }
	 fprintf(fid,"\n");	
  }
  fclose(fid);
  k = 0;
  for( jl = 0; jl < Nloop; jl++ ) {
  	ind = nr1*Nloop + jl;
  	umin = min(umin,g[ind]);
  	umax = max(umax,g[ind]);
  	umean += g[ind];
  	k++;
  }
  umean /= k;
  	
  printf("Nrad = %i, Nloop = %i, Kr = %i, Kl = %i, subsampling = %i\n",Nrad,Nloop,Kr,Kl,SSF);
  printf("At the limit cycle: umin = %.6e, umax = %.6e, umean = %.6e\n",umin,umax,umean);

  return 0;
}  
