/* Ordered Line Integral Method OLIM-MID for finding the quasi-potential in 3D */
/* Computes the quasi-potential with respect to the equilibrium point that must be
a mesh point. Its index is stored in the variable Iindex */
/* This version of olim3D does not do local factoring */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "findQmatrix.h"
#define PI 3.141592653589793
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-12
#define NX 1001
#define NY 1001
#define NZ 1001
#define K 20
#define E3 0.333333333333333
#define e6 0.166666666666666
// Lorenz system parameters
#define sigma  10.0
#define beta  2.666666666666667
#define rho 24.4
#define SSF 2 // subsampling factor: the output file will have 501^3 numbers if SSF = 2


struct myvector {
  double x;  
  double y;
  double z;
};

struct myvector2 {
	double a1;
	double a2;
};	

struct neibcount {
	long na;
	long nc;
};	

struct index3 {
	long i;
	long j;
	long k;
};

struct mymatrix {
  double a11;
  double a12;
  double a21;
  double a22;
};

struct mymatrix3 {
	double a11;
	double a12;
	double a13;
	double a21;
	double a22;
	double a23;
	double a31;
	double a32;
	double a33;
};	

struct mysol {
  double u; // min value
  double a; // minimizer
  char c;   // indicator of success
};  

struct mysol2 {
  double a1; // minimizer component 1
  double a2; // minimizer component 2
  double u;  // min value
  char c;    // indicator of success
};  
	
struct sol_info {
	char type; // solution type:	1 = 1ptupdate, 2 = 2ptupdate, 3 = 3ptupdate, 0 = initialization, 'n' = never reached
	long ind0; // ind0 gives the best one-pt-update
	double U1pt; // the best one-point update value of U 
};
	
struct interpdata {
    double fx0;
    double fx1;
    double fx2;
    double fx3;
    double fx4;
    double fx5;
    double fx6;
    double fx7;
    double fy0;
    double fy1;
    double fy2;
    double fy3;
    double fy4;
    double fy5;
    double fy6;
    double fy7;
    double fz0;
    double fz1;
    double fz2;
    double fz3;
    double fz4;
    double fz5;
    double fz6;
    double fz7;
    
};

struct matrixS3 {
	double a11;
	double a12;
	double a13;
	double a22;
	double a23;
	double a33;
};	
	



int main(void);
struct myvector myfield(struct myvector x); /* B */
void param(void);
void olim(void);
struct mysol triangle_update(long ind,long ind0,long ind1);
double one_pt_update(long ind,long ind0);
void addtree(long ind); /* adds a node to the binary tree
                        of the "considered" points */
void updatetree(long ind); /* updates the binary tree */
void deltree(void); /* deletes the root of the binary tree */
struct myvector a_times_vec(struct myvector v,double a);
struct myvector cross_product(struct myvector a,struct myvector b);
double dot_product(struct myvector a,struct myvector b);
struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b);
struct myvector vec_difference(struct myvector v1,struct myvector v2);
struct myvector vec_sum(struct myvector v1,struct myvector v2);
struct myvector getpoint(long ind);
double length_vec(struct myvector x);
struct index3 getindex3(long ind);
struct mysol hybrid_nonlin_solver(double u0,double u1,
            struct myvector x0,struct myvector x1,struct myvector b0,struct myvector b1,struct myvector x);
double myfun(double s,double du,struct myvector b0,
             struct myvector X01,struct myvector B10,struct myvector xmx0);
             
double geometric_action_line(struct myvector x0, struct myvector x1);
double init(struct myvector x);
double triple_product(struct myvector a,struct myvector b,struct myvector c);
struct myvector vec_lin_comb3(struct myvector x0,struct myvector x1,struct myvector x2,
                              double s,double t);
double lin_comb3(double f0,double f1,double f2,double s,double t);
double geometric_action_line_simp(struct myvector x0, struct myvector x1);
struct mysol2 nonlin_solver2(char KKT,double a1,double a2,double u0,double u1,double u2,
		struct myvector x0,struct myvector x1,struct myvector x2,
		struct myvector b0,struct myvector b1,struct myvector b2,struct myvector x);


double norm2squared(struct myvector x);
long get_neii_index( long d );
long far_neighbors_index_list( long *farlist );
double LinearQpot3D( struct myvector x );


void ShootMAP();
long get_nearest_meshpoint_index(struct myvector p);
struct interpdata get_interpdata(struct myvector p,long ind);
struct myvector myrhs(struct myvector p,struct interpdata gdata,long ind);

/***************************************/

const long nx1 = NX - 1,ny1 = NY - 1,nz1 = NZ - 1,NXY = NX*NY,NXYZ = NX*NY*NZ,nx2 = NX - 2,ny2 = NY - 2,nz2 = NZ - 2;
long KK;
long count = 0; 
double hmax,hx,hy,hz,h,dmaxsquared;
char ms[NX*NY*NZ]; /* 0 = 'Unknown', 1 = 'Considered', 2 = "in Accepted Front", 3 = "Accepted" */
double Uexact[NX*NY*NZ], g[NX*NY*NZ]; /* function to be computed */
long pos[NX*NY*NZ]; /* pos(index of mesh pt) = position in binary tree */
long tree[NX*NY*NZ]; /* tree(position in the tree) = index of mesh pt */
struct sol_info solinfo[NX*NY*NZ];
struct myvector x_ShootMAP;
long Iindex; // the index of the equilibrium point with respect to which the quasi-potential will be computed
long N1ptu,N2call,N3call,N3ptu,N2ptu;


long Nfar; // the number of Considered points in the far neighborhood
long *farlist; // the list of index shifts for defining points in the far neighborhood

// variables for the potential
double XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX; // define the computational domain

// variables for shooting an instanton;
struct myvector *path, x_ipoint;
long Npath;
long Npathmax = 20*max(max(NX,NY),NZ);
struct myvector x_ShootMAP; // the point from which the instanton is shot backward to the attractor

// definitions for exploring the nearest neighborhood of Accepted points
// the number of admissible simplex updates for minimal triangle update
long nsimplex[18] = {8,8,8,8,8,8,2,2,2,2,2,2,2,2,2,2,2,2};
long isimplex[18][8] = { {1,3,4,5,6,9,10,11},
						 {0,2,4,5,6,7,14,17},
						 {1,3,4,5,7,8,12,13},
						 {0,2,4,5,8,9,15,16},
						 {0,1,2,3,11,12,14,15},
						 {0,1,2,3,10,13,14,16},
						 {0,1,-1,-1,-1,-1,-1,-1},	
						 {1,2,-1,-1,-1,-1,-1,-1},	
						 {2,3,-1,-1,-1,-1,-1,-1},	
						 {0,3,-1,-1,-1,-1,-1,-1},	
						 {0,5,-1,-1,-1,-1,-1,-1},	
						 {0,4,-1,-1,-1,-1,-1,-1},	
						 {2,4,-1,-1,-1,-1,-1,-1},	
						 {2,5,-1,-1,-1,-1,-1,-1},	
						 {1,4,-1,-1,-1,-1,-1,-1},	
						 {3,4,-1,-1,-1,-1,-1,-1},	
						 {3,5,-1,-1,-1,-1,-1,-1},	
						 {1,5,-1,-1,-1,-1,-1,-1},	
};
// const long NXp1,NXm1,NXYp1,NXYm1,NXYpNX,NXYmNX,NXYpNXp1,NXYpNXm1,NXYmNXp1,NXYmNXm1;
const long    NXp1 = NX + 1,
    NXm1 = NX - 1,
    NXYp1 = NXY + 1,
    NXYm1 = NXY - 1,
    NXYpNX = NXY + NX,
    NXYmNX = NXY - NX,
    NXYpNXp1 = NXY + NX + 1,
    NXYpNXm1 = NXY + NX - 1,
    NXYmNXp1 = NXY - NX + 1,
    NXYmNXm1 = NXY - NX - 1;


struct mymatrix3 Qmatrix;

/**************************************/
struct myvector myfield(struct myvector x) {
    struct myvector v;
    
    v.x = sigma*(x.y - x.x);
    v.y = x.x*(rho - x.z) - x.y;
    v.z = x.x*x.y - beta*x.z;
    
    return v;
}

/*************************************/
/*************************************/

void param() {
    long ind;
    double r1 = rho - 1, sq, rad; 
    
    sq = sqrt(beta*r1); // sqrt(beta*(rho - 1)
    
    KK = K*K;
    
	Iindex = nx1/2 + NX*ny1/2 + NXY*nz1/2;
	if( rho > 1.0 ) {
		// positive equilibrium point of Lorenz'63
		x_ipoint.x = sq;
		x_ipoint.y = sq;
		x_ipoint.z = r1;
	}
	else {
		// the origin is asymptotically stable equilibrium
		x_ipoint.x = 0.0;
		x_ipoint.y = 0.0;
		x_ipoint.z = 0.0;
	}
		
	rad =  1.5*length_vec(x_ipoint); // for rho = 24.4
	// rad = 12.0; // for rho = 20.0
	// rad = 13.0;// for rho = 15
	
	XMIN = x_ipoint.x - rad; XMAX = x_ipoint.x + rad;
	YMIN = x_ipoint.y - rad; YMAX = x_ipoint.y + rad;
	ZMIN = x_ipoint.z - rad; ZMAX = x_ipoint.z + rad;
	hx = (XMAX - XMIN)/nx1;
	hy = (YMAX - YMIN)/ny1;
	hz = (ZMAX - ZMIN)/nz1;
	x_ipoint = getpoint(Iindex);
	printf("x_ipoint: Iindex = %li, %.4e, %.4e, %.4e\n",Iindex,x_ipoint.x,x_ipoint.y,x_ipoint.z);


    Qmatrix = findQmatrix(sigma,beta,rho); // set up struct matrixS3 Qmatrix
    printf("in param()\n");
    h = sqrt(hx*hx + hy*hy + hz*hz);
    hmax = max(max(hx,hy),hz);
    dmaxsquared = max(KK*hmax*hmax,hmax*hmax*3.0) + 2.0e-16;
    for( ind=0; ind<NXYZ; ind++ ) {
        ms[ind] = 0;
        g[ind] = INFTY;
        solinfo[ind].type = 'n';
        solinfo[ind].U1pt = INFTY;
    }
    // compute the exact solution if it is available
    
}

/************************************/

void ipoint() {
  long ind,ind0,n,m;
    const long neii[26] = {1, NX, -1, -NX, NXY, -NXY, // 0 -- 5: l1 distance = 1 neighbors
    			NXp1, NXm1, -NXp1, -NXm1, // 6 -- 9: l1 distance = 2, xy-plane
    			-NXYm1, NXYp1, NXYm1, -NXYp1, // 10 -- 13: l1 distance = 2, xz-plane
    			NXYpNX, NXYmNX, -NXYpNX, -NXYmNX, // 14 -- 17: l1 distance = 2, yz-plane
    			NXYpNXp1, NXYpNXm1,NXYmNXp1, NXYmNXm1, // 18 -- 21: l1 distance = 3, top four corners
    			-NXYmNXm1, -NXYmNXp1, -NXYpNXm1, -NXYpNXp1 // 22 -- 25: l1 distance = 3, bottom four corners    						 
    		}; 
  struct myvector x;
  
  	printf("in ipoint\n");
  
	g[Iindex] = 0.0;
	solinfo[Iindex].type = 0;
	ms[Iindex] = 3;
	for( n = 0; n < 26; n++ ) {
		ind = Iindex + neii[n];
		x = getpoint(ind);
		g[ind] = LinearQpot3D(vec_difference(x,x_ipoint));
		ms[ind] = 2; 
		solinfo[ind].type = 0;
	}
	for( n = 0; n < 26; n++ ) {
		ind0 = Iindex + neii[n];
		for(  m = 0; m < 26; m++ ) {
			ind = ind0 + neii[m];
			if( ms[ind] == 0 ) {
				ms[ind] = 1;
				x = getpoint(ind);
				g[ind] = LinearQpot3D(vec_difference(x,x_ipoint));
				solinfo[ind].type = 0;
				addtree(ind);
			}
		}
	}
				
}

/********************************************/

double LinearQpot3D( struct myvector x ) {

   return Qmatrix.a11*x.x*x.x + Qmatrix.a22*x.y*x.y + Qmatrix.a33*x.z*x.z 
    	+ 2.0*(Qmatrix.a12*x.x*x.y + Qmatrix.a13*x.x*x.z + Qmatrix.a23*x.y*x.z);
}


/**********************************************/
/*** ordered line longegral method ***/

void olim(void) {
    long k,m,n,ind,inew,ind0,ind1,ind2,neii_index,ic;
    double gtemp,gold;
    long Naf,Nc,AFneib[26],NCneib[26]; 
    long NAC = 0; // the number of Accepted points
//     long Nfar; // the number of Considered points in the far neighborhood
//     long *farlist; // the list of index shifts for defining points in the far neighborhood
    struct mysol sol;
    struct mysol2 sol2;
    struct myvector vec,b0,b1,b2,v0,v1,v2,vnew;
    struct index3 pnew;
    char chh;
    /* neighbor's indices  with L_inf distance <= 1 and L1 distance <= 2*/
    const long neii[26] = {1, NX, -1, -NX, NXY, -NXY, // 0 -- 5: l1 distance = 1 neighbors
    						 NXp1, NXm1, -NXp1, -NXm1, // 6 -- 9: l1 distance = 2, xy-plane
    						 -NXYm1, NXYp1, NXYm1, -NXYp1, // 10 -- 13: l1 distance = 2, xz-plane
    						 NXYpNX, NXYmNX, -NXYpNX, -NXYmNX, // 14 -- 17: l1 distance = 2, yz-plane
    						 NXYpNXp1, NXYpNXm1,NXYmNXp1, NXYmNXm1, // 18 -- 21: l1 distance = 3, top four corners
    						 -NXYmNXm1, -NXYmNXp1, -NXYpNXm1, -NXYpNXp1 // 22 -- 25: l1 distance = 3, bottom four corners    						 
    						 }; 
    
    printf("in olim()\n");
    
    
    while( count > 0 ) {
        inew = tree[1];
        vnew = getpoint(inew);
        pnew = getindex3(inew);
        
         
    	ms[inew] = 2;
        deltree();
        NAC++;
        
        if( pnew.i<=1 || pnew.i>=nx2 || pnew.j<=1 || pnew.j>= ny2 || pnew.k <= 1 || pnew.k >= nz2 || g[inew] >= INFTY-1) {
            printf("The boundary is reached:\n");
            printf("%li accepted points,(%li,%li,%li) is accepted, g=%.4e\n",NAC,pnew.i,pnew.j,pnew.k,g[inew]);
            break; /* quit if we reach the boundary of the computational domain */
        }
        
        
        /* Inspect the neighbors of the new Accepted point */
        // Find nearest neighbors to be shifted to Accepted
        // List nearest neighbors in N18 neighborhood that are AF to use for updates
        // List Unknown neighbors that will be new considered
        Naf = 0;
        Nc = 0;
        for( k = 0; k < 18; k++ ) {
            ind1 = inew + neii[k]; // neighbor of the new Accepted point
            // 		  printf("ind1 = %li, ms = %li\n",ind1,ms[ind1]);
            // update AcceptedFront
            if( ms[ind1] == 2 ) {
                m = 0;
                for( n = 0; n < 26; n++ ) {
                    ind0 = ind1 + neii[n];
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
                //vec = getpoint(ind);
                NCneib[Nc] = ind1;
                Nc++;
            }
        }
        for( k = 18; k < 26; k++ ) {
            ind1 = inew + neii[k]; // neighbor of the new Accepted point
            // 		  printf("ind1 = %li, ms = %li\n",ind1,ms[ind1]);
            // update AcceptedFront
            if( ms[ind1] == 2 ) {
                m = 0;
                for( n = 0; n < 26; n++ ) {
                    ind0 = ind1 + neii[n];
                    if( ms[ind0] < 2 ) m++;
                }
                if( m == 0 ) { /* ind1 has no considered neighbors */
                    ms[ind1] = 3;
                }
            }
            else if( ms[ind1] == 0 ) { // the neighbor ind1 will be a new Considered point
                //vec = getpoint(ind);
                NCneib[Nc] = ind1;
                Nc++;
            }
        }
//          printf("Naf = %li, Nc = %li\n",Naf,Nc);
        
        //----- UPDATE EXISTING CONSIDERED POINTS IN THE FAR NEIGHBORHOOD OF THE NEW AF POINT
        /* update considered points */
        
 	for( k = 0; k < Nfar; k++ )  {
		ind = inew + farlist[k];
		if( ind < NXYZ && ind >= 0 ) {
			if( ms[ind] == 1 ) {
				if( solinfo[ind].type > 0 ) {
					vec = getpoint(ind); 
					if( norm2squared(vec_difference(vec,vnew)) <= dmaxsquared ) {
// 						Vflag = (norm2squared(vec_difference(vec,x_ipoint)) < Rfac2) ? 'y' : 'n'; // for local factoring
						gold = g[ind];
						gtemp = one_pt_update(ind,inew);
						if( gtemp < g[ind] ) {
							g[ind] = gtemp;			
							solinfo[ind].type = 1;
							solinfo[ind].ind0 = inew;
							solinfo[ind].U1pt = gtemp;
						}
						else if( gtemp < solinfo[ind].U1pt ) { // the smallest 1pt update value
							solinfo[ind].ind0 = inew;
							solinfo[ind].U1pt = gtemp;
						}	
		

				// 		p = getindex3(ind);
						b0 = myfield(vec_lin_comb(vnew,vec,0.5,0.5));
		
						//  CASE 1: if the new AF point does not give the minimum of one-point update
						if( inew != solinfo[ind].ind0 ) { // there is only one triangle update to exercise
							ind1 = solinfo[ind].ind0;
				// 			p1 = getindex3(ind1);
							neii_index = get_neii_index(inew - ind1);
							if( neii_index >= 0 ) { // ind0,ind1,ind form a triangle
								v1 = getpoint(ind1);
								b1 = myfield(vec_lin_comb(v1,vec,0.5,0.5));
 								sol = hybrid_nonlin_solver(g[ind1],g[inew],v1,vnew,b1,b0,vec); // this order is important
								if( sol.c == 'y' ) {  // triangle update is successful
									if( sol.u < g[ind] ) {
										g[ind] = sol.u;
										solinfo[ind].type = 2; // solinfo[ind].ind0 remains ind1
									}
									// Do 3ptu whenever 2ptu succeeds	
									for( m = 0; m < nsimplex[neii_index]; m++ ) {
										ind2 = ind1 + neii[isimplex[neii_index][m]];
										if( ms[ind2] == 2 ) {
											v2 = getpoint(ind2);
											b2 = myfield(vec_lin_comb(v2,vec,0.5,0.5));
 											sol2 = nonlin_solver2('y',sol.a,0.0,g[ind1],g[inew],g[ind2],v1,vnew,v2,b1,b0,b2,vec);
											if( sol2.c == 'y' && sol2.u < g[ind] ) {
												g[ind] = sol2.u;
												solinfo[ind].type = 3;
											}
										}
									}
								
								} // if( sol.c == 'y' )
							} // if( neii_index >= 0 ) 
						}	// if( inew != solinfo[ind].ind0 )  -- if the new AF point does not give min of one-point update				
							
						// CASE 2: if the new AF point gives minimum of one-point update			
						else {  // if( inew != solinfo[ind].ind0 ) // try all possible triangle updates
	// 						sflag = 'n';
							for( m = 0; m < Naf; m++ ) {
								ind1 = AFneib[m];
								v1 = getpoint(ind1);
								b1 = myfield(vec_lin_comb(v1,vec,0.5,0.5));
 								sol = hybrid_nonlin_solver(g[inew],g[ind1],vnew,v1,b0,b1,vec); // this order is important
								if( sol.c == 'y' ) {  // triangle update is successful
									if( sol.u < g[ind] ) {
										g[ind] = sol.u;
										solinfo[ind].type = 2; // solinfo[ind].ind0 remains ind1
									}

									neii_index = get_neii_index(ind1 - inew);
									// PRECAUTION
									if( neii_index < 0 ) {
										printf("Error in CASE 2: inew = %li, ind1 = %li, ms[inew] = %i, ms[ind1] = %i\n",inew,ind1,ms[inew],ms[ind1]);
										exit(1);
									}	
									for( n = 0; n < nsimplex[neii_index]; n++ ) {
										ind2 = inew + neii[isimplex[neii_index][n]];  // unew used to be ind0
									
										if( ms[ind2] == 2 ) {
											// To avoid repetitions of simplex updates
											chh = 'n'; ic = m + 1;
											while( chh == 'n' && ic < Naf ) {
												if( ind2 == AFneib[ic] ) chh = 'y';
												ic++;
											}
											if( chh == 'y' ) {	   									
												v2 = getpoint(ind2);
												b2 = myfield(vec_lin_comb(v2,vec,0.5,0.5));
												sol2 = nonlin_solver2('y',sol.a,0.0,g[inew],g[ind1],g[ind2],vnew,v1,v2,b0,b1,b2,vec);
												if( sol2.c == 'y' && sol2.u < g[ind] ) {
													g[ind] = sol2.u;
													solinfo[ind].type = 3;
												}
											}
										}
									}	
								}
							} // for( m = 0; m < Naf; m++ ) 
						} // end else
					if( g[ind] < gold ) updatetree(ind);
					} // end if( norm2squared(vec_difference(vec,vnew)) <= dmaxsquared )
				} // end if( solinfo[ind].type > 0 )
		    } // end if( ms[ind] == 1 )	
		} // end if( ind < NXYZ && ind >= 0 )
	} // end for( k = 0; k < Nfar; k++ ) 		

// 	printf("count = %li\n",count);
				
	// UPDATE NEW CONSIDERED POINTS			
	for( k = 0; k < Nc;	k++ ) {
		ind = NCneib[k];
		vec = getpoint(ind);
		ms[ind] = 1;
		vec = getpoint(ind);
		
		// find minimal one-point update
		for( n = 0; n < Nfar; n++ ) {
			ind0 = ind + farlist[n];
			if( ind0 < NXYZ && ind0 >= 0 )  {
				if( ms[ind0] == 2 ) {
					v0 = getpoint(ind0); 
					if( norm2squared(vec_difference(vec,v0)) <= dmaxsquared ) {
						gtemp = one_pt_update(ind,ind0);
						if( gtemp < g[ind] ) {
							g[ind] = gtemp;			
							solinfo[ind].type = 1;
							solinfo[ind].ind0 = ind0;
							solinfo[ind].U1pt = gtemp;
						}
					}
				}
			}
		} // for( n = 0; n < Nfar; n++ ) 
		// do triangle update with the minimizer of one-point update
		ind0 = solinfo[ind].ind0;
		v0 = getpoint(ind0);
		b0 = myfield(vec_lin_comb(v0,vec,0.5,0.5));
		// find minimal triangle update
		for( m = 0; m < 18; m++ ) {
			ind1 = ind0 + neii[m];
			if( ms[ind1] == 2 ) {
				v1 = getpoint(ind1);
				b1 = myfield(vec_lin_comb(v1,vec,0.5,0.5));
 				sol = hybrid_nonlin_solver(g[ind0],g[ind1],v0,v1,b0,b1,vec); // this order is important
				if( sol.c == 'y' ) {  // triangle update is successful
					if( sol.u < g[ind] ) {
						g[ind] = sol.u;
						solinfo[ind].type = 2; // solinfo[ind].ind0 remains ind1
					}
					// try simplex update whenever triangle update succeeds
					neii_index = get_neii_index(ind1 - ind0);
					for( n = 0; n < nsimplex[neii_index]; n++ ) {
						ind2 = ind0 + neii[isimplex[neii_index][n]];						
						if( ms[ind2] == 2 ) {
						
							   // To avoid repetitions of simplex updates
						       chh = 'n'; ic = m + 1;
						       while( chh == 'n' && ic < Naf ) {
						       	   if( ind2 == AFneib[ic] ) chh = 'y';
						       	   ic++;
						       }	   

						
							if( chh == 'y' ) {
								v2 = getpoint(ind2);
								b2 = myfield(vec_lin_comb(v2,vec,0.5,0.5));
								sol2 = nonlin_solver2('y',sol.a,0.0,g[ind0],g[ind1],g[ind2],v0,v1,v2,b0,b1,b2,vec);
								if( sol2.c == 'y' && sol2.u < g[ind] ) {
									g[ind] = sol2.u;
									solinfo[ind].type = 3;
								}
							}
						}
					}
					
					
					
				}
			}
		} // for( m = 0; m < 18; m++ )	
		addtree(ind);
// 		bb = myfield(vec);
// 		pp = getindex3(ind);
	} // end for( k = 0; k < Nc; k++ )	
	
// 	printf("count = %li\n",count);
  } // end while
}			
				

/**********************************************/
//     const long neii[18] = {1, NX, -1, -NX, NXY, -NXY, // 0 -- 5: l1 distance = 1 neighbors
//     						 1 + NX, -1 + NX, -1 - NX, 1 - NX, // 6 -- 9: l1 distance = 2, xy-plane
//     						 1 - NXY, 1 + NXY, -1 + NXY, -1 - NXY, // 10 -- 13: l1 distance = 2, xz-plane
//     						 NX + NXY, -NX + NXY, -NX - NXY, NX - NXY}; // 14 -- 17: l1 distance = 2, yz-plane


long get_neii_index( long d ) {	

	switch(d) {
		case 1:
			return 0;
			break;
		case NX:
			return 1;
			break;
		case -1:
			return 2;
			break;
		case -NX:
			return 3;
			break;
		case NXY:
			return 4;
			break;
		case -NXY:
			return 5;
			break;
		case NXp1:
			return 6;
			break;
		case NXm1:
			return 7;
			break;
		case -NXp1:
			return 8;
			break;
		case -NXm1:
			return 9;	
			break;
		case -NXYm1:
			return 10;
			break;
		case NXYp1:
			return 11;
			break;
		case NXYm1:
			return 12;
			break;
		case -NXYp1:
			return 13;
			break;
		case NXYpNX:
			return 14;
			break;
		case NXYmNX:
			return 15;
			break;													
		case -NXYpNX:
			return 16;
			break;
		case -NXYmNX:
			return 17;
			break;
		default:
			return -1;
			break;
	}	
					
}








/*********************************************/

double one_pt_update(long ind,long ind0) {
    struct myvector x,x0;
    double gtemp;
    
    x = getpoint(ind);
    x0 = getpoint(ind0);
    gtemp = g[ind0] + geometric_action_line(x0,x);
    
    N1ptu++;
    //  fprintf(fup,"%li\t%.6e\t%li\t%li\t%li\t%li\t%.4e\n",ind,gtemp,1,ind0,-1,-1,length(l.x,l.y));
    
    return gtemp;
}
/*-------------*/

double geometric_action_line(struct myvector x0, struct myvector x1) {
    struct myvector l,b;
    
    l = vec_difference(x1,x0);
    b = myfield(vec_lin_comb(x0,x1,0.5,0.5));
    
    return length_vec(b)*length_vec(l) - dot_product(b,l);
}

double geometric_action_line_simp(struct myvector x0, struct myvector x1) {
    struct myvector l,b0,bm,b1;
    
    l = vec_difference(x1,x0);
    b0 = myfield(x0);
    bm = myfield(vec_lin_comb(x0,x1,0.5,0.5));
    b1 = myfield(x1);
    
    return ((length_vec(b0)+4.0*length_vec(bm)+length_vec(b1))*length_vec(l) 
    		- dot_product(vec_lin_comb(vec_lin_comb(b0,b1,1.0,1.0),bm,1.0,4.0),l))/6.0;
}


/*-------------*/

struct myvector getpoint(long ind) {
    struct myvector l;
    
    l.x = hx*((ind%NXY)%NX) + XMIN;
    l.y = hy*((ind%NXY)/NX) + YMIN;
    l.z = hz*(ind/NXY) + ZMIN;
    return l;
}



/***** N o n l i n e a r   1D    s o l v e r *****/


struct mysol hybrid_nonlin_solver(double u0,double u1,
            struct myvector x0,struct myvector x1,struct myvector b0,struct myvector b1,struct myvector x) {
    double a = 0.0, b = 1.0, c, fa, fb, fc, d, fd, dd, df, dm, ds, t,du = u1 - u0;
    struct myvector X01,B10,xs,xmx0;
    struct mysol sol;
    double NONLIN_TOL = 1.0e-6;
    long iter = 0, itermax = 100;
    
    X01 = vec_difference(x0,x1);
    B10 = vec_difference(b1,b0);
    xmx0 = vec_difference(x,x0);
    
    c = a;
    fa = myfun(a,du,b0,X01,B10,xmx0);
    fb = myfun(b,du,b0,X01,B10,xmx0);
    fc = fa;
    
    N2call++;
    
    if( (fa > 0 && fb > 0 ) || (fa < 0 && fb < 0) ) {
        //	 root is not bracketed
        sol.c = 'n';
        sol.u = INFTY;
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
        fd = myfun(d,du,b0,X01,B10,xmx0);
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
    sol.a = b;
    xs = vec_difference(x0,a_times_vec(X01,b));
    sol.u = u0 + b*du + geometric_action_line(xs,x);
    
    N2ptu++;
    
    return sol;
}

/*--------------*/

double myfun(double s,double du,struct myvector b0,
             struct myvector X01,struct myvector B10,struct myvector xmx0) {
    double ls,lbs,bat;
    struct myvector xmxs,bs;
    
    xmxs = vec_sum(xmx0,a_times_vec(X01,s));
    bs = vec_sum(b0,a_times_vec(B10,s));
    ls = length_vec(xmxs);
    lbs = length_vec(bs);
    bat = lbs/ls;
    
    if( bat > TOL ) {
        return du + dot_product(xmxs,X01)*bat + dot_product(bs,B10)/bat
        - dot_product(X01,bs) - dot_product(B10,xmxs);
    }
    else {
        return du - dot_product(B10,xmxs);
    }
    
}


/**********************************************/

struct mysol2 nonlin_solver2(char KKT,double a1,double a2,double u0,double u1,double u2,
		struct myvector x0,struct myvector x1,struct myvector x2,
		struct myvector b0,struct myvector b1,struct myvector b2,struct myvector x) {
		
	struct mysol2 sol;
	double bat,ibat,Hdet,ls,lbs;
	int iter = 0, iter_max = 10;
	struct myvector2 Udiff,Xv,Bw,Bxmxs,grad,p;
	struct myvector X01,X02,B10,B20,v,w,xmxs,xs,bs;
	double H11,H12,H22; // entries of the Hessian matrix
	double BX11,BX12,BX22,XX11,XX12,XX22,BB11,BB12,BB22;
	double safefac = 1.0 - 1e-12,fac1,fac2,fac3,slen;
	double Ntol = 1e-12;
	
	N3call++;
	
	B10 = vec_difference(b1,b0);
	B20 = vec_difference(b2,b0);
	X01 = vec_difference(x0,x1);
	X02 = vec_difference(x0,x2);
	
	Udiff.a1 = u1 - u0;
	Udiff.a2 = u2 - u0;
	
	BX11 = 2.0*dot_product(B10,X01);
	BX12 = dot_product(B20,X01) + dot_product(B10,X02);
	BX22 = 2.0*dot_product(B20,X02);
	
	XX11 = dot_product(X01,X01);
	XX12 = dot_product(X01,X02);
	XX22 = dot_product(X01,X02);
	
	BB11 = dot_product(B10,B10);
	BB12 = dot_product(B10,B20);
	BB22 = dot_product(B20,B20);
	
	
// 	if( IND == 521 && IND0 == 653 && IND1 == 654 && IND2 == 533 ) ch = 'y';

	xs = vec_difference(x0,vec_lin_comb(X01,X02,a1,a2)); // x0 - a1*X01 - a2*X02
	xmxs = vec_difference(x,xs);
	ls = length_vec(xmxs);
	v = a_times_vec(xmxs,1.0/ls); // unit vector in the direction of xmxs
	Xv.a1 = dot_product(X01,v);
	Xv.a2 = dot_product(X02,v);
	bs = vec_sum(b0,vec_lin_comb(B10,B20,a1,a2)); // b0 + a1*(b1-b0) + a2*(b2-b0)
	lbs = length_vec(bs);
	Bxmxs.a1 = dot_product(B10,xmxs);
	Bxmxs.a2 = dot_product(B20,xmxs);
	if( lbs < Ntol )  {
		grad.a1 = Udiff.a1 - Bxmxs.a1;
		grad.a2 = Udiff.a2 - Bxmxs.a2;
	}
	else {
		w = a_times_vec(bs,1.0/lbs);
		Bw.a1 = dot_product(B10,w);
		Bw.a2 = dot_product(B20,w);
		grad.a1 = Udiff.a1 + ls*Bw.a1 + lbs*Xv.a1 - Bxmxs.a1 - dot_product(X01,bs);
		grad.a2 = Udiff.a2 + ls*Bw.a2 + lbs*Xv.a2 - Bxmxs.a2 - dot_product(X02,bs);
	}		
	
// 	if( ch == 'y' ) {
// 		printf("iter = %i: a1 = %.4e, a2 = %.4e,grad = [%.4e,%.4e]\n",iter,a1,a2,grad.a1,grad.a2);
// 	}	
	// perform the KKT test
	if( KKT == 'y' ) { // do KKT test: check if the initial point is a local solution
		if( grad.a2 >= 0.0 ) { // then the minimum on the boundary is a local solution
			sol.c = 'n';
			return sol;
		}	
	}
	
	N3ptu++;
	
	while( iter < iter_max && max(fabs(grad.a1),fabs(grad.a2)) > Ntol ) {
		bat = lbs/ls; 
		H11 = bat*(XX11 - Xv.a1*Xv.a1) - BX11;
		H12 = bat*(XX12 - Xv.a1*Xv.a2) - BX12;
		H22 = bat*(XX22 - Xv.a2*Xv.a2) - BX22;
		if( lbs >= Ntol ) {
			ibat = 1.0/bat;
			H11 += 2.0*Xv.a1*Bw.a1 + (BB11 - Bw.a1*Bw.a1)*ibat;
			H12 += (Xv.a1*Bw.a2 + Xv.a2*Bw.a1) + (BB12 - Bw.a1*Bw.a2)*ibat;
			H22 += 2.0*Xv.a2*Bw.a2 + (BB22 - Bw.a2*Bw.a2)*ibat;
		}
		// check if H is positive definite
		Hdet = H11*H22 - H12*H12;
		if( Hdet > 1e-12 && H11 + H22 > 0.0 ) { // H is pos def
			// find Newton direction 
			p.a1 = -(grad.a1*H22 - grad.a2*H12)/Hdet;
			p.a2 = -(H11*grad.a2 - H12*grad.a1)/Hdet;
		}	
		else {
			// take steepest descend direction
			p.a1 = -grad.a1;
			p.a2 = -grad.a2;
		}
		// find step length and make a step
		fac1 = (p.a1 < 0.0) ? -a1/p.a1 : INFTY;
		fac2 = (p.a2 < 0.0) ? -a2/p.a2 : INFTY;
		fac3 = (p.a1 + p.a2 > 0.0) ? (1.0 - a1 - a2)/(p.a1 + p.a2) : INFTY;
		
		slen = safefac*min(1,min(min(fac1,fac2),fac3)); // step length
// 		if( ch == 'y' ) {
// 			printf("fac1 = %.4e, fac2 = %.4e, fac3 = %.4e\n",fac1,fac2,fac3);
// 			printf("iter = %i: a1 = %.4e, a2 = %.4e,grad = [%.4e,%.4e], slen = %.4e,p = [%.4e,%.4e]\n",
// 				iter,a1,a2,grad.a1,grad.a2,slen,p.a1,p.a2);
// 		}	
		a1 += slen*p.a1; //if( a1 < 0.0 ) a1 = 0.0;
		a2 += slen*p.a2; //if( a2 < 0.0 ) a2 = 0.0;
		//if( a1 + a2 > 1.0 ) { a1 /= (a1 + a2); a2 /= (a1 + a2); }	
		
		if( a1 < 0.0 || a2 < 0.0 || a1 + a2 > 1.0 ) {
			printf("NONLIN_SOLVER2: a1 = %.4e, a2 = %.4e\n",a1,a2);
// 			printf("IND = %i, IND0 = %i, IND1 = %i, IND2 = %i\n",IND,IND0,IND1,IND2);
			printf("fac1 = %.4e, fac2 = %.4e, fac3 = %.4e\n",fac1,fac2,fac3);
			printf("iter = %i: a1 = %.4e, a2 = %.4e,grad = [%.4e,%.4e], slen = %.4e,p = [%.4e,%.4e]\n",
				iter,a1,a2,grad.a1,grad.a2,slen,p.a1,p.a2);
			exit(1);
		}	
			
		xs = vec_difference(x0,vec_lin_comb(X01,X02,a1,a2)); // x0 - a1*X01 - a2*X02
		xmxs = vec_difference(x,xs);
		ls = length_vec(xmxs);
		v = a_times_vec(xmxs,1.0/ls); // unit vector in the direction of xmxs
		Xv.a1 = dot_product(X01,v);
		Xv.a2 = dot_product(X02,v);
		bs = vec_sum(b0,vec_lin_comb(B10,B20,a1,a2)); // b0 + a1*(b1-b0) + a2*(b2-b0)
		lbs = length_vec(bs);
		Bxmxs.a1 = dot_product(B10,xmxs);
		Bxmxs.a2 = dot_product(B20,xmxs);
		if( lbs < Ntol )  {
			grad.a1 = Udiff.a1 - Bxmxs.a1;
			grad.a2 = Udiff.a2 - Bxmxs.a2;
		}
		else {
			w = a_times_vec(bs,1.0/lbs);
			Bw.a1 = dot_product(B10,w);
			Bw.a2 = dot_product(B20,w);
			grad.a1 = Udiff.a1 + ls*Bw.a1 + lbs*Xv.a1 - Bxmxs.a1 - dot_product(X01,bs);
			grad.a2 = Udiff.a2 + ls*Bw.a2 + lbs*Xv.a2 - Bxmxs.a2 - dot_product(X02,bs);
		}		
		iter++;
	}
	if( max(fabs(grad.a1),fabs(grad.a2)) < Ntol  ) {
		sol.c = 'y';
		sol.a1 = a1;
		sol.a2 = a2;
		sol.u = u0 + a1*Udiff.a1 + a2*Udiff.a2 + lbs*ls - dot_product(xmxs,bs);
	}	
	else {
		sol.c = 'n';
	}	
	return sol;	
}		
		
/**********************************************/


double length3(double x,double y,double z) {
    return sqrt(x*x+y*y+z*z);
}

//--------------------

double length_vec(struct myvector x) {
    return sqrt(norm2squared(x));
}

//--------------------

double norm2squared(struct myvector x) {
    return x.x*x.x + x.y*x.y + x.z*x.z;
}

	
/***********/
double dot_product(struct myvector a,struct myvector b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

double triple_product(struct myvector a,struct myvector b,struct myvector c) {
    
    return a.x*(b.y*c.z - c.y*b.z) - b.x*(a.y*c.z - c.y*a.z) + c.x*(a.y*b.z - b.y*a.z);
}

/*****************/

struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b) {
    struct myvector v;
    
    v.x = a*v1.x + b*v2.x;
    v.y = a*v1.y + b*v2.y;
    v.z = a*v1.z + b*v2.z;
    return v;
}

struct myvector vec_lin_comb3(struct myvector x0,struct myvector x1,struct myvector x2,
                              double s,double t) {
    
    return vec_lin_comb(x0,vec_lin_comb(x1,x2,t,1.0 - t),s,1.0 - s);
}

double lin_comb3(double f0,double f1,double f2,double s,double t) {
    
    return f0*s + (f1*t + f2*(1.0-t))*(1 - s);
}

/*****************/

struct myvector cross_product(struct myvector a,struct myvector b) {
    struct myvector v;
    
    v.x = a.y*b.z - a.z*b.y;
    v.y = a.z*b.x - a.x*b.z;
    v.z = a.x*b.y - a.y*b.x;
    
    return v;
}
/*****************/

struct myvector vec_difference(struct myvector v1,struct myvector v2) {
    struct myvector v;
    
    v.x = v1.x - v2.x;
    v.y = v1.y - v2.y;
    v.z = v1.z - v2.z;
    return v;
}
/*****************/

struct myvector vec_sum(struct myvector v1,struct myvector v2) {
    struct myvector v;
    
    v.x = v1.x + v2.x;
    v.y = v1.y + v2.y;
    v.z = v1.z + v2.z;
    return v;
}
/*****************/
struct myvector a_times_vec(struct myvector v,double a) {
    struct myvector av;
    
    av.x = a*v.x;
    av.y = a*v.y;
    av.z = a*v.z;
    
    return av;
}



/****************************************/

struct index3 getindex3(long ind) {
	struct index3 m;
	long n = ind%NXY;
	
	m.i = n%NX;
	m.j = n/NX;
	m.k = ind/NXY;
	
	return m;
}
	
//*********************************************
long far_neighbors_index_list( long *farlist ) {
	long Nfar = 0;
	long i,j,k,jmax,kmax;
	
	for( i = -K; i <= K; i++ ) {
		jmax = ceil(sqrt(KK - i*i));
		for( j = -jmax; j <= jmax; j++ ) {
			kmax = ceil(sqrt(KK - min(i*i + j*j,KK)));
			for( k = -kmax; k <= kmax; k++ ) {
				if( i != 0 || j != 0 || k != 0 ) {
					farlist[Nfar] = i + j*NX + k*NXY;
					Nfar++;
				}
			}
		}
	}
	return Nfar;				
				
}

/****************************************/
/**************************************************************/
/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(long ind) {
    long loc, ptemp;
    long indp, indc;
    char ch;
    
    count++;
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

void updatetree(long ind) {
    long loc, lcc;
    double g0;
    
    g0=g[ind];
    loc=pos[ind];
    while( loc > 1 && g0 < g[tree[loc/2]] ) {
        tree[loc]=tree[loc/2];
        pos[tree[loc]]=loc;
        loc=loc/2;
        tree[loc]=ind;
        pos[tree[loc]]=loc;
    }
    lcc=count;
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
    long loc, ptemp, ind, lcc, ic, ic1, ic2, mind;
    char chd, ch='n';;
    
    mind=tree[1];
    pos[tree[1]]=0;
    tree[1]=tree[count];
    pos[tree[1]]=1;
    count--;
    loc=1;
    ind=tree[1];
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
        if( (g[ind]) > (g[ic]) ) {chd='l'; if(ch=='y') printf("left\n");}
        else chd='n';
    }
    else chd='n';
    while( chd != 'n' ) {
        ptemp=pos[ind];
        pos[ind]=pos[ic];
        tree[loc]=ic;
        pos[ic]=ptemp;
        tree[lcc]=ind;
        loc=lcc;
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
            if(ch=='y') printf("child: loc(%li)=%li, t1=%.12e\n",ic1,lcc,g[ic1]);
            if( (g[ind]) > (g[ic]) ) { chd='l';if(ch=='y') printf("left\n");}
            else chd='n';
        }
        else chd='n';
    } /* end while( chd != 'n' ) */
}

/********** shoot MAP **************/
void ShootMAP() {
    // shoots a MAP for Maier-Stein
    long j,indprev = -1,ind;
    struct myvector y,k1,k2,k3,k4;
    struct interpdata gdata;
    double dt = 0.1*h,dt2;
    
    path = (struct myvector *)malloc(Npathmax*sizeof(struct myvector));
    
    dt2 = 0.5*dt;
    path[0] = x_ShootMAP;
    j = 0;
    while( length_vec(vec_difference(x_ipoint,path[j])) > h ) {
        y = path[j];
        ind = get_nearest_meshpoint_index(y);
        if( ind != indprev ) gdata = get_interpdata(y,ind);
        k1 = myrhs(y,gdata,ind);
        indprev = ind;
        
        y = vec_sum(path[j],a_times_vec(k1,dt2));
        ind = get_nearest_meshpoint_index(y);
        if( ind != indprev ) gdata = get_interpdata(y,ind);
        k2 = myrhs(y,gdata,ind);
        indprev = ind;
        
        y = vec_sum(path[j],a_times_vec(k2,dt2));
        ind = get_nearest_meshpoint_index(y);
        if( ind != indprev ) gdata = get_interpdata(y,ind);
        k3 = myrhs(y,gdata,ind);
        indprev = ind;
        
        y = vec_sum(path[j],a_times_vec(k3,dt));
        ind = get_nearest_meshpoint_index(y);
        if( ind != indprev ) gdata = get_interpdata(y,ind);
        k4 = myrhs(y,gdata,ind);
        indprev = ind;
        
        path[j+1] = vec_sum(path[j],a_times_vec(vec_sum(vec_sum(k1,a_times_vec(vec_sum(k2,k3),2.0)),k4),h*e6));
        j++;
        printf("j = %li, path: %.4e, %.4e, %.4e, distance: %.4e\n",
        		j,path[j].x,path[j].y,path[j].z,length_vec(vec_difference(x_ipoint,path[j])));
        if( j >= Npathmax ) {
            printf("Failed to compute MAP: the maximal allowed path length is reached: j = %li\n",j);
            return;
        }
    }
    path[j] = x_ipoint;
    j++;
    Npath = j;
    printf("Npath = %li\n",Npath);
}


/*************************************************/
long get_nearest_meshpoint_index(struct myvector p) {
    long ix,iy,iz;
    
    ix = floor((p.x - XMIN)/hx);
    iy = floor((p.y - YMIN)/hy);
    iz = floor((p.z - ZMIN)/hz);
    return ix + iy*NX + iz*NXY;
}

/*************************************************/

struct interpdata get_interpdata(struct myvector p,long ind) {
    struct interpdata gdata;
    
    gdata.fx0 = 0.5*(g[ind + 1] - g[ind - 1])/hx; // dg/dx at ind
    gdata.fx1 = 0.5*(g[ind + 2] - g[ind])/hx; // dg/dx at ind + 1
    gdata.fx2 = 0.5*(g[ind + 1 + NX] - g[ind - 1 + NX])/hx; // dg/dx ind + NX
    gdata.fx3 = 0.5*(g[ind + 2 + NX] - g[ind + NX])/hx; // dg/dx ind + NX + 1
    gdata.fx4 = 0.5*(g[ind + 1 + NXY] - g[ind - 1 + NXY])/hx;
    gdata.fx5 = 0.5*(g[ind + 2 + NXY] - g[ind + NXY])/hx;
    gdata.fx6 = 0.5*(g[ind + 1 + NXY + NX] - g[ind - 1 + NXY + NX])/hx;
    gdata.fx7 = 0.5*(g[ind + 2 + NXY + NX] - g[ind + NXY + NX])/hx;
    
    gdata.fy0 = 0.5*(g[ind + NX] - g[ind - NX])/hy; // dg/dy at ind
    gdata.fy1 = 0.5*(g[ind + NX + 1] - g[ind - NX + 1])/hy; //dg/dy at ind + 1
    gdata.fy2 = 0.5*(g[ind + 2*NX] - g[ind])/hy; // dg/dy at ind + NX
    gdata.fy3 = 0.5*(g[ind + 2*NX + 1] - g[ind + 1])/hy; // dg/dy at ind + NX + 1
    gdata.fy4 = 0.5*(g[ind + NXY + NX] - g[ind + NXY -NX])/hy;
    gdata.fy5 = 0.5*(g[ind + NXY + NX + 1] - g[ind + NXY -NX + 1])/hy;
    gdata.fy6 = 0.5*(g[ind + NXY + 2*NX] - g[ind + NXY])/hy;
    gdata.fy7 = 0.5*(g[ind + NXY + 2*NX + 1] - g[ind + NXY + 1])/hy;
    
    gdata.fz0 = 0.5*(g[ind + NXY] - g[ind - NXY])/hz;
    gdata.fz1 = 0.5*(g[ind + NXY + 1] - g[ind - NXY + 1])/hz;
    gdata.fz2 = 0.5*(g[ind + NX + NXY] - g[ind + NX - NXY])/hz;
    gdata.fz3 = 0.5*(g[ind + NX + NXY + 1] - g[ind + NX - NXY + 1])/hz;
    gdata.fz4 = 0.5*(g[ind + 2*NXY] - g[ind])/hz;
    gdata.fz5 = 0.5*(g[ind + 2*NXY + 1] - g[ind + 1])/hz;
    gdata.fz6 = 0.5*(g[ind + NX + 2*NXY] - g[ind + NX])/hz;
    gdata.fz7 = 0.5*(g[ind + NX + 2*NXY + 1] - g[ind + NX + 1])/hz;
    
    return gdata;
}


/*************************************************/
struct myvector myrhs(struct myvector p,struct interpdata gdata,long ind) {
    double ax,ay,ax1,ay1,gnorm,az,az1;
    long i,j,k;
    
    struct myvector gu,b;
    
    i = (ind%NXY)%NX;
    j = (ind%NXY)/NX;
    k = ind/NXY;
    
    if( i > NX - 3 || j > NY - 3 || k > NZ - 3 || i < 2 || j < 2 || k < 2) {
        printf("In myrhs: y = (%.4e,%.4e,%.4e), (i,j) = (%li,%li,%li)\n",p.x,p.y,p.z,i,j,k);
        printf("Too close to the boundary\n");
        exit(1);
    }
    
    ax = (p.x - XMIN - i*hx)/hx; ax1 = 1.0 - ax;
    ay = (p.y - YMIN - j*hy)/hy; ay1 = 1.0 - ay;
    az = (p.z - ZMIN - k*hz)/hz; az1 = 1.0 - az;
    
    gu.x = -(((gdata.fx0*ax1 + gdata.fx1*ax)*ay1 + (gdata.fx2*ax1 + gdata.fx3*ax)*ay)*az1 + 
    		((gdata.fx4*ax1 + gdata.fx5*ax)*ay1 + (gdata.fx6*ax1 + gdata.fx7*ax)*ay)*az);
    gu.y = -(((gdata.fy0*ax1 + gdata.fy1*ax)*ay1 + (gdata.fy2*ax1 + gdata.fy3*ax)*ay)*az1 + 
    		((gdata.fy4*ax1 + gdata.fy5*ax)*ay1 + (gdata.fy6*ax1 + gdata.fy7*ax)*ay)*az);
    gu.z = -(((gdata.fz0*ax1 + gdata.fz1*ax)*ay1 + (gdata.fz2*ax1 + gdata.fz3*ax)*ay)*az1 + 
    		((gdata.fz4*ax1 + gdata.fz5*ax)*ay1 + (gdata.fz6*ax1 + gdata.fz7*ax)*ay)*az);
    b = myfield(p);
    gu = vec_difference(gu,b);
    
    gnorm = length_vec(gu);
    gu = a_times_vec(gu,1.0/gnorm);
    
    
    return gu;
}


/********************************************************/
/*** main ***/

 int main() {
  long i,j,k,ind,m; 
  double umax;
  clock_t CPUbegin;
  double cpu;
  FILE *fg;
  char fname[100];
  double q1,q2,q3,qc2,qc3;
  
 	  sprintf(fname,"LorenzQpot_rho%.2f.txt",rho);

	  printf("rho = %.2f\n",rho);
	  count = 0;
	  umax = 0.0;
	  N1ptu = 0;
	  N2ptu = 0;
	  N3ptu = 0;
	  N3call = 0;
	  N2call = 0;
	  param();
	  ipoint();
	  k = (2*K+1);
	  farlist = (long *)malloc(k*k*k*sizeof(long));
	  Nfar = far_neighbors_index_list(farlist);
	  CPUbegin=clock();
	  olim();
	  cpu = (clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
	  printf("N = %i: cputime of olim() = %g\n",NX,cpu); 
	  m = 0;
	  fg = fopen(fname,"w");
	  for( k=0; k<NZ; k+=SSF ) {
		 for( i=0; i<NX; i+=SSF ) {
			for( j=0; j<NY; j+=SSF ) {
				ind = k*NXY + j*NY + i;
				if( ms[ind] < 2 ) g[ind] = INFTY;
				else {
					umax = max(umax,g[ind]);
					m++;
				}
				fprintf(fg,"%.6e\n",g[ind]);
			}
		 }
	  }
	  fclose(fg);
	  printf("#finite data = %li,umax = %.4e,  NX = %i, NY = %i, NZ = %i, K = %i\n",
			m,umax,NX,NY,NZ,K);

 	  sprintf(fname,"parameters_rho%.2f.txt",rho);
	  
	  fg = fopen(fname,"w");
	  fprintf(fg,"%.12e\n%.12e\n%.12e\n%.12e\n%.12e\n%.12e\n%i\n%i\n%i\n",
	  		XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,NX,NY,NZ);
	  fprintf(fg,"%.14e\n%.14e\n%.14e\n",sigma,beta,rho);
	  fprintf(fg,"%.14e\n%.14e\n%.14e\n",x_ipoint.x,x_ipoint.y,x_ipoint.z);
	  fprintf(fg,"%i\n%i\n",K,SSF);
	  fclose(fg);		
					
	  q1 = N1ptu;
	  q2 = N2ptu;
	  q3 = N3ptu;
	  qc2 = N2call;
	  qc3 = N3call;		
	  printf("Per mesh point: <#1ptupdates> = %.2f\t<#2ptupdates> = %.2f\n<#3ptupdates> = %.2f\t<#2call> = %.2f\t<#3call> = %.2f\n",
			q1/m,q2/m,q3/m,qc2/m,qc3/m);
		
  return 0;
}  
