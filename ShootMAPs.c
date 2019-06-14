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
#define E3 0.333333333333333
#define e6 0.166666666666666
// Lorenz system parameters
#define sigma  10.0
#define beta  2.666666666666667
#define rho 20.0

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
struct myvector myrhs_lin(struct myvector p);
/***************************************/

long NX,NY,NZ,nx1,ny1,nz1,NXY,NXYZ,nx2,ny2,nz2;
long K, SSF;
double hx,hy,hz,h;
double *g; /* function to be computed */
struct myvector x_ShootMAP;
long Iindex; // the index of the equilibrium point with respect to which the quasi-potential will be computed

// variables for the potential
double XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX; // define the computational domain

// variables for shooting an instanton;
struct myvector *path, x_ipoint;
long Npath;
long Npathmax;
struct myvector x_ShootMAP; // the point from which the instanton is shot backward to the attractor

struct mymatrix3 Qmatrix, Rmatrix, J, Jrev;


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
    double r1 = rho - 1, sq; 
    
    sq = sqrt(beta*r1); // sqrt(beta*(rho - 1)
    if( rho > 1.0 ) {
		// positive equilibrium point of Lorenz'63
		x_ipoint.x = sq;
		x_ipoint.y = sq;
		x_ipoint.z = r1;
		J.a11 = -sigma;	J.a12 = sigma; J.a13 = 0.0;
		J.a21 = 1.0; J.a22 = -1.0; J.a23 = -sq;
		J.a31 = sq; J.a32 = sq; J.a33 = -beta;

	}
	else {
		// the origin is asymptotically stable equilibrium
		x_ipoint.x = 0.0;
		x_ipoint.y = 0.0;
		x_ipoint.z = 0.0;
		J.a11 = -sigma;	J.a12 = sigma; J.a13 = 0.0;
		J.a21 = -rho; J.a22 = -1.0; J.a23 = 0.0;
		J.a31 = 0.0; J.a32 = 0.0; J.a33 = -beta;
	}	


	nx1 = NX - 1;
	ny1 = NY - 1;
	nz1 = NZ - 1;
	NXY = NX*NY;
	NXYZ = NX*NY*NZ;
	nx2 = NX - 2;
	ny2 = NY - 2;
	nz2 = NZ - 2;
	
	hx = (XMAX - XMIN)/nx1;
	hy = (YMAX - YMIN)/ny1;
	hz = (ZMAX - ZMIN)/nz1;
	h = sqrt(hx*hx + hy*hy + hz*hz);
	Npathmax = 1000*max(max(NX,NY),NZ);
	
	Qmatrix = findQmatrix(sigma,beta,rho); // set up struct matrixS3 Qmatrix
	Rmatrix = matrix_sum(J,Qmatrix);
	Jrev = matrix_sum(Rmatrix,Qmatrix);
	print_matrix(Jrev,"Jrev");

	
	printf("x_ipoint: Iindex = %li, %.4e, %.4e, %.4e\n",Iindex,x_ipoint.x,x_ipoint.y,x_ipoint.z);
    
}




/*-------------*/

struct myvector getpoint(long ind) {
    struct myvector l;
    
    l.x = hx*((ind%NXY)%NX) + XMIN;
    l.y = hy*((ind%NXY)/NX) + YMIN;
    l.z = hz*(ind/NXY) + ZMIN;
    return l;
}





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
	
/********** shoot MAP **************/
void ShootMAP() {
    long j,indprev = -1,ind;
    struct myvector y,k1,k2,k3,k4,z;
    struct interpdata gdata;
    double dt = 0.1*h,dt2;
    double Rswitch = 1.0;
    
    path = (struct myvector *)malloc(Npathmax*sizeof(struct myvector));
    
    dt2 = 0.5*dt;
    path[0] = x_ShootMAP;
    j = 0;
    while( length_vec(vec_difference(x_ipoint,path[j])) > Rswitch ) {
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
        
        path[j+1] = vec_sum(path[j],a_times_vec(vec_sum(vec_sum(k1,a_times_vec(vec_sum(k2,k3),2.0)),k4),dt*e6));
        j++;
        printf("j = %li, path: %.4e, %.4e, %.4e, distance: %.4e, step = %.4e\n",
        		j,path[j].x,path[j].y,path[j].z,length_vec(vec_difference(x_ipoint,path[j])),length_vec(vec_difference(path[j],path[j-1])));
        if( j >= Npathmax || path[j].x < XMIN || path[j].x > XMAX || path[j].y < YMIN || path[j].y > YMAX || path[j].z < ZMIN || path[j].z > ZMAX ) {
            printf("Failed to compute MAP: the maximal allowed path length is reached: j = %li\n",j);
             Npath = j;
            return;
        }
    }
    dt = 0.01*h;
    dt2 = 0.5*dt;
    z = vec_difference(path[j],x_ipoint);
    while( length_vec(z) > h ) {
        y = z;
        k1 = myrhs_lin(y);
        
        y = vec_sum(z,a_times_vec(k1,dt2));
        k2 = myrhs_lin(y);
        
        y = vec_sum(z,a_times_vec(k2,dt2));
        k3 = myrhs_lin(y);
        
        y = vec_sum(z,a_times_vec(k3,dt));
        k4 = myrhs_lin(y);
        
        z = vec_sum(z,a_times_vec(vec_sum(vec_sum(k1,a_times_vec(vec_sum(k2,k3),2.0)),k4),dt*e6));
        path[j + 1] = vec_sum(z,x_ipoint);
        j++;
        printf("j = %li, path: %.4e, %.4e, %.4e, distance: %.4e, step = %.4e\n",
        		j,path[j].x,path[j].y,path[j].z,length_vec(vec_difference(x_ipoint,path[j])),length_vec(vec_difference(path[j],path[j-1])));
        if( j >= Npathmax || path[j].x < XMIN || path[j].x > XMAX || path[j].y < YMIN || path[j].y > YMAX || path[j].z < ZMIN || path[j].z > ZMAX ) {
            printf("Failed to compute MAP: the maximal allowed path length is reached: j = %li\n",j);
             Npath = j;
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

/*************************************************/

struct myvector myrhs_lin(struct myvector p) {
	struct myvector rhs;
	
	rhs = a_times_vec(matrix_vector(Jrev,p),-1.0);
	rhs = a_times_vec(rhs,1.0/length_vec(rhs));
	
	return rhs;
}


/********************************************************/
/*** main ***/

 int main() {
  long i,j,k,ind,m; 
  double qpot,x0,y0,z0,a1,a2,a3;
  FILE *fg,*fid;
  char fname[100];
  
  
  
 	  sprintf(fname,"parameters_rho%.2f.txt",rho);
  
  	  fg = fopen(fname,"r");
	  fscanf(fg,"%le\n%le\n%le\n%le\n%le\n%le\n%li\n%li\n%li\n",
	  		&XMIN,&XMAX,&YMIN,&YMAX,&ZMIN,&ZMAX,&NX,&NY,&NZ);
	  printf("XMIN = %.4e\nXMAX = %.4e\nYMIN = %.4e\nYMAX = %.4e\nZMIN = %.4e\nZMAX = %4e\nNX = %li\nNY = %li\nNZ = %li\n",
	  		XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,NX,NY,NZ);
	  fscanf(fg,"%le\n%le\n%le\n%le\n%le\n%le\n%li\n%li\n",
	  		&a1,&a2,&a3,&x0,&y0,&z0,&K,&SSF);
	  printf("sigma = %.1f\nbeta = %.4e\nrho = %.2f\nxe = %.4e\nye = %.4e\nze = %.4e\nK = %li\nSSF = %li\n",
	  		sigma,beta,rho,x0,y0,z0,K,SSF);
	  fclose(fg);		
	  printf("done reading parameters\n");

	  // take subsampling into account
	  NX = (NX - 1)/SSF + 1;			
	  NY = (NY - 1)/SSF + 1;			
	  NZ = (NZ - 1)/SSF + 1;			

 	  sprintf(fname,"LorenzQpot_rho%.2f.txt",rho);

	  printf("rho = %.2f\n",rho);
	  param();
	  
	  printf("h = %.4e\n",h);

	  g = (double*)malloc(NXYZ*sizeof(double));
	  
	  m = 0;
	  fg = fopen(fname,"r");
	  if( fg == NULL ) {
	  	  printf("Cannot find the file %s\n",fname);
	  	  exit(1);
	  }
	  printf("Reading data for the quasi-potential ...");
	  for( k=0; k<NZ; k++ ) {
		 for( i=0; i<NX; i++ ) {
			for( j=0; j<NY; j++ ) {
				ind = k*NXY + j*NY + i;
				fscanf(fg,"%le\n",&qpot);
				g[ind] = qpot;
			}
		 }
	  }
	  fclose(fg);
	  printf(" done\n");

	       
	  sprintf(fname,"MAP_rho%.2f.txt",rho);
	  fid = fopen(fname,"w");
	  if( rho > 1 && rho < 13 ) {
		  x_ShootMAP.x = 0.0;
		  x_ShootMAP.y = 0.1;
		  x_ShootMAP.z = 0.0;
	  }
	  if( fabs(rho - 15.0) < 1e-14 ) {
		  x_ShootMAP.x = 1.023962e+01;
		  x_ShootMAP.y = 5.972989e+00;
		  x_ShootMAP.z = 2.294375e+01;
	  }
	  if( fabs(rho - 20.0) < 1e-14 ) {
		  x_ShootMAP.x = 2.930613311969992;
		  x_ShootMAP.y = 3.999148384027791;
		  x_ShootMAP.z = 10.674805059437325;
	  }
	  
	  
	  ShootMAP();
	  for( j = 0; j < Npath; j++ ) {
		  fprintf(fid,"%.4e\t%.4e\t%.4e\n",path[j].x,path[j].y,path[j].z);
	  }
	  fclose(fid);



  return 0;
}  
