#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "findQmatrix.h"
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
// Lorenz system parameters
// #define sigma  10.0
// #define beta  2.666666666666667
// #define rho 24.4




struct myvector {
    double x;
    double y;
    double z;
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

struct mymatrix2 {
	double a11;
	double a12;
	double a21;
	double a22;
};	


struct matrixS3 {
	double a11;
	double a12;
	double a13;
	double a22;
	double a23;
	double a33;
};	

struct mymatrix3 findQmatrix(double sigma,double beta,double rho);
double wilkinson(double a,double b,double *coef); // Wilkinson's hybrid bisection/secant nonlinear solver
double polyval(double x,double *coef);
double length_vec(struct myvector x);
double norm2squared(struct myvector x);
struct myvector a_times_vec(struct myvector v,double a);
double dot_product(struct myvector a,struct myvector b);
struct mymatrix3 make_matrix(struct myvector w1,struct myvector w2,struct myvector w3);
struct myvector cross_product(struct myvector a,struct myvector b);
struct myvector matrix_vector(struct mymatrix3 A,struct myvector x);
struct mymatrix2 qpot2D(struct mymatrix2 M);
struct mymatrix3 find_ONB_in_unstable_manifold(struct mymatrix3 M);
struct mymatrix3 matrix_product(struct mymatrix3 A,struct mymatrix3 B);
struct myvector vec_difference(struct myvector v1,struct myvector v2);
void print_matrix(struct mymatrix3 U,char *name);
struct mymatrix3 matrix_sum(struct mymatrix3 A,struct mymatrix3 B);


//-----------
double polyval(double x,double *coef){
	double pval;
	int i;
	
	pval = 1.0;
	for( i = 0; i < 3; i++ ) {
		pval *= x;
		pval += coef[i];
	}
	
	return pval;	
}
//--------------

double wilkinson(double a,double b,double *coef) {
    double c, fa, fb, fc, d, fd, dd, df, dm, ds,t;
    double NONLIN_TOL = 4e-16;
    int iter = 0, itermax = 100;
        
    c = a;
    fa = polyval(a,coef);
    fb = polyval(b,coef);
    fc = fa;
    
    if( (fa > 0 && fb > 0 ) || (fa < 0 && fb < 0) ) {
        //	 root is not bracketed
        printf("wilkinson: root is not bracketed: a == %.4e, b = %.4e, fa = %.4e, fb = %.4e\n",a,b,fa,fb);
        exit(1);
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
        fd = polyval(d,coef);
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
        
    return b;
}


//--------------------
// 
// double length_vec(struct myvector x) {
//     return sqrt(norm2squared(x));
// }
// 
// //--------------------
// 
// double norm2squared(struct myvector x) {
//     return x.x*x.x + x.y*x.y + x.z*x.z;
// }
// 
// /*****************/
// struct myvector a_times_vec(struct myvector v,double a) {
//     struct myvector av;
//     
//     av.x = a*v.x;
//     av.y = a*v.y;
//     av.z = a*v.z;
//     
//     return av;
// }
// 
// double dot_product(struct myvector a,struct myvector b) {
//     return a.x*b.x+a.y*b.y+a.z*b.z;
// }
// 
// struct myvector cross_product(struct myvector a,struct myvector b) {
//     struct myvector v;
//     
//     v.x = a.y*b.z - a.z*b.y;
//     v.y = a.z*b.x - a.x*b.z;
//     v.z = a.x*b.y - a.y*b.x;
//     
//     return v;
// }

struct myvector matrix_vector(struct mymatrix3 A,struct myvector v) {
	struct myvector y;
	
	y.x = A.a11*v.x + A.a12*v.y + A.a13*v.z;
	y.y = A.a21*v.x + A.a22*v.y + A.a23*v.z;
	y.z = A.a31*v.x + A.a32*v.y + A.a33*v.z;
	
	return y;
}	

struct mymatrix2 qpot2D(struct mymatrix2 M) {
	struct mymatrix2 Q;
	double aux1,aux2,aux1sq,aux2sq,a,b;
	
	aux1 = M.a11 + M.a22; aux1sq = aux1*aux1;
	aux2 = M.a21 - M.a12;  aux2sq = aux2*aux2;
	a = aux1sq/(aux1sq + aux2sq);
	b = aux1*aux2/(aux1sq + aux2sq);
	Q.a11 = -(a*M.a11 + b*M.a21);
	Q.a12 = -(a*M.a12 + b*M.a22);
	Q.a21 = Q.a12;
	Q.a22 = -(a*M.a22 - b*M.a12);

	return Q;
}	

// struct myvector vec_difference(struct myvector v1,struct myvector v2) {
//     struct myvector v;
//     
//     v.x = v1.x - v2.x;
//     v.y = v1.y - v2.y;
//     v.z = v1.z - v2.z;
//     return v;
// }

struct mymatrix3 find_ONB_in_unstable_manifold(struct mymatrix3 M) {
	struct mymatrix3 P;
	double det,a1 = M.a12,a2 = M.a13,b1,b2;
	double c11 = M.a22 - M.a11, c12 = M.a32, c21 = M.a23, c22 = M.a33 - M.a11;
	struct myvector w1,w2,w3;
	// The plan:
	// M is of the form [M.a11,[a1,a2];[0;0],C] 
	// M[[b1,b2];I] = [[b1,b2];I]C
	// solve for b1,b2
	// orthogonalize [[b1,b2];I] and complete with the third vector to get 
	det = c11*c22 - c12*c21; // det(C^\top  - M11*I)
	b1 = (a1*c22 - c12*a2)/det;
	b2 = (c11*a2 - c21*a1)/det;
	w1.x = b1; w1.y = 1.0; w1.z = 0.0;
	w1 = a_times_vec(w1,1.0/length_vec(w1)); // normalize w1;
	w2.x = b2; w2.y = 0.0; w2.z = 1.0;
	w2 = vec_difference(w2,a_times_vec(w1,dot_product(w1,w2)));
	w2 = a_times_vec(w2,1.0/length_vec(w2)); // normalize w3;
	w3 = cross_product(w1,w2);
	w3 = a_times_vec(w3,1.0/length_vec(w3)); // normalize w3;
	P = make_matrix(w1,w2,w3);
	
	return P;
}

struct mymatrix3 matrix_product(struct mymatrix3 A,struct mymatrix3 B) {
	struct mymatrix3 C;
	
	C.a11 = A.a11*B.a11 + A.a12*B.a21 + A.a13*B.a31;
	C.a12 = A.a11*B.a12 + A.a12*B.a22 + A.a13*B.a32;
	C.a13 = A.a11*B.a13 + A.a12*B.a23 + A.a13*B.a33;
	
	C.a21 = A.a21*B.a11 + A.a22*B.a21 + A.a23*B.a31;
	C.a22 = A.a21*B.a12 + A.a22*B.a22 + A.a23*B.a32;
	C.a23 = A.a21*B.a13 + A.a22*B.a23 + A.a23*B.a33;

	C.a31 = A.a31*B.a11 + A.a32*B.a21 + A.a33*B.a31;
	C.a32 = A.a31*B.a12 + A.a32*B.a22 + A.a33*B.a32;
	C.a33 = A.a31*B.a13 + A.a32*B.a23 + A.a33*B.a33;
	
	return C;
	
}	

struct mymatrix3 transpose(struct mymatrix3 A) {
	struct mymatrix3 C;
	
	C.a11 = A.a11; C.a22 = A.a22; C.a33 = A.a33;
	C.a12 = A.a21; C.a13 = A.a31; C.a23 = A.a32;
	C.a21 = A.a12; C.a31 = A.a13; C.a32 = A.a23;
	
	return C;
}

struct mymatrix3 make_matrix(struct myvector w1,struct myvector w2,struct myvector w3) {
	struct mymatrix3 P;
	
	// P = [w1,w2,w3]
	P.a11 = w1.x; P.a12 = w2.x; P.a13 = w3.x;
	P.a21 = w1.y; P.a22 = w2.y; P.a23 = w3.y;
	P.a31 = w1.z; P.a32 = w2.z; P.a33 = w3.z;

	return P;
}

struct mymatrix3 matrix_sum(struct mymatrix3 A,struct mymatrix3 B) {
	struct mymatrix3 C;
	
	C.a11 = A.a11 + B.a11;
	C.a12 = A.a12 + B.a12;
	C.a13 = A.a13 + B.a13;
	
	C.a21 = A.a21 + B.a21;
	C.a22 = A.a22 + B.a22;
	C.a23 = A.a23 + B.a23;

	C.a31 = A.a31 + B.a31;
	C.a32 = A.a32 + B.a32;
	C.a33 = A.a33 + B.a33;

	return C;
}



void print_matrix(struct mymatrix3 U,char *name ) {
	printf("Matrix %s:\n",name);
	printf("%.14e\t%.14e\t%.14e\n%.14e\t%.14e\t%.14e\n%.14e\t%.14e\t%.14e\n",
		U.a11,U.a12,U.a13,U.a21,U.a22,U.a23,U.a31,U.a32,U.a33);	
}

//--------------
struct mymatrix3 findQmatrix(double sigma,double beta,double rho) {
	// coefficients of the monic characteristic polynomial of Lorenz Jacobian
	double *coef;
	double a,b,fa,fb,root,x1,x2,y1,y2,z,tr;
	int i, imax = 100;
	double aux, sq = sqrt(beta*(rho - 1.0));
	struct myvector v,o1,o2,a1,a2,w1,w2;
	struct mymatrix3 J,J1,J2,V,R,U,O,Psi;
	struct mymatrix2 M,Q,WUW,WRW;
	double lo;
	
	if( rho > 1.0 ) {
		coef = (double *)malloc(3.0*sizeof(double));
		coef[0] = beta + sigma + 1.0;
		coef[1] = beta*(sigma + rho);
		coef[2] =  2.0*beta*sigma*(rho - 1.0);
		// bracket the root of the characteristic polynomial
		a = -20.0;
		b = -19.0;
		fa = polyval(a,coef);
		fb = polyval(b,coef);
		i = 0;
		while( fa*fb > 0  && i < imax ) {
			a = b;
			fa = fb;
			b += 1.0;
			fb = polyval(b,coef);
			i++;
		}
		if( i == imax ) {
			printf("A root of the characteristic polynomial	is not found\n");
			exit(1);
		}
		// once the root is bracketed, find it
		root = wilkinson(a,b,coef);	
		// find the corresponding eigenvector
		// it will be the last row of the matrix product
		// putting J^\top into the echelon form
		aux = 1.0 - (1.0 + root)*(sigma + root)/sigma;
		v.x = sq/aux;
		v.y = sq*(sigma + root)/(sigma*aux);
		v.z = 1.0;
		v = a_times_vec(v,1.0/length_vec(v));
		// complete v to an orthogonal basis
		do {
			o1.x = rand();
			o1.y = rand();
			o1.z = rand();
			o1 = a_times_vec(o1,1.0/length_vec(o1));
			o1 = vec_difference(o1,a_times_vec(v,dot_product(o1,v)));
			lo = length_vec(o1);
		}
		while( lo < 1e-3 );
		o1 = a_times_vec(o1,1.0/lo);
		o2 = cross_product(v,o1);
		o2 = a_times_vec(o2,1.0/length_vec(o2));
		// Put J to a Schur form
		J.a11 = -sigma;	J.a12 = sigma; J.a13 = 0.0;
		J.a21 = 1.0; J.a22 = -1.0; J.a23 = -sq;
		J.a31 = sq; J.a32 = sq; J.a33 = -beta;
		// Schur form of J is J1 = O^\top J O, where O = [v, o1, o2]
		// the first column of J1 is [root;0;0]
		O = make_matrix(v,o1,o2);
		J1 = matrix_product(transpose(O),matrix_product(J,O));
	}
	else {
		J.a11 = -sigma;	J.a12 = sigma; J.a13 = 0.0;
		J.a21 = -rho; J.a22 = -1.0; J.a23 = 0.0;
		J.a31 = 0.0; J.a32 = 0.0; J.a33 = -beta;
		O.a11 = 0.0; O.a12 = 0.0; O.a13 = 1.0;
		O.a21 = 0.0; O.a22 = 1.0; O.a23 = 0.0;
		O.a31 = 1.0; O.a32 = 0.0; O.a33 = 0.0;
		J1 = matrix_product(transpose(O),matrix_product(J,O));
	}
	// find the quasi-potential decomposition of the submatrix J1(2:3,2:3)
	M.a11 = J1.a22; M.a12 = J1.a23; M.a21 = J1.a32; M.a22 = J1.a33;
	Q = qpot2D(M);
	// define J2 by reverting Q
	J2.a11 = J1.a11; J2.a12 = J1.a12; J2.a13 = J1.a13;
	J2.a21 = 0.0; J2.a22 = 2.0*Q.a11 + M.a11; J2.a23 = 2.0*Q.a12 + M.a12;
	J2.a31 = 0.0; J2.a32 = 2.0*Q.a21 + M.a21; J2.a33 = 2.0*Q.a22 + M.a22;
	// find on ONB in the unstable manifold of J2
	Psi = find_ONB_in_unstable_manifold(J2);
	// V = false quasi-potential matrix [0,0,0;[[0;0],Q]];
	// R = false rotational matrix [J1(1,1:3);[[0;0],M + Q]];
	// The quasi-potential matrix for J1 and V restricted to unstable subspace of J2 must coincide
	V.a11 = 0.0; V.a12 = 0.0; V.a13 = 0.0;
	V.a21 = 0.0; V.a22 = Q.a11; V.a23 = Q.a12;
	V.a31 = 0.0; V.a32 = Q.a21; V.a33 = Q.a22;
	R.a11 = J1.a11; R.a12 = J1.a12; R.a13 = J1.a13;
	R.a21 = 0.0; R.a22 = Q.a11 + M.a11; R.a23 = Q.a12 + M.a12;
	R.a31 = 0.0; R.a32 = Q.a21 + M.a21; R.a33 = Q.a22 + M.a22;
	w1.x = Psi.a11; w1.y = Psi.a21; w1.z = Psi.a31;
	w2.x = Psi.a12; w2.y = Psi.a22; w2.z = Psi.a32;
	// compute WUW = V restricted to the unstable subspace of J2
	a1 = matrix_vector(V,w1); a2 = matrix_vector(V,w2);
	WUW.a11 = dot_product(w1,a1); WUW.a12 = dot_product(w1,a2);
	WUW.a21 = WUW.a12; WUW.a22 = dot_product(w2,a2);
	// compute WRW = R restricted to the unstable subspace of J2
	a1 = matrix_vector(R,w1); a2 = matrix_vector(R,w2);
	WRW.a11 = dot_product(w1,a1); WRW.a12 = dot_product(w1,a2);
	WRW.a21 = dot_product(w2,a1); WRW.a22 = dot_product(w2,a2);
	// set up a system of linear equations
	J2 = matrix_product(transpose(Psi),matrix_product(J1,Psi));
	tr = -WRW.a11 - WRW.a22;
	z = tr - J2.a33;
	aux = z + tr;
	M.a11 = WUW.a11 + WRW.a11 + aux;
	M.a12 = WUW.a12 + WRW.a21;
	M.a21 = WUW.a21 + WRW.a12;
	M.a22 = WUW.a22 + WRW.a22 + aux;
	y1 = -WUW.a11*J2.a13 - WUW.a12*J2.a23 - z*J2.a31;
	y2 = -WUW.a21*J2.a13 - WUW.a22*J2.a23 - z*J2.a32;
	aux = M.a11*M.a22 - M.a21*M.a12;
	x1 = (y1*M.a22 - y2*M.a12)/aux;
	x2 = (M.a11*y2 - M.a21*y1)/aux;
	// write the quasi-potential matrix in the current coordinates
	U.a11 = WUW.a11; U.a12 = WUW.a12; U.a13 = x1;
	U.a21 = WUW.a21; U.a22 = WUW.a22; U.a23 = x2;
	U.a31 = x1; U.a32 = x2; U.a33 = z;
	// return to the original coordinates: U = (O*Psi)*U*(O*Psi)';
	Psi = matrix_product(O,Psi);
	U = matrix_product(Psi,matrix_product(U,transpose(Psi)));
	print_matrix(U,"U");
	// verification
	R = matrix_sum(J,U);
	O = matrix_product(U,R);
	V = matrix_sum(O,transpose(O));
	print_matrix(V,"U*R + transpose(U*R)");
		
	return U;
}

