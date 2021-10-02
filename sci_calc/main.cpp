//
//  main.cpp
//  sci_calc
//
//  Created by 张子越 on 2020/9/24.
//  Copyright © 2020 张子越. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <unistd.h>

using namespace std;

typedef long double LD;
typedef long long LL;

const LD Pi=3.1415926535897932384626433832795028841971693993751;
const LD g=9.8;
const LD R=8.3144621;
const LD atm=101325;
const LD RES=1e-8;

const int MPL=5;//maxpolylevel

const int Solve_Depth=8;

//---------------------------------------------IO
void priLD(LD x){
	printf("%.4Lf\n",x);
}//defult .4
void priLD(LD x,int w){
	if (w>9) {
		printf("LD too long\n");
		return;
	}
	static char c[3];
	c[0]='0'+w;
	c[1]='\0';
	string wei=c;
	string s="%."+wei+"Lf\n";
	//printf("%s",s.c_str());
	printf(s.c_str(),x);
}
//---------------------------------------------/IO


//---------------------------------------------math
int ksm(int x,int p){
	int u=x,ans=1;
	while (p) {
		if (p&1) {
			ans*=u;
		}
		u*=u;
		p>>=1;
	}
	return ans;
}

LD facto(int x){
	LD tot=1;
	for (int i=1; i<=x; tot*=i++);
	return tot;
}
LD sqr(LD x){
	return x*x;
}
LD fen(LD x, LD y){
	return y/x;
}
LD frac(LD x,LD y){
	return x/y;
}
LD Volume_ball(LD r){
	return (fen(3, 4)*Pi*r*r*r);
}
LD SufArea_ball(LD r){
	return 4*Pi*r*r;
}
LD Area_circ(LD r){
	return Pi*r*r;
}
LD Cycl_circ(LD r){
	return 2*Pi*r;
}

LD func(LD x){
	x=pow(-0.4, x)/(x*facto(x));
	return x;
}
LD funcpri(LD x){
	LD esp=1e-6;
	static LD l,r;
	l=(func(x)-func(x-esp))/esp;
	r=(func(x+esp)-func(x))/esp;
	if (abs(l-r)>1) {
		exit(3);
	}
	printf("f'(%.4Lf)=%.4Lf\n",x,(l+r)/2);
	return (l+r)/2;
}
LD solvefunc(LD x0){
	static LD t;
	for (int i=0; i<Solve_Depth; ++i) {
		t=funcpri(x0);
		if (abs(t)<1e-5) printf("Waring: f' too small\n");
		x0=x0-(func(x0)/t);
		
	}
	return x0;
}
LD RiemannS(LD l,LD r,LD c,LD n){
	LD d=(r-l)/n;
	LD ans=0;
	for (int i=0; i<n; ++i) {
		ans+=func(l+i*d+c*d)*d;
	}
	return ans;
}
LD SimpsonInt(LD l,LD r,LD k,LD n){
	if (k!=n/2) {
		priLD(233.9999);
	}
	LD d=(r-l)/n;
	LD ans=func(l)*d,x=l;
	LD dk=(r-l)/k;
	for (int i=0; i<k; ++i) {
		ans+=(4*func(x+d+dk*i)+2*func(x+dk*(i+1)))*d;
	}
	ans-=func(r)*d;
	ans/=3;
	return ans;
}//2k=n
LD Biconst(LD n,int k){
	LD ans=1/facto(k);
	for (int i=0; i<k; ++i) {
		ans*=n-i;
	}
	return ans;
}
LD randNormal(LD mu=0, LD sigma=1){
	LD x_1,x_2,z;
	x_1=random();
	x_2=random();
	z=sqrt(log(x_1)*-2)*cos(2*Pi*x_2);
	return z*sigma+mu;
}
LD sign(LD x){
	return x/abs(x);
}
bool not0(LD x){
	return (abs(x)<RES)?0:1;
}
//---------------------------------------------/math

//---------------------------------------------/matrix
const int MATSIZE=8;
struct matrix {
	unsigned n,m;
	LD a[MATSIZE][MATSIZE];
};
matrix matmake(int n,int m, LD a[]){
	if (n>MATSIZE | m>MATSIZE) {
		printf("The matrix is too big.\n");
		exit(-1);
	}
	matrix mat;
	mat.n=n;
	mat.m=m;
	for (int i=0; i<n; ++i) {
		for (int j=0; j<m; ++j) {
			mat.a[i][j]=a[i*m+j];
		}
	}
	return mat;
}
matrix matmake(int n,int m){
	if (n>MATSIZE | m>MATSIZE) {
		printf("The matrix is too big.\n");
		exit(-1);
	}
	matrix mat;
	mat.n=n;
	mat.m=m;
	for (int i=0; i<n; ++i) {
		for (int j=0; j<m; ++j) {
			mat.a[i][j]=0;
		}
	}
	return mat;
}
matrix i_matrix(int n){
	matrix ans=matmake(n, n);
	for (int i=0; i<n; ++i) {
		ans.a[i][i]=1;
	}
	return ans;
}
matrix diag(int n,LD num[]){
	matrix ans=i_matrix(n);
	for (int i=0; i<n; ++i) {
		ans.a[i][i]=num[i];
	}
	return ans;
}
matrix copysize(matrix mat){
	matrix ans;
	ans.n=mat.n;
	ans.m=mat.m;
	return ans;
}
void view(matrix mat){
	printf("%d %d\n",mat.n,mat.m);
	printf("\\begin{bmatrix}\n");
	for (int i=0; i<mat.n; ++i) {
		for (int j=0; j<mat.m; ++j) {
			printf("%8.3Lf",mat.a[i][j]);
			if (j!=mat.m-1){
				printf("&");
			}
		}
		printf("\\\\\n");
	}
	printf("\\end{bmatrix}\n");
}
matrix operator+(matrix a, matrix b){
	matrix ans;
	if (!(a.n==b.n)|!(a.m==b.m)) {
		printf("Sizes do not match in matrix add.\n");
		exit(1);
	}
	ans=copysize(a);
	for (int i=0; i<a.n; ++i) {
		for (int j=0; j<a.m; ++j) {
			ans.a[i][j]=a.a[i][j]+b.a[i][j];
		}
	}
	return ans;
}
matrix operator*(matrix a, int c){
	matrix ans=copysize(a);
	for (int i=0;i<a.n;++i){
		for (int j=0;j<a.m;++j){
			ans.a[i][j]=a.a[i][j]*c;
		}
	}
	return ans;
}
matrix operator*(int c, matrix a){
	matrix ans=copysize(a);
	for (int i=0;i<a.n;++i){
		for (int j=0;j<a.m;++j){
			ans.a[i][j]=a.a[i][j]*c;
		}
	}
	return ans;
}
matrix operator*(LD c, matrix a){
	matrix ans=copysize(a);
	for (int i=0;i<a.n;++i){
		for (int j=0;j<a.m;++j){
			ans.a[i][j]=a.a[i][j]*c;
		}
	}
	return ans;
}
matrix operator*(double c, matrix a){
	matrix ans=copysize(a);
	for (int i=0;i<a.n;++i){
		for (int j=0;j<a.m;++j){
			ans.a[i][j]=a.a[i][j]*c;
		}
	}
	return ans;
}
matrix operator*(matrix a, LD c){
	matrix ans=copysize(a);
	for (int i=0;i<a.n;++i){
		for (int j=0;j<a.m;++j){
			ans.a[i][j]=a.a[i][j]*c;
		}
	}
	return ans;
}
matrix operator*(matrix a, double c){
	matrix ans=copysize(a);
	for (int i=0;i<a.n;++i){
		for (int j=0;j<a.m;++j){
			ans.a[i][j]=a.a[i][j]*c;
		}
	}
	return ans;
}
matrix operator-(matrix a, matrix b){
	return a + b*(-1);
}
matrix trans(matrix a){
	matrix ans=matmake(a.m, a.n);
	for (int i=0; i<a.n; ++i) {
		for (int j=0; j<a.m; ++j) {
			ans.a[j][i]=a.a[i][j];
		}
	}
	return ans;
}
matrix solve(matrix mat){
	matrix t=copysize(mat);
	int nn=mat.n;
	int mm=mat.m;
	for (int i=0; i<nn; ++i) {
		for (int j=0; j<mm; ++j) {
			t.a[i][j]=mat.a[i][j];
		}
	}
	int head=0;
	LD u;
	for (int i=0; i<nn; ++i) {
		if (abs(t.a[i][head])>RES) {
			for (int j=mm-1; j>=0; --j) {
				t.a[i][j]/=t.a[i][head];
			}//row normalize
			printf("%d %d\n",i,head);
//			view(t);
			for (int k=0; k<nn; ++k) if (k!=i) {
				u=t.a[k][head];
				printf("%d %Lf\n",k,u);
				for (int j=0; j<mm; ++j) {
					t.a[k][j]-=u*t.a[i][j];
//					view(t);
				}
//				view(t);
			}//kill var
//			view(t);
			head = 0;
		}else{
			++head;
			--i;
			if (head==mm) {
				break;
			}
		}
	}
	return t;
}
matrix leftproduct(matrix a,matrix b){
	matrix ans=matmake(a.n, b.m);
	int r=a.m;
	if (r!=b.n) {
		printf("error! The size do not match!");
		exit(1);
	}
	for (int i=0; i<a.n; ++i) {
		for (int j=0; j<b.m; ++j) {
			for (int k=0; k<r; ++k) {
				ans.a[i][j]+=a.a[i][k]*b.a[k][j];
			}
		}
	}
	return ans;
}
matrix operator*(matrix a,matrix b){
	return leftproduct(a, b);
}
matrix ksm(matrix x,int p){
	if (x.n!=x.m) {
		printf("error matrix for ksm.\n");
		exit(1);
	}
	matrix u=x;
	matrix ans=i_matrix(x.n);
	while (p) {
		if (p&1) {
			ans=ans*u;
		}
		u = u*u;
		p>>=1;
	}
	return ans;
}
void row_swap(matrix &a,int x,int y){
	for (int j=0; j<a.m; ++j) {
		swap(a.a[x][j],a.a[y][j]);
	}
}
void col_swap(matrix &a,int x,int y){
	for (int i=0; i<a.m; ++i) {
		swap(a.a[i][x],a.a[i][y]);
	}
}
matrix inve(matrix a){
	if (a.n!=a.m) {
		printf("The matrix is not a square matrix.");
	}
	matrix ans=i_matrix(a.n);
	int nn=a.n;
	LD u;
	bool done[MATSIZE];
	memset(done, 0, MATSIZE);
	for (int i=0; i<nn; ++i) {
		if (abs(a.a[i][i])<RES) {
			if (i==nn-1) {
				printf("The matrix is not invertable.\n");
				return i_matrix(nn);
			}
			for (int k=i+1; k<nn; ++k) {
				if (abs(a.a[k][i])>RES) {
					row_swap(a, i, k);
					row_swap(ans, i, k);
				}else if (k==nn-1){
					printf("The matrix is not invertable.\n");
					return i_matrix(nn);
				}
			}
		}
		u=a.a[i][i];
		for (int j=0; j<nn; ++j) {
			a.a[i][j]/=u;
			ans.a[i][j]/=u;
		}
//		view(a);
		for (int k=0; k<nn; ++k) if (k!=i) {
			u=a.a[k][i];
			for (int j=0; j<nn; ++j) {
				ans.a[k][j]-=u*ans.a[i][j];
				a.a[k][j]-=u*a.a[i][j];
			}
		}
//		view(a);
	}
	return ans;
}
matrix elem(int n,int i,int j,LD beta){
	matrix ans=i_matrix(n);
	ans.a[i][j]=beta;
	return ans;
}
matrix rev_elem(int n,int i,int j,LD beta){
	return elem(n, i, j, -beta);
}
//---------------------------------------------/matrix


//---------------------------------------------poly
struct int_ord_poly {
	int a[2*MPL+2];//^0 at a[MPL]
};

LD value_intcpoly(LD x,int_ord_poly p){
	LD px=1,s=0;
	for (int i=MPL; i<MPL*2+1; ++i,px*=x) {
		s+=p.a[i]*px;
	}
	if (x==0) {
		exit(2);
	}
	px=1/x;
	for (int i=MPL-1; i; --i,px/=x) {
		s+=p.a[i]*px;
	}
	return s;
}

int iord(int x){
	return x-MPL;
}

void write_intcpoly(int_ord_poly p){
	bool b=0;
	for (int i=2*MPL; i>=0; --i) {
		if (p.a[i]) {
			if (b) {
				if (p.a[i]>0) putchar('+');
			}
			b=1;
			if (i==MPL) {
				printf("%d ",p.a[i]);
				continue;
			}
			if (p.a[i]!=1) {
				printf("%dx^{%d} ",p.a[i],iord(i));
			}else{
				printf("x^{%d} ",iord(i));
			}
		}
	}
	putchar('\n');
}

int_ord_poly pe;
int_ord_poly diff_intcpoly(int_ord_poly p0){
	int_ord_poly pa=pe;
	for (int i=0; i<MPL*2; ++i) {
		pa.a[i]=(iord(i)+1)*p0.a[i+1];
	}
	return pa;
}

LD solve(int_ord_poly f,LD seed){
	int_ord_poly fp=diff_intcpoly(f);
	static LD t;
	for (int i=0; i<Solve_Depth; ++i) {
		t=value_intcpoly(seed, fp);
		if (abs(t)< 1e5) printf("Waring:f'' too small\n");
		seed=seed-(value_intcpoly(seed, f)/t);
	}
	return seed;
}
//---------------------------------------------/poly


//---------------------------------------------complex 3D
struct complex3 {
	LD x,y,z;
};
LD mag(complex3 a){
	return sqrt(sqr(a.x)+sqr(a.y)+sqr(a.z));
}
void view(complex3 a){
	printf("%Lf\n\t%Lf\n\t%Lf\n\t%Lf\n",mag(a),a.x,a.y,a.z);
}
complex3 operator-(complex3 a, complex3 b){
	complex3 ans;
	ans.x=a.x-b.x;
	ans.y=a.y-b.y;
	ans.z=a.z-b.z;
	return ans;
}
complex3 operator*(LD c, complex3 a){
	a.x*=c;
	a.y*=c;
	a.z*=c;
	return a;
}
LD dotp(complex3 a, complex3 b){
	return a.x*b.x+a.y*b.y+b.z*a.z;
}
complex3 crossp(complex3 a, complex3 b){
	complex3 ans;
	ans.x=a.y*b.z-a.z*b.y;
	ans.y=a.z*b.x-a.x*b.z;
	ans.z=a.x*b.y-a.y*b.x;
	return ans;
}
//---------------------------------------------/complex 3D


//---------------------------------------------complex 2D
struct complex2{
	LD x,y;
	LD rad,len;
};
LD dotp(complex2 a,complex2 b){
	return a.x*b.x+a.y*b.y;
};
LD crossp(complex2 a,complex2 b){
	return a.x*b.y-b.x*a.y;
};
void init_p(complex2 &c){
	c.x=c.len*cos(c.rad);
	c.y=c.len*sin(c.rad);
};
void init_d(complex2 &c){
	c.len=sqrt(dotp(c,c));
	c.rad=atan(c.y/c.x);
};
complex2 make_d(LD x,LD y){
	static complex2 c;
	c.x=x;
	c.y=y;
	init_d(c);
	return c;
}
complex2 make_p(LD l,LD rad){
	static complex2 c;
	c.len=l;
	c.rad=rad;
	init_p(c);
	return c;
}
LD mag(complex2 c){
	return sqrt(c.x*c.x+c.y*c.y);
}
void view(complex2 a){
	printf("%Lf\n\t%Lf\n\t%Lf\n",mag(a),a.x,a.y);
	printf("\t\t%Lf\n",a.rad);
}
complex2 operator+(complex2 a,complex2 b){
	a.x+=b.x;
	a.y+=b.y;
	init_d(a);
	return a;
};
complex2 operator+=(complex2 &a,complex2 b){
	a.x+=b.x;
	a.y+=b.y;
	init_d(a);
	return a;
}
complex2 operator-(complex2 a,complex2 b){
	a.x-=b.x;
	a.y-=b.y;
	init_d(a);
	return a;
};
complex2 operator-=(complex2 &a,complex2 b){
	a.x-=b.x;
	a.y-=b.y;
	init_d(a);
	return a;
}
complex2 operator*(complex2 a,complex2 b){
	a.len*=b.len;
	a.rad+=b.rad;
	init_p(a);
	return a;
};//complex product
complex2 operator*(complex2 a,LD l){
	a.x*=l;
	a.y*=l;
	init_d(a);
	return a;
}
complex2 operator*(complex2 a,int l){
	a.x*=l;
	a.y*=l;
	init_d(a);
	return a;
}
complex2 operator/(complex2 a,LD l){
	a.x/=l;
	a.y/=l;
	init_d(a);
	return a;
}
complex2 operator/(complex2 a,int l){
	a.x/=l;
	a.y/=l;
	init_d(a);
	return a;
}
complex2 operator*(LD l,complex2 a){
	a.x*=l;
	a.y*=l;
	init_d(a);
	return a;
}
complex2 operator*(int l,complex2 a){
	a.x*=l;
	a.y*=l;
	init_d(a);
	return a;
}
complex2 operator/(LD l,complex2 a){
	a.x/=l;
	a.y/=l;
	init_d(a);
	return a;
}
complex2 operator/(int l,complex2 a){
	a.x/=l;
	a.y/=l;
	init_d(a);
	return a;
}
complex2 operator*=(complex2 &a,complex2 b){
	a.len*=b.len;
	a.rad+=b.rad;
	init_p(a);
	return a;
}
LD operator^(complex2 a,complex2 b){
	return crossp(a, b);
}//cross product
LD operator&(complex2 a,complex2 b){
	return dotp(a, b);
}//dot product
complex2 operator/(complex2 a,complex2 b){
	a.rad-=b.rad;
	a.len/=b.rad;
	init_p(a);
	return a;
}
complex2 operator/=(complex2 &a,complex2 b){
	a.rad-=b.rad;
	a.len/=b.rad;
	init_p(a);
	return a;
}
complex2 operator*=(complex2 &a,LD l){
	a.x*=l;
	a.y*=l;
	init_d(a);
	return a;
}
complex2 operator*=(complex2 &a,int l){
	a.x*=l;
	a.y*=l;
	init_d(a);
	return a;
}
complex2 operator/=(complex2 &a,LD l){
	a.x/=l;
	a.y/=l;
	init_d(a);
	return a;
}
complex2 operator/=(complex2 &a,int l){
	a.x/=l;
	a.y/=l;
	init_d(a);
	return a;
}
complex2 ksm(complex2 x,int p){
	if (!p) {
		exit(2);
	}
	complex2 u=x,ans=x;
	p-=1;
	while (p) {
		if (p&1) {
			ans*=u;
		}
		u*=u;
		p>>=1;
	}
	return ans;
}
//----------------------------------------------/complex 2D


//----------------------------------------------system: centre of mass
struct pt_2D {
	LD m;
	complex2 r;
};

struct sys2D_m_cen {
	LD M,ri0,ri;
	complex2 r,rM;
};

void init(sys2D_m_cen s){
	s.M=0;
	s.ri=s.ri0=0;
	s.r.len=0;
	init_p(s.r);
	s.rM=s.r;
}
void add_pt(sys2D_m_cen &s, pt_2D p){
	s.rM+=p.r*p.m;
	s.M+=p.m;
	s.r=s.rM/s.M;
}
void operator+=(sys2D_m_cen &s, pt_2D p){
	add_pt(s, p);
}
//----------------------------------------------/system: centre of mass


//----------------------------------------------angles
LD rad_to_deg(LD rad){
	rad/=Pi;
	rad*=180;
	return rad;
}
LD deg_to_rad(LD d){
	d/=180;
	d*=Pi;
	return d;
}
LD Theta;
LD sin(){
	return sin(Theta);
}
LD cos(){
	return cos(Theta);
}
//----------------------------------------------/angles


//----------------------------------------------PHY
LD phy_M,phy_Ekt;
LD v_from_Ekt(){
	return sqrt(2*phy_Ekt/phy_M);
}

LD v1f(LD v1,LD m1,LD v2,LD m2){
	return (m1-m2)/(m1+m2)*v1+2*m2/(m1+m2)*v2;
}
LD v2f(LD v1,LD m1,LD v2,LD m2){
	return 2*m1/(m1+m2)*v1+(m2-m1)/(m1+m2)*v2;
}
complex2 ec_v1f(complex2 v1,LD m1,complex2 v2,LD m2){
	return (m1-m2)/(m1+m2)*v1+2*m2/(m1+m2)*v2;
}
complex2 ec_v2f(complex2 v1,LD m1,complex2 v2,LD m2){
	return 2*m1/(m1+m2)*v1+(m2-m1)/(m1+m2)*v2;
}
LD ec_v1f(LD v1,LD m1,LD v2,LD m2){
	return (m1-m2)/(m1+m2)*v1+2*m2/(m1+m2)*v2;
}
LD ec_v2f(LD v1,LD m1,LD v2,LD m2){
	return 2*m1/(m1+m2)*v1+(m2-m1)/(m1+m2)*v2;
}

//----------------------------------------------/PHY

//----------------------------------------------/STA

void array_add_1(int **a,int i,int x,int m){
	for (int j = 0; j<m; ++j) {
		a[i][j]+=x;
	}
}
void array_add_0(int **a,int j,int x,int n){
	for (int i = 0; i<n; ++i) {
		a[i][j]+=x;
	}
}
LD mpdfx[10],mpdfy[10];
LD mu_x,mu_y,sigmasqr_x,sigmasqr_y,cov_xy;

//----------------------------------------------/STA


int_ord_poly p,_dp;
complex2 ca,cb,cc,ct,f1,f2,f3;
LD lt1,lt2;
LD gx,gy;
LD phy_k,phy_h;
pt_2D p1,p2,p3;
sys2D_m_cen sys1;
LD ans;
int a[300][300];

LD ff(LD x){
	return 3*log(x-3)+log(x+1);
}

void Q(){
	int id=0;
	matrix a= matmake(3, 4, (LD[]){
		1,-1,2,1,
		2,1,1,8,
		1,1,0,5
	});
	matrix b= matmake(2, 2, (LD[]){
		1,0.5,
		0.5,-4
	});
	printf("----\tQ%d:\t\t----\n",id);
	view(a);
	view(b);
	
	printf("==  END  ==\n\n");
}

void Q1(){
	int id=1;
	LD h=0.5;
	matrix a= matmake(2, 2, (LD[]){
		1,2,
		1,-2
	});
	matrix b= matmake(2, 2, (LD[]){
		h,h,
		h,-h
	});
	printf("----\tQ%d:\t\t----\n",id);
	view(a);
	view(b);
	view(ksm(a*b, 10));
	printf("==  END  ==\n\n");
}
void Q2(){
	int id=2;
	matrix a= matmake(1, 2, (LD[]){
		1,1
	});
	a=trans(a);
	matrix b= matmake(1, 2, (LD[]){
		1,-1
	});
	b=trans(b);
	printf("----\tQ%d:\t\t----\n",id);
	matrix c[5];
	view(c[1]=a*trans(a));
	view(c[2]=a*trans(b));
	view(c[3]=b*trans(a));
	view(c[4]=b*trans(b));
	view((c[1]+c[2]+c[3]+c[4])*0.25);
	view((c[1]-c[2]+c[3]-c[4])*0.25);
	view((c[1]+c[2]-c[3]-c[4])*0.25);
	view((c[1]-c[2]-c[3]+c[4])*0.25);
	printf("==END==\n\n");
}
void Q3(){
	int id=0;
	matrix u= matmake(4, 5, (LD[]){
		1,0,2,0,-1,
		0,1,3,0,-2,
		0,0,0,1,5,
		0,0,0,0,0
	});
	matrix x0= matmake(5, 1, (LD[]){
		3,2,0,2,0
	});
	matrix a= matmake(4, 5, (LD[]){
		2,-1,1,-2,-10,
		1,2,8,-1,-10,
		-3,3,3,3,12,
		-2,1,-1,4,20
	});
	printf("----\tQ%d:\t\t----\n",id);
	view(x0);
	view(u);
	view(u*x0);// =P^{-1}\vec b
	view(a);
	view(a*x0);
//	view(solve(a));
	printf("==  END  ==\n\n");
}
void Q4(){
	int id=4;
	matrix a= i_matrix(2);
	matrix b1= matmake(2, 2, (LD[]){
		2,0,
		0,1
	});
	matrix b2= matmake(2, 2, (LD[]){
		1,1,
		0,1
	});
	matrix b3= matmake(2, 2, (LD[]){
		1,0,
		3,1
	});
	matrix b4= matmake(2, 2, (LD[]){
		1,0,
		0,1
	});
	printf("----\tQ%d:\t\t----\n",id);
	view(b1);
	view(b2);
	view(b3);
	view(b4*b3*b2*b1*a);
	printf("==  END  ==\n\n");
}
void Q5(){
	int id=5;
	matrix a= matmake(3, 3, (LD[]){
		3,4,1,
		-2,0,1,
		1,2,2
	});
	matrix b= diag(3, (LD[]){1,-2,4});
	printf("----\tQ%d:\t\t----\n",id);
	view(a);
	view(b);
	view(b*a);
	view(a*b);
	printf("==  END  ==\n\n");
}

int main() {
//	LD h=0.5;
//	matrix a= matmake(3, 4, (LD[]){
//		1,-1,2,1,
//		2,1,1,8,
//		1,1,0,5
//	});
//	matrix b= matmake(2, 2, (LD[]){
//		1,0.5,
//		0.5,-4
//	});
//	view(a);
//	view(solve(a));
//	view(b);
//	matrix c=leftproduct(a+b, a-b);
//	view(c);
//	view(a*b*a);
//	view(b*a);
//	view(a*b);
	
//	Q1();
//	Q2();
//	Q3();
//	Q4();
//	Q5();
	
	return 0;
}
