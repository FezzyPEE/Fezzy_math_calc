//
//  fezzy_matrix.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_matrix_h
#define fezzy_matrix_h

const int MATSIZE=8;
struct matrix {
	unsigned n,m,rk=-1;
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
	printf("$%d\\times%d$ matrix:\n",mat.n,mat.m);
	printf("$$\n\\begin{bmatrix}\n");
	for (int i=0; i<mat.n; ++i) {
		for (int j=0; j<mat.m; ++j) {
			printf("%8.3Lf",mat.a[i][j]);
			if (j!=mat.m-1){
				printf("&");
			}
		}
		printf("\\\\\n");
	}
	printf("\\end{bmatrix}\n$$\n");
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
matrix matget(matrix mat,int x,int y,int n,int m){
	matrix ans = matmake(n, m);
	for (int i=x; i<x+n; ++i) {
		for (int j=y; j<y+m; ++j) {
			ans.a[i-x][j-y]=mat.a[i][j];
		}
	}
	return ans;
}
matrix matbind_r(matrix mup, matrix mdown){
	if (mup.m!=mdown.m) {
		printf("!!!Error binding.!!!\n");
		return mup;
	}
	int un=mup.n,dn=mdown.n,mm=mup.m;
	matrix ans=matmake(un+dn, mm);
	for (int i=0; i<un; ++i) {
		for (int j=0; j<mm; ++j) {
			ans.a[i][j]=mup.a[i][j];
		}
	}
	for (int i=0; i<dn; ++i) {
		for (int j=0; j<mm; ++j) {
			ans.a[un+i][j]=mdown.a[i][j];
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
//			printf("%d %d\n",i,head);
//			view(t);
			for (int k=0; k<nn; ++k) if (k!=i) {
				u=t.a[k][head];
//				printf("%d %Lf\n",k,u);
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
				return solve(a);
			}
			for (int k=i+1; k<nn; ++k) {
				if (abs(a.a[k][i])>RES) {
					row_swap(a, i, k);
					row_swap(ans, i, k);
				}else if (k==nn-1){
					printf("The matrix is not invertable.\n");
					return solve(a);
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
int matrank(matrix &mat){
//	if (mat.rk!=-1) {
//		return mat.rk;
//	}
	matrix a=solve(mat);
	int i;
	bool b=0;
	for (i=0; i<a.n; ++i) {
		b=0;
		for (int j=0; j<a.m; ++j) {
			b|=(a.a[i][j]>RES);
		}
		if (!b) {
			break;
		}
	}
	mat.rk=i;
	return i;
}
matrix column_space(matrix mat){
	matrank(mat);
	if (mat.m<=mat.rk) {
		return matmake(1, 1, (LD[]){0});
	}
	return matbind_r(-1*matget(solve(mat), 0, mat.rk, mat.rk, mat.m-mat.rk), i_matrix(mat.m-mat.rk));
}
matrix row_basis(matrix mat){
	matrix ans=matmake(0, mat.m),tmp;
	for (int i=0; i<mat.n; ++i) {
		tmp=matbind_r(ans, matget(mat, i, 0, 1, mat.m));
		if (matrank(tmp)>matrank(ans)) {
			ans=tmp;
		}
	}
	return ans;;
}


#endif /* fezzy_matrix_h */
