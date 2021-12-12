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
char sep='\n';
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
//	printf("$%d\\times%d$ matrix:\n",mat.n,mat.m);
	printf("\\begin{bmatrix}\n");
	for (int i=0; i<mat.n; ++i) {
		for (int j=0; j<mat.m; ++j) {
			printf("%8.4Lf",mat.a[i][j]);
			if (j!=mat.m-1){
				printf("&");
			}
		}
		printf("\\\\\n");
	}
	printf("\\end{bmatrix}");
	putchar(sep);
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
//	printf("matget: ");
//	view(ans);
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
matrix matbind_c(matrix mleft, matrix mright){
	if (mleft.n!=mright.n) {
		printf("!!!Error binding.!!!\n");
		return mleft;
	}
	int le=mleft.m,ri=mright.m,nn=mleft.n;
	matrix ans=matmake(nn, le+ri);
	for (int j=0; j<le; ++j) {
		for (int i=0; i<nn; ++i) {
			ans.a[i][j]=mleft.a[i][j];
		}
	}
	for (int j=0; j<ri; ++j) {
		for (int i=0; i<nn; ++i) {
			ans.a[i][le+j]=mright.a[i][j];
		}
	}
//	printf("mbind: ");
//	view(ans);
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
matrix null_space(matrix mat){
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
}//finding independent rows
matrix colum_basis(matrix mat){
	return trans(row_basis(trans(mat)));
}
matrix transi_matrix(matrix tar, matrix ori){
	return inve(tar)*ori;
}
void view_by_colum(matrix mat){
	sep=',';
	for (int j=0; j<mat.m; ++j) {
		if (j==mat.m-1) {
			sep='\n';
		}
		view(matget(mat, 0, j, mat.n, 1));
	}
	prl();
}
void view_by_row(matrix mat){
	sep=',';
	for (int i=0; i<mat.n; ++i) {
		if (i==mat.n-1) {
			sep='\n';
		}
		view(matget(mat, i, 0, 1, mat.m));
	}
	prl();
}
LD det(matrix a){
	LD ans=1;
	int index[MATSIZE];
	for (int i=0; i<MATSIZE; ++i) {
		index[i]=i;
	}
	if (a.n!=a.m) {
		printf("The matrix is not a square matrix.");
	}
	int nn=a.n;
	LD u;
	for (int i=0; i<nn; ++i) {
		if (abs(a.a[i][i])<RES) {
			if (i==nn-1) {
				return 0;
			}
			for (int k=i+1; k<nn; ++k) {
				if (abs(a.a[k][i])>RES) {
					row_swap(a, i, k);
					swap(index[i], index[k]);
				}else if (k==nn){
					view(a);
					return 0;
				}
			}
		}
		u=a.a[i][i];
		ans*=u;
		for (int j=0; j<nn; ++j) {
			a.a[i][j]/=u;
		}
//		view(a);
		for (int k=0; k<nn; ++k) if (k!=i) {
			u=a.a[k][i];
			for (int j=0; j<nn; ++j) {
				a.a[k][j]-=u*a.a[i][j];
			}
		}
//		view(a);
	}
	view(a);
	return ans*(pow(-1, inve_pair(nn,index)));
}
LD cofactor(matrix mat, int i, int j){
	matrix a=matbind_c(matbind_r(matget(mat, 0, 0, i, j),
								 matget(mat, i+1, 0, mat.n-i-1, j)),
					   matbind_r(matget(mat, 0, j+1, i, mat.m-j-1),
								 matget(mat, i+1, j+1, mat.n-i-1, mat.m-j-1)));
	return det(a);
}
matrix adj(matrix mat){
	matrix ans=copysize(mat);
	for (int i=0; i<mat.n; ++i) {
		for (int j=0; j<mat.m; ++j) {
			ans.a[i][j]=cofactor(mat, i, j);
			if ((i+j)%2) ans.a[i][j]*=-1;
		}
	}
	return trans(ans);
}
void cramers(matrix mat,matrix b){
	if (b.m!=1||b.n!=mat.m) {
		printf("Incorrect vec b in Cramers.\n");
		exit(1);
	}
	matrix temp;
	LD deta=det(mat);
	if (deta==0) {
		printf("Singular matrix in Cramers.\n");
	}
	for (int i=0; i<mat.m; ++i) {
		temp=matbind_c(matget(mat, 0, 0, mat.n, i),b);
		temp=matbind_c(temp, matget(mat, 0, i+1, mat.n, mat.m-1-i));
		printf("detA_%2d=%8.4Lf, b_%2d=%8.4Lf,\n",i+1,det(temp),i+1,det(temp)/deta);
	}
}
void format_det(matrix a){
	printf("$$\ndet(");
	view(a);
	printf(")=%.4Lf\n$$\n",det(a));
}
matrix least_square(matrix a, matrix b){
//	robust needed here.
	return inve(trans(a)*a)*trans(a)*b;
}
matrix P(matrix a){
	return a*inve(trans(a)*a)*trans(a);
}

#endif /* fezzy_matrix_h */
