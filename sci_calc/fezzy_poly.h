//
//  poly.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef poly_h
#define poly_h

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

#endif /* poly_h */
