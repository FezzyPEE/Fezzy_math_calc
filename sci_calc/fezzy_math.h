//
//  fezzy_math.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_math_h
#define fezzy_math_h


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

#endif /* fezzy_math_h */
