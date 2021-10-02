//
//  fezzy_complex.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_complex_h
#define fezzy_complex_h

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

#endif /* fezzy_complex_h */
