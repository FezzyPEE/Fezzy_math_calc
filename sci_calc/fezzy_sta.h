//
//  fezzy_sta.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_sta_h
#define fezzy_sta_h

void array_add_r(int **a,int i,int x,int m){
	for (int j = 0; j<m; ++j) {
		a[i][j]+=x;
	}
}
void array_add_c(int **a,int j,int x,int n){
	for (int i = 0; i<n; ++i) {
		a[i][j]+=x;
	}
}

const int SEQ_LEN=20;
struct seq {
	int len;
	LD a[SEQ_LEN];
};
seq seq_make(int len,LD a[]){
	seq ans;
	ans.len=len;
	for (int i=0; i<len; ++i) {
		ans.a[i]=a[i];
	}
	return ans;
}
LD seq_mu(seq s){
	LD ans=0;
	for (int i=0; i<s.len; ++i) {
		ans+=s.a[i];
	}
	ans/=s.len;
	return ans;
}
LD seq_sigma2(seq x,seq y){
	if (x.len!=y.len) {
		printf("Error in seq_sigma.\n");
		exit(1);
	}
	LD mu1=seq_mu(x),mu2=seq_mu(y);
	int n=x.len;
	LD ans=0;
	for (int i=0; i<n; ++i) {
		ans+=(x.a[i]-mu1)*(y.a[i]-mu2);
	}
	ans/=n;
	return ans;
}
LD seq_sigma2(seq s){
	return seq_sigma2(s, s);
}

#endif /* fezzy_sta_h */
