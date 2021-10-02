//
//  fezzy_sta.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_sta_h
#define fezzy_sta_h

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

#endif /* fezzy_sta_h */
