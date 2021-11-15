//
//  fezzy_lib.h
//  sci_calc
//
//  Created by 张子越 on 2021/11/15.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_lib_h
#define fezzy_lib_h

int inve_pair(int n,int a[]){
	int ans=0;
	for (int i=0; i<n; ++i) {
		for (int j=i; j<n; ++j) {
			if (a[i]>a[j]) {
				++ans;
			}
		}
	}
	return ans;
}


#endif /* fezzy_lib_h */
