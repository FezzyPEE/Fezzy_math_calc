//
//  fezzy_phy.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_phy_h
#define fezzy_phy_h

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

#endif /* fezzy_phy_h */
