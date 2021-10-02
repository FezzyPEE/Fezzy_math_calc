//
//  fezzy_angles.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_angles_h
#define fezzy_angles_h

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


#endif /* fezzy_angles_h */
