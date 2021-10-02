//
//  phy_sys_com.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef phy_sys_com_h
#define phy_sys_com_h

struct pt_2D {
	LD m;
	complex2 r;
};

struct sys2D_m_cen {
	LD M,ri0,ri;
	complex2 r,rM;
};

void init(sys2D_m_cen s){
	s.M=0;
	s.ri=s.ri0=0;
	s.r.len=0;
	init_p(s.r);
	s.rM=s.r;
}
void add_pt(sys2D_m_cen &s, pt_2D p){
	s.rM+=p.r*p.m;
	s.M+=p.m;
	s.r=s.rM/s.M;
}
void operator+=(sys2D_m_cen &s, pt_2D p){
	add_pt(s, p);
}

#endif /* phy_sys_com_h */
