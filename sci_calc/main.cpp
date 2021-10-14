//
//  main.cpp
//  sci_calc
//
//  Created by 张子越 on 2020/9/24.
//  Copyright © 2020 张子越. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <unistd.h>
#include <vector>

using namespace std;

typedef long double LD;
typedef long long LL;

const LD Pi=3.1415926535897932384626433832795028841971693993751;
const LD g=9.8;
const LD R=8.3144621;
const LD atm=101325;
const LD RES=1e-8;

const int MPL=5;//maxpolylevel

const int Solve_Depth=8;

//---------------------------------------------IO
#include "fezzy_io.h"
//---------------------------------------------/IO


//---------------------------------------------math
#include "fezzy_math.h"
//---------------------------------------------/math

//---------------------------------------------/matrix
#include "fezzy_matrix.h"
//---------------------------------------------/matrix


//---------------------------------------------poly
#include "fezzy_poly.h"
//---------------------------------------------/poly


//---------------------------------------------complex 3D
//---------------------------------------------/complex 3D
#include "fezzy_complex.h"
//---------------------------------------------complex 2D
//----------------------------------------------/complex 2D


//----------------------------------------------system: centre of mass
#include "fezzy_phy_sys_com.h"
//----------------------------------------------/system: centre of mass


//----------------------------------------------angles
#include "fezzy_angles.h"
LD Theta;
LD sin(){
	return sin(Theta);
}
LD cos(){
	return cos(Theta);
}
//----------------------------------------------/angles

//----------------------------------------------PHY
LD phy_M,phy_Ekt;
#include "fezzy_phy.h"
//----------------------------------------------/PHY

//----------------------------------------------/STA
#include "fezzy_sta.h"
LD mpdfx[10],mpdfy[10];
LD mu_x,mu_y,sigmasqr_x,sigmasqr_y,cov_xy;
//----------------------------------------------/STA


int_ord_poly p,_dp;
complex2 ca,cb,cc,ct,f1,f2,f3;
LD lt1,lt2;
LD gx,gy;
LD phy_k,phy_h;
pt_2D p1,p2,p3;
sys2D_m_cen sys1;
LD ans;
int a[300][300];

#include "homework.h"

int main() {
//	LD h=0.5;
//	matrix a= matmake(3, 4, (LD[]){
//		1,-1,2,1,
//		2,1,1,8,
//		1,1,0,5
//	});
//	matrix b= matmake(2, 2, (LD[]){
//		1,0.5,
//		0.5,-4
//	});
//	view(a);
//	view(solve(a));
//	view(b);
//	matrix c=leftproduct(a+b, a-b);
//	view(c);
//	view(a*b*a);
//	view(b*a);
//	view(a*b);
	
	
	Q1();
	
	return 0;
}
