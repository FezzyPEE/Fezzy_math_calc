//
//  homework.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef homework_h
#define homework_h
void Q(){
	int id=0;
//	matrix a= matmake(3, 4, (LD[]){
//		1,-1,2,1,
//		2,1,1,8,
//		1,1,0,5
//	});
//	matrix b= matmake(2, 2, (LD[]){
//		1,0.5,
//		0.5,-4
//	});
	printf("## Q%d\n\n",id);
//	view(a);
//	view(b);
	
	printf("\n\n");
}


void Q1(){
	complex3 c,n,h;
	c = (complex3){0.579315,       0.000000,       0.000000};
	n = (complex3){2.766538,      -0.000000,      -0.000000};
	h = (complex3){-1.444769,      -0.000000,      -0.010000};
	LD aa=6*7/dis(c, n)+6*1/dis(c, h)+7*1/dis(h, n);
	n = (complex3){  0.203886,      -0.000000,      -0.000000};
	h = (complex3){ -1.686773,       0.000000,       0.000000};
	c = (complex3){  2.427750,       0.000000,      -0.010000};
	LD bb=6*7/dis(c, n)+6*1/dis(c, h)+7*1/dis(h, n);
	
	priLD(-aa*627.509474-58591.13);
	priLD(-bb*627.509474-58574.5);
}





#endif /* homework_h */
