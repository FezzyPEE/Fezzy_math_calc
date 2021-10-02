//
//  fezzy_io.h
//  sci_calc
//
//  Created by 张子越 on 2021/10/3.
//  Copyright © 2021 张子越. All rights reserved.
//

#ifndef fezzy_io_h
#define fezzy_io_h

void priLD(LD x){
	printf("%.4Lf\n",x);
}//defult .4
void priLD(LD x,int w){
	if (w>9) {
		printf("LD too long\n");
		return;
	}
	static char c[3];
	c[0]='0'+w;
	c[1]='\0';
	string wei=c;
	string s="%."+wei+"Lf\n";
	//printf("%s",s.c_str());
	printf(s.c_str(),x);
}

#endif /* fezzy_io_h */
