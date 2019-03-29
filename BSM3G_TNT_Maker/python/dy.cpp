#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <math.h>
#include <fstream>

using namespace std;

int main(void)
{
 double TT=1213.288;
 double errTT=20.931;
 double tW=75.932;
 double errtW=3.620;
 double Other=1093.559;
 double errOther=23.773;
 double DY=66172.955;
 double errDY=365.029;
 double Data=76389.000;
 double errData=276.386;
 double D, errD, R, errR;

 D=Data-TT-tW-Other;
 errD=pow((errData*errData+errTT*errTT+errtW*errtW+errOther*errOther),0.5);
 R=D/DY;
 errR=R*pow(pow((errD/D),2)+pow((errDY/DY),2),0.5);
 cout<<"R="<<R<<"+-"<<errR<<"\n";

return 0;
}
