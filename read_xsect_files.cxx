#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include "global.h"

using namespace std;

//This subroutine reads cross sections from files and fills out the following GLOBAL arrays:
//FOR RIPANI CROSS SECTIONS:
//W[17], Q2[3]
//S12_ARR[12][17],  S23_ARR[12][17],  THETA_ARR[6], ALPHA_ARR[6]
//SIGMA_ARR[6][3][17][12][12][6][6] <-> SIGMA_ARR[flag_sigma][q2_bin][w_binj][s23_bin][s12_bin][theta_bin][alpha_bin],
//where the flag_sigma corresponds to 0 - sigma_t, 1 - sigma_l, 2 - sigma_c2f, 3 - sigma_s2f, 4 - sigma_cf, 5 - sigma_sf
//----------------
//FOR GOLOVACH PHOTOPRODUCTION CROSS SECTIONS 
//W_ARR_GOL[30];
//S12_ARR_GOL[16][30];
//S23_ARR_GOL[16][30];
//THETA_ARR_GOL[14]; 
//ALPHA_ARR_GOL[14];
//SIGMA_ARR_GOL[30][16][16][14][14]; - correspond for sigma_t

void read_xsect_files(){
Float_t S12_val, S23_val, TH_val, ALP_val;
Float_t s12_min, s12_max,s23_min, s23_max, th_min, th_max, alp_min, alp_max;
Float_t dm12,dm23;
Float_t ds12_tmp,ds23_tmp, dalpha_tmp;
Float_t dalpha, ds12, ds23;
Float_t th_l, th_r, dtheta, dtheta_tmp;
string file_names[51];
string file_names_gol[30];
string file_names_fed[56];

string file_names_rip2[11];

Float_t Xsect_int,cross_sect,cross_sect_t,cross_sect_l,eps_l_rip2;
Short_t wbin, q2bin;
//Float_t th_l, th_r,dm12,dm23;
//files' names (Ripani_cr_sects)
file_names[0] = "ripani_cr_sect/rip_4diffsec_065_14125.dat";
file_names[1] = "ripani_cr_sect/rip_4diffsec_065_14375.dat";
file_names[2] = "ripani_cr_sect/rip_4diffsec_065_14625.dat";
file_names[3] = "ripani_cr_sect/rip_4diffsec_065_14875.dat";
file_names[4] = "ripani_cr_sect/rip_4diffsec_065_15125.dat";
file_names[5] = "ripani_cr_sect/rip_4diffsec_065_15375.dat";
file_names[6] = "ripani_cr_sect/rip_4diffsec_065_15625.dat";
file_names[7] = "ripani_cr_sect/rip_4diffsec_065_15875.dat";
file_names[8] = "ripani_cr_sect/rip_4diffsec_065_16125.dat";
file_names[9] = "ripani_cr_sect/rip_4diffsec_065_16375.dat";
file_names[10] = "ripani_cr_sect/rip_4diffsec_065_16625.dat";
file_names[11] = "ripani_cr_sect/rip_4diffsec_065_16875.dat";
file_names[12] = "ripani_cr_sect/rip_4diffsec_065_17125.dat";
file_names[13] = "ripani_cr_sect/rip_4diffsec_065_17375.dat";
file_names[14] = "ripani_cr_sect/rip_4diffsec_065_17625.dat";
file_names[15] = "ripani_cr_sect/rip_4diffsec_065_17875.dat";
file_names[16] = "ripani_cr_sect/rip_4diffsec_065_18125.dat";
file_names[17] = "ripani_cr_sect/rip_4diffsec_095_14125.dat";
file_names[18] = "ripani_cr_sect/rip_4diffsec_095_14375.dat";
file_names[19] = "ripani_cr_sect/rip_4diffsec_095_14625.dat";
file_names[20] = "ripani_cr_sect/rip_4diffsec_095_14875.dat";
file_names[21] = "ripani_cr_sect/rip_4diffsec_095_15125.dat";
file_names[22] = "ripani_cr_sect/rip_4diffsec_095_15375.dat";
file_names[23] = "ripani_cr_sect/rip_4diffsec_095_15625.dat";
file_names[24] = "ripani_cr_sect/rip_4diffsec_095_15875.dat";
file_names[25] = "ripani_cr_sect/rip_4diffsec_095_16125.dat";
file_names[26] = "ripani_cr_sect/rip_4diffsec_095_16375.dat";
file_names[27] = "ripani_cr_sect/rip_4diffsec_095_16625.dat";
file_names[28] = "ripani_cr_sect/rip_4diffsec_095_16875.dat";
file_names[29] = "ripani_cr_sect/rip_4diffsec_095_17125.dat";
file_names[30] = "ripani_cr_sect/rip_4diffsec_095_17375.dat";
file_names[31] = "ripani_cr_sect/rip_4diffsec_095_17625.dat";
file_names[32] = "ripani_cr_sect/rip_4diffsec_095_17875.dat";
file_names[33] = "ripani_cr_sect/rip_4diffsec_095_18125.dat";
file_names[34] = "ripani_cr_sect/rip_4diffsec_130_14125.dat";
file_names[35] = "ripani_cr_sect/rip_4diffsec_130_14375.dat";
file_names[36] = "ripani_cr_sect/rip_4diffsec_130_14625.dat";
file_names[37] = "ripani_cr_sect/rip_4diffsec_130_14875.dat";
file_names[38] = "ripani_cr_sect/rip_4diffsec_130_15125.dat";
file_names[39] = "ripani_cr_sect/rip_4diffsec_130_15375.dat";
file_names[40] = "ripani_cr_sect/rip_4diffsec_130_15625.dat";
file_names[41] = "ripani_cr_sect/rip_4diffsec_130_15875.dat";
file_names[42] = "ripani_cr_sect/rip_4diffsec_130_16125.dat";
file_names[43] = "ripani_cr_sect/rip_4diffsec_130_16375.dat";
file_names[44] = "ripani_cr_sect/rip_4diffsec_130_16625.dat";
file_names[45] = "ripani_cr_sect/rip_4diffsec_130_16875.dat";
file_names[46] = "ripani_cr_sect/rip_4diffsec_130_17125.dat";
file_names[47] = "ripani_cr_sect/rip_4diffsec_130_17375.dat";
file_names[48] = "ripani_cr_sect/rip_4diffsec_130_17625.dat";
file_names[49] = "ripani_cr_sect/rip_4diffsec_130_17875.dat";
file_names[50] = "ripani_cr_sect/rip_4diffsec_130_18125.dat";
//files' names (Golovach_cr_sects)

file_names_gol[0] = "golovach_cr_sect/5diffsec161.dat";
file_names_gol[1] = "golovach_cr_sect/5diffsec164.dat";
file_names_gol[2] = "golovach_cr_sect/5diffsec166.dat";
file_names_gol[3] = "golovach_cr_sect/5diffsec169.dat";
file_names_gol[4] = "golovach_cr_sect/5diffsec171.dat";
file_names_gol[5] = "golovach_cr_sect/5diffsec174.dat";
file_names_gol[6] = "golovach_cr_sect/5diffsec176.dat";
file_names_gol[7] = "golovach_cr_sect/5diffsec179.dat";
file_names_gol[8] = "golovach_cr_sect/5diffsec181.dat";
file_names_gol[9] = "golovach_cr_sect/5diffsec184.dat";
file_names_gol[10] = "golovach_cr_sect/5diffsec186.dat";
file_names_gol[11] = "golovach_cr_sect/5diffsec189.dat";
file_names_gol[12] = "golovach_cr_sect/5diffsec191.dat";
file_names_gol[13] = "golovach_cr_sect/5diffsec194.dat";
file_names_gol[14] = "golovach_cr_sect/5diffsec196.dat";
file_names_gol[15] = "golovach_cr_sect/5diffsec199.dat";
file_names_gol[16] = "golovach_cr_sect/5diffsec201.dat";
file_names_gol[17] = "golovach_cr_sect/5diffsec204.dat";
file_names_gol[18] = "golovach_cr_sect/5diffsec206.dat";
file_names_gol[19] = "golovach_cr_sect/5diffsec209.dat";
file_names_gol[20] = "golovach_cr_sect/5diffsec211.dat";
file_names_gol[21] = "golovach_cr_sect/5diffsec214.dat";
file_names_gol[22] = "golovach_cr_sect/5diffsec219.dat";
file_names_gol[23] = "golovach_cr_sect/5diffsec224.dat";
file_names_gol[24] = "golovach_cr_sect/5diffsec228.dat";
file_names_gol[25] = "golovach_cr_sect/5diffsec234.dat";
file_names_gol[26] = "golovach_cr_sect/5diffsec239.dat";
file_names_gol[27] = "golovach_cr_sect/5diffsec244.dat";
file_names_gol[28] = "golovach_cr_sect/5diffsec249.dat";
file_names_gol[29] = "golovach_cr_sect/5diffsec254.dat";

//files' names (Fedotov_cr_sects)
file_names_fed[0] =  "fedotov_cr_sect/fedotov_4diffsec_0225_15125.dat";
file_names_fed[1] =  "fedotov_cr_sect/fedotov_4diffsec_0225_15375.dat";
file_names_fed[2] =  "fedotov_cr_sect/fedotov_4diffsec_0225_15625.dat";
file_names_fed[3] =  "fedotov_cr_sect/fedotov_4diffsec_0225_15875.dat";
file_names_fed[4] =  "fedotov_cr_sect/fedotov_4diffsec_0275_13125.dat";
file_names_fed[5] =  "fedotov_cr_sect/fedotov_4diffsec_0275_13375.dat";
file_names_fed[6] =  "fedotov_cr_sect/fedotov_4diffsec_0275_13625.dat";
file_names_fed[7] =  "fedotov_cr_sect/fedotov_4diffsec_0275_13875.dat";
file_names_fed[8] =  "fedotov_cr_sect/fedotov_4diffsec_0275_14125.dat";
file_names_fed[9] =  "fedotov_cr_sect/fedotov_4diffsec_0275_14375.dat";
file_names_fed[10] =  "fedotov_cr_sect/fedotov_4diffsec_0275_14625.dat";
file_names_fed[11] =  "fedotov_cr_sect/fedotov_4diffsec_0275_14875.dat";
file_names_fed[12] =  "fedotov_cr_sect/fedotov_4diffsec_0275_15125.dat";
file_names_fed[13] =  "fedotov_cr_sect/fedotov_4diffsec_0275_15375.dat";
file_names_fed[14] =  "fedotov_cr_sect/fedotov_4diffsec_0275_15625.dat";
file_names_fed[15] =  "fedotov_cr_sect/fedotov_4diffsec_0325_13125.dat";
file_names_fed[16] =  "fedotov_cr_sect/fedotov_4diffsec_0325_13375.dat";
file_names_fed[17] =  "fedotov_cr_sect/fedotov_4diffsec_0325_13625.dat";
file_names_fed[18] =  "fedotov_cr_sect/fedotov_4diffsec_0325_13875.dat";
file_names_fed[19] =  "fedotov_cr_sect/fedotov_4diffsec_0325_14125.dat";
file_names_fed[20] =  "fedotov_cr_sect/fedotov_4diffsec_0325_14375.dat";
file_names_fed[21] =  "fedotov_cr_sect/fedotov_4diffsec_0325_14625.dat";
file_names_fed[22] =  "fedotov_cr_sect/fedotov_4diffsec_0325_14875.dat";
file_names_fed[23] =  "fedotov_cr_sect/fedotov_4diffsec_0325_15125.dat";
file_names_fed[24] =  "fedotov_cr_sect/fedotov_4diffsec_0325_15375.dat";
file_names_fed[25] =  "fedotov_cr_sect/fedotov_4diffsec_0425_13125.dat";
file_names_fed[26] =  "fedotov_cr_sect/fedotov_4diffsec_0425_13375.dat";
file_names_fed[27] =  "fedotov_cr_sect/fedotov_4diffsec_0425_13625.dat";
file_names_fed[28] =  "fedotov_cr_sect/fedotov_4diffsec_0425_13875.dat";
file_names_fed[29] =  "fedotov_cr_sect/fedotov_4diffsec_0425_14125.dat";
file_names_fed[30] =  "fedotov_cr_sect/fedotov_4diffsec_0425_14375.dat";
file_names_fed[31] =  "fedotov_cr_sect/fedotov_4diffsec_0425_14625.dat";
file_names_fed[32] =  "fedotov_cr_sect/fedotov_4diffsec_0425_14875.dat";
file_names_fed[33] =  "fedotov_cr_sect/fedotov_4diffsec_0425_15125.dat";
file_names_fed[34] =  "fedotov_cr_sect/fedotov_4diffsec_0475_13125.dat";
file_names_fed[35] =  "fedotov_cr_sect/fedotov_4diffsec_0475_13375.dat";
file_names_fed[36] =  "fedotov_cr_sect/fedotov_4diffsec_0475_13625.dat";
file_names_fed[37] =  "fedotov_cr_sect/fedotov_4diffsec_0475_13875.dat";
file_names_fed[38] =  "fedotov_cr_sect/fedotov_4diffsec_0475_14125.dat";
file_names_fed[39] =  "fedotov_cr_sect/fedotov_4diffsec_0475_14375.dat";
file_names_fed[40] =  "fedotov_cr_sect/fedotov_4diffsec_0475_14625.dat";
file_names_fed[41] =  "fedotov_cr_sect/fedotov_4diffsec_0475_14875.dat";
file_names_fed[42] =  "fedotov_cr_sect/fedotov_4diffsec_0525_13125.dat";
file_names_fed[43] =  "fedotov_cr_sect/fedotov_4diffsec_0525_13375.dat";
file_names_fed[44] =  "fedotov_cr_sect/fedotov_4diffsec_0525_13625.dat";
file_names_fed[45] =  "fedotov_cr_sect/fedotov_4diffsec_0525_13875.dat";
file_names_fed[46] =  "fedotov_cr_sect/fedotov_4diffsec_0525_14125.dat";
file_names_fed[47] =  "fedotov_cr_sect/fedotov_4diffsec_0525_14375.dat";
file_names_fed[48] =  "fedotov_cr_sect/fedotov_4diffsec_0525_14625.dat";
file_names_fed[49] =  "fedotov_cr_sect/fedotov_4diffsec_0525_14875.dat";
file_names_fed[50] =  "fedotov_cr_sect/fedotov_4diffsec_0575_13125.dat";
file_names_fed[51] =  "fedotov_cr_sect/fedotov_4diffsec_0575_13375.dat";
file_names_fed[52] =  "fedotov_cr_sect/fedotov_4diffsec_0575_13625.dat";
file_names_fed[53] =  "fedotov_cr_sect/fedotov_4diffsec_0575_13875.dat";
file_names_fed[54] =  "fedotov_cr_sect/fedotov_4diffsec_0575_14125.dat";
file_names_fed[55] =  "fedotov_cr_sect/fedotov_4diffsec_0575_14375.dat";

file_names_rip2[0] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_18375.dat";
file_names_rip2[1] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_18625.dat";
file_names_rip2[2] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_18875.dat";
file_names_rip2[3] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_19125.dat";
file_names_rip2[4] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_19375.dat";
file_names_rip2[5] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_19625.dat";
file_names_rip2[6] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_19875.dat";
file_names_rip2[7] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_20125.dat";
file_names_rip2[8] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_20375.dat";
file_names_rip2[9] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_20625.dat";
file_names_rip2[10] =  "rip_q2_130_w_18_21/wgt18_4diffsec_130_20875.dat";

//Define theta, alpha and Q2 arrays for Ripani

THETA_ARR[0] = 0.;
THETA_ARR[1] = 0.6343185;
THETA_ARR[2] = 1.258637;
THETA_ARR[3] = 1.882956;
THETA_ARR[4] = 2.507274;
THETA_ARR[5] = 3.141593;


ALPHA_ARR[0] = 0.;
ALPHA_ARR[1] = 1.262637;
ALPHA_ARR[2] = 2.515274;
ALPHA_ARR[3] = 3.767911;
ALPHA_ARR[4] = 5.020548;
ALPHA_ARR[5] = 6.283186;

Q2_ARR[0] = 0.65;
Q2_ARR[1] = 0.95;
Q2_ARR[2] = 1.30;


//Float_t dtheta = (THETA_ARR[5] - THETA_ARR[0])/5.;



//loop over files for Ripani
for (Short_t i=0; i<=50; i++) {

wbin = i - Int_t(i/17)*17;
q2bin = Int_t(i/17);
//Define W array
W_ARR[i - Int_t(i/17)*17] = 1.4125 + 0.025*(i - Int_t(i/17)*17);

//cout<<i <<" "<< Int_t(i/17) << " "<<i - Int_t(i/17)*17<< "\n";
//cout << "Reading RIPANI diff cross sections for Q^2 = "<< Q2_ARR[Int_t(i/17)]<< " GeV^2, W = "<< W_ARR[i - Int_t(i/17)*17] << " GeV \n";

string dummy,xsect;

string file=file_names[i];
ifstream input(file.c_str());
if(input.is_open()){

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
S12_ARR[is12-1][wbin] = atof(xsect.c_str());
getline(input,xsect);
S23_ARR[is23-1][wbin] = atof(xsect.c_str());
getline(input,dummy);
getline(input,dummy);

getline(input,xsect);
//sigma_t
SIGMA_ARR[0][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_l
SIGMA_ARR[1][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_c2f
SIGMA_ARR[2][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_s2f
SIGMA_ARR[3][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_cf
SIGMA_ARR[4][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_sf
SIGMA_ARR[5][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
EPSILON_L_RIPANI[q2bin][wbin] = atof(xsect.c_str());
getline(input,dummy);



};
};
};
};
};
input.close();

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {



for (Short_t j=0;j<6;j++){


//SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR[11][wbin]-S12_ARR[0][wbin])*(S23_ARR[11][wbin]-S23_ARR[0][wbin]);

SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*2.*M_PI*2.*M_PI;


if ((q2bin==2)&&(wbin==15)) SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.02;
if ((q2bin==2)&&(wbin==16)) SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.95;
};

};
};
};
};

};//end loop-ripani  - i
///---------------------------------------------------------------------------------------
//Define theta, alpha arrays for Golovach

THETA_ARR_GOL[0] = 0.;
THETA_ARR_GOL[1] = 0.2501225;
THETA_ARR_GOL[2] =0.4902450;
THETA_ARR_GOL[3] =0.7303675;
THETA_ARR_GOL[4] =0.9704900;
THETA_ARR_GOL[5] =1.210613;
THETA_ARR_GOL[6] =1.450735;
THETA_ARR_GOL[7] =1.690858;
THETA_ARR_GOL[8] =1.930980;
THETA_ARR_GOL[9] =2.171103;
THETA_ARR_GOL[10] =2.411225;
THETA_ARR_GOL[11] =2.651348;
THETA_ARR_GOL[12] =2.891470;
THETA_ARR_GOL[13] =3.141593;

ALPHA_ARR_GOL[0] = 0.;
ALPHA_ARR_GOL[1] = 0.4917835;
ALPHA_ARR_GOL[2] = 0.9735669;
ALPHA_ARR_GOL[3] = 1.455350;
ALPHA_ARR_GOL[4] = 1.937134;
ALPHA_ARR_GOL[5] = 2.418917;
ALPHA_ARR_GOL[6] = 2.900701;
ALPHA_ARR_GOL[7] = 3.382484;
ALPHA_ARR_GOL[8] = 3.864268;
ALPHA_ARR_GOL[9] = 4.346051;
ALPHA_ARR_GOL[10] = 4.827835;
ALPHA_ARR_GOL[11] = 5.309618;
ALPHA_ARR_GOL[12] = 5.791402;
ALPHA_ARR_GOL[13] = 6.283185;


//loop over files for Golovach
for (Short_t i=0; i<=29; i++) {


//Define W array
if ((i>=0)&&(i<=21)) W_ARR_GOL[i] = 1.6125 + 0.025*i;

if ((i>=22)&&(i<=29)) W_ARR_GOL[i] = 2.1875 + 0.05*(i-22);

//cout<<i <<" "<< Int_t(i/17) << " "<<i - Int_t(i/17)*17<< "\n";
//cout << "Reading GOLOVACH diff cross sections for W = "<< W_ARR_GOL[i] << " GeV \n";

string dummy,xsect;

string file=file_names_gol[i];
ifstream input(file.c_str());
if(input.is_open()){

for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
//S12_val = atof(xsect.c_str());
//cout << "qqq" << 
S12_ARR_GOL[is12-1][i] = atof(xsect.c_str());
getline(input,xsect);
//S23_val = atof(xsect.c_str());
S23_ARR_GOL[is23-1][i] = atof(xsect.c_str());
getline(input,dummy);
//TH_val = atof(dummy.c_str());
getline(input,dummy);
//ALP_val = atof(dummy.c_str());
//cout << i-51 << " "<<W_ARR_GOL[i-51] << " "<< S23_ARR_GOL[is23-1][i-51]  << "\n";
getline(input,xsect);
//sigma_t

SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
getline(input,dummy);


/*
if ((is23==1)&&(itheta==1)&&(ialpha==1)&&(is12==1)) s12_min = S12_val;
if ((is23==1)&&(itheta==1)&&(ialpha==1)&&(is12==16)) s12_max = S12_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) s23_min = S23_val;
if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==16)) s23_max = S23_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) th_min = TH_val;
if ((is12==1)&&(itheta==14)&&(ialpha==1)&&(is23==1)) th_max = TH_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) alp_min = ALP_val;
if ((is12==1)&&(itheta==1)&&(ialpha==14)&&(is23==1)) alp_max = ALP_val;
*/
};
};
};
};
};
input.close();


//cout << "qqq" << s12_min << "\n";


for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {



//SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_GOL[15][i]-S12_ARR_GOL[0][i])*(S23_ARR_GOL[15][i]-S23_ARR_GOL[0][i]);

SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*2.*M_PI*2.*M_PI;
//SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*2.*M_PI;

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.932;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.928;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.98;

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.92*1.006*1.000051212;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.98*1.00773713;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.91*0.993;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.97*1.0129949;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.986;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.03;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.03*1.019;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.1*1.05*0.9696;



//----------------------------------------
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*0.98;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*0.95;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.05*0.95;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.04*0.95;

if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.06*0.95;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.04*0.95;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*0.95;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.13*0.95;


if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.01*0.95;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.06*0.95;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.99*0.95;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.0105*0.95;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.02*0.95;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.08*0.95;
//---------------------

if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.03*0.95;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.03*1.01;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.16*1.02;

if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.02*0.95;

if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.01*1.04*0.95;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.1*1.0*0.954;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*1.07*0.95;

if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*0.95;
};
};
};
};

for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.8;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.91;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.94;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.1;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.3;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.3;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.15;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.15;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.6;


if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.91;
//if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.1;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.3;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.35;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.15;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.15;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.2;


if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.78;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.1;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.3;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.3;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.2;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.2;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.2;


if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.62;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.97;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1]*1.05;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.2;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.35;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.3;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.;


if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*1.25;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.9;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.25;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.35;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.55;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.5;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.5;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.5;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*1.3;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*1.1;
//if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.9;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.5;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.4;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.6;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.7;

if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.85;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.9;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.8;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.05;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.1;

if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.87;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.83;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.9;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.85;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.85;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*0.9;
//-----------------------------------------------
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.75;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.81;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.9;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.7;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.03;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.8;



if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.52;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.58;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.64;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1]*0.88;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.95;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.8;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.03;

if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.85;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.85;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*0.9;


if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.515;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.54;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.69;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1]*0.92;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.8;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.07;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.9;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.9;




if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.43;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.6;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.8;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1]*0.97;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.76;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.71;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.07;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.84;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.78;




if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.43;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.57;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.83;

if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1]*0.97;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.73;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.05;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.96;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.96;






if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.43;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.6;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.8;

if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.9;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.75;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.07;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.94;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.96;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.68;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.58;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.74;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.85;








if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.48;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.4;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.65;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.7;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*0.9;




if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.25;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.44;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.55;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1]*0.85;
//if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.78;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.2;
//if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.96;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.96;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.2;



if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.83;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.53;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.87;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.6;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.1;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.47;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.83;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.7;




if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.233;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.48;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;


if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.21;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.48;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.78;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;


if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.2;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.3;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.65;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.3;


//-------------------------------------------------


if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.28;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.55;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;

if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.75;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.2;



if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.4;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.58;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;

if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.8;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.5;


if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.25;



if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.5;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.75;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.8;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.6;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.2;


if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.76;
//if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.97;
//if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.8;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.4;




if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.95;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.9;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.4;


if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*1.07;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.95;
//if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.9;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.3;


if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.81;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;

if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.6;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.3;



if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*1.23;
//if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.93;


if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.27;

};
};
};


for (Int_t is23 = 1; is23 <=16; is23++) {

for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {

SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.94;


if (i==0) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.9;
if (i==0) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.95;
if (i==0) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.95;
if (i==0) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.75;

if (i==1) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.9;
if (i==1) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.92;
if (i==1) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.92;
if (i==1) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.9;


if (i==2) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.8;
if (i==2) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;

if (i==3) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.75;
if (i==3) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.1;
if (i==3) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.15;
if (i==3) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.15;
if (i==3) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==3) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.05;



if (i==4) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.75;
if (i==4) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;
if (i==4) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.8;
if (i==4) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.05;
if (i==4) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.1;


if (i==5) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.75;
if (i==5) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;
if (i==5) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.8;
if (i==5) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.94;
if (i==5) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.98;
if (i==5) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.85;

if (i==6) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.55;
if (i==6) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.85;
if (i==6) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.8;
if (i==6) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.1;
if (i==6) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.15;


if (i==7) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.55;
if (i==7) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.85;
if (i==7) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.78;
if (i==7) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*0.95;
if (i==7) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.8;


//if ((i>=8)&&(i<=12)) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.95;

//--------------------------------------------------------
if ((i>=8)&&(i<=12)) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.95;

if (i==8) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.25;
if (i==8) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==8) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.1;
if (i==8) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.3;
if (i==8) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.3;
if (i==8) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.34;

if (i==8) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.34;
if (i==8) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.26;
if (i==8) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*1.1;
if (i==8) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==8) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.5;

if (i==9) SIGMA_ARR_GOL[i][is23-1][0][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][0][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.5;
if (i==9) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.56;
if (i==9) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;

if (i==9) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.13;
if (i==9) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==9) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.55;
if (i==9) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.55;
if (i==9) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.54;

if (i==9) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.4;
if (i==9) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.12;
if (i==9) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.52;
if (i==9) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.37;


if (i==10) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==10) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.55;
if (i==10) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==10) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.15;
if (i==10) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.27;
if (i==10) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.27;
if (i==10) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.29;
if (i==10) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;

if (i==10) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.95;
if (i==10) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.95;
if (i==10) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.95;
if (i==10) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==10) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==11) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==11) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.68;
if (i==11) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==11) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.14;
if (i==11) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.32;
if (i==11) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.37;
if (i==11) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.36;
if (i==11) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.17;



if (i==11) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.9;
if (i==11) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==11) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==12) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==12) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==12) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==12) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==12) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.39;
if (i==12) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.42;
if (i==12) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==12) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.17;




if (i==12) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==12) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


//if (i==13) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.93;

if (i==13) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==13) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==13) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==13) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==13) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.43;
if (i==13) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.47;
if (i==13) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.42;
if (i==13) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.13;

if (i==13) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==13) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.6;
if (i==13) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.5;
if (i==13) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.45;
if (i==13) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.27;


if (i==14) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.7;
if (i==14) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][0][ialpha-1]*0.93;
if (i==14) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][0][ialpha-1]*0.97;


if (i==14) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==14) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==14) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==14) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==14) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.39;
if (i==14) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.42;
if (i==14) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==14) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.14;

if (i==14) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.05;
//if (i==14) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==14) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.6;
if (i==14) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.9;
if (i==14) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.3;
if (i==14) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.2;



if (i==15) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.65;
if (i==15) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.75;


if (i==15) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==15) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==15) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==15) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==15) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.39;
if (i==15) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.42;
if (i==15) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==15) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.17;

if (i==15) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.05;
//if (i==15) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.8;
if (i==15) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.6;
if (i==15) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==15) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.2;






if (i==16) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.8;

if (i==16) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==16) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==16) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==16) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.23;
if (i==16) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==16) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.52;
if (i==16) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==16) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==16) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.2;
//if (i==15) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==16) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.8;
if (i==16) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.6;
if (i==16) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==16) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.7;




if (i==17) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.83;
if (i==17) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][0][ialpha-1]*0.83;


if (i==17) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==17) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==17) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==17) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.23;
if (i==17) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==17) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.49;
if (i==17) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==17) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==17) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.2;
//if (i==15) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.8;
if (i==17) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.6;
if (i==17) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==17) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.5;



if (i==18) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.85;
if (i==18) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][0][ialpha-1]*0.85;

if (i==18) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==18) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.7;
if (i==18) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.7;
if (i==18) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.7;



if (i==18) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==18) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==18) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==18) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.23;
if (i==18) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==18) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.5;
if (i==18) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==18) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==18) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.2;
if (i==18) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.9;
if (i==18) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.6;
if (i==18) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.25;
if (i==18) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.9;



if (i==19) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==19) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.5;
if (i==19) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.8;

if (i==19) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.23;
if (i==19) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==19) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==19) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==19) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==19) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.72;
if (i==19) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.55;
if (i==19) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.4;
if (i==19) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.37;


if (i==20) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==20) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.4;
if (i==20) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;

if (i==20) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.35;
if (i==20) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==20) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==20) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==20) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==20) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.72;
if (i==20) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.55;
if (i==20) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.4;
if (i==20) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.37;







if (i==21) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.85;
if (i==21) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][0][ialpha-1]*0.85;
if (i==21) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][0][ialpha-1]*0.95;

if ((i>=14)&&(i<=21)) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.9;


if (i==21) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.85;
if (i==21) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.9;
if (i==21) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==21) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.85;


if (i==21) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==21) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.4;
if (i==21) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;

if (i==21) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.35;
if (i==21) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.4;
if (i==21) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==21) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==21) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==21) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.95;
//if (i==21) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.72;
if (i==21) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*1.25;
if (i==21) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.3;
if (i==21) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*1.1;


//---------------------------------------------------------------------------


if (i==22) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==22) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.35;

if (i==22) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.45;
if (i==22) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==22) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.4;
if (i==22) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.3;
if (i==22) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.4;
if (i==22) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.2;
if (i==22) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.1;
if (i==22) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==22) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.15;
if (i==22) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.1;





if (i==23) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==23) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==23) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.45;
if (i==23) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;

if (i==23) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.6;
if (i==23) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.8;


if (i==23) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.6;
if (i==23) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==23) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.2;
if (i==23) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==23) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.15;
if (i==23) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.1;



if (i==24) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==24) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==24) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==24) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==24) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.7;

if (i==24) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.15;
if (i==24) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.48;
if (i==24) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.8;


if (i==24) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.7;
if (i==24) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==24) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.2;
if (i==24) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==24) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.15;
if (i==24) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.1;




if (i==25) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==25) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==25) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==25) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.25;


if (i==25) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==25) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==25) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.1;
if (i==25) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.7;
if (i==25) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.9;


if (i==25) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.7;
if (i==25) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.3;
if (i==25) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.1;
if (i==25) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*0.82;
if (i==25) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.05;
if (i==25) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.85;
if (i==25) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.65;




if (i==26) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==26) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==26) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==26) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.15;

if (i==26) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==26) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==26) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.1;
if (i==26) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.9;
if (i==26) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.8;


if (i==26) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.5;
if (i==26) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.32;
if (i==26) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==26) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.05;
if (i==26) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.1;
if (i==26) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.16;
if (i==26) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.98;

if (i==26) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.15;



if (i==27) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==27) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==27) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==27) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;

if (i==27) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==27) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==27) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.67;
if (i==27) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*2.;
if (i==27) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.75;


if (i==27) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.5;
if (i==27) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.42;
if (i==27) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==27) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.1;
if (i==27) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.1;
if (i==27) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.12;
if (i==27) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.98;

if (i==27) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.08;


if (i==28) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.75;

if (i==28) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==28) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==28) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==28) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;

if (i==28) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==28) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==28) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;
if (i==28) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.7;
if (i==28) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*2.05;


if (i==28) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.75;
if (i==28) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.73;
if (i==28) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.49;
if (i==28) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.73;
if (i==28) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.6;
if (i==28) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.1;
if (i==28) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*1.45;

if (i==28) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.5;



//if (i==27) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==29) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.5;
if (i==29) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==29) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==29) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==29) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==29) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
//if (i==29) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;





if (i==29) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==29) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==29) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.;
if (i==29) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.7;
if (i==29) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.9;


if (i==29) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.65;
if (i==29) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.65;
if (i==29) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.46;
if (i==29) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.54;
if (i==29) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*3.02;
if (i==29) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.85;
if (i==29) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*1.3;

if (i==29) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*2.;




};
};
};




for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {


if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.1;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.1;


if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.2;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.2;

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.1;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.1;


if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.1;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.11;

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.05;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.05;


if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.1;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.1;


if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.15;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.15;

if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;


if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.3;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.3;

if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.2;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.2;



if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.1;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.1;

if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;


if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.3;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.3;

if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.2;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.2;


if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.9;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.9;

if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.9;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.9;

if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;

if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.1;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.1;




if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.9;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.9;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.8;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.8;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.75;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.75;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.7;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.7;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.9;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.9;

if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.85;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.85;

if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;

if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.75;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.75;

//--------------------------------------------
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.45;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.45;

if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.35;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.35;




if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.96;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.78;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;

if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;


if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.78;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.96;


if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.65;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.65;

if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.45;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.45;

if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.03;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.78;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.75;


if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;


if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.75;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.79;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.03;



if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;

if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;

if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.9;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.9;


if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.85;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.85;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.92;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.89;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.89;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.92;


if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;


if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.92;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.92;


if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;


if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;



if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.05;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.05;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.95;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.95;



if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;


if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;

if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.9;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.9;




//if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.95;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.92;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.03;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.92;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.92;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.03;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.92;
//if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.95;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;

if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.1;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.1;

if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.9;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.9;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.9;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.9;



if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.03;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.03;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.08;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.08;


if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.85;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.85;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.95;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.95;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;


if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.65;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.65;


if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.65;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.65;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.75;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.75;






if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.08;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.25;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.87;


if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.08;

if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.85;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.85;

if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.25;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.87;


if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;


if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.75;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.75;

if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;



if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.08;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.08;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.2;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.2;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.9;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.9;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.7;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.7;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.7;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.7;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.9;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.9;


if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.9;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.9;



if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.2;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.2;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;





if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.1;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.1;



if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.5;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.5;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;





if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.1;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.1;




if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.5;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.5;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;



if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.1;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.1;


if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.5;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.5;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.7;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.7;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


//----------------------------------------------------------

if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.7;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.7;

if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;



if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.25;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.3;

if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.4;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.4;

if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;



if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.35;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.4;

if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.5;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.5;

if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;




if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.35;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.35;

if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.5;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.5;

if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;




if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.35;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.35;

if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.5;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.5;

if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;



if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;


if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.55;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.55;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.6;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.6;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.65;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.65;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.75;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.75;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.9;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.9;


if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.35;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.35;


if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.5;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.5;

if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.3;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.3;

if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.15;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.15;


if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.15;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.15;



if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.55;
if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.55;

if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;

if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.9;
if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.9;



};
};
};


for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {



if (i==0) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.1;
if (i==0) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.2;
if (i==0) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.7;
if (i==0) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.8;
if (i==0) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.05;

if (i==1) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.9;

if (i==1) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.2;
if (i==1) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.1;
if (i==1) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.85;
if (i==1) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.7;
if (i==1) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==1) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.05;
if (i==1) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.2;
if (i==1) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.2;



if (i==2) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==2) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==2) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*0.9;

if (i==2) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.8;
if (i==2) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.9;
if (i==2) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.1;


if (i==3) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.8;
if (i==3) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==3) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.8;
if (i==3) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.93;
if (i==3) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.05;
if (i==3) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.1;

if (i==4) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==4) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==4) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.75;
if (i==4) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.95;


if (i==5) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.95;
if (i==5) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.05;
if (i==5) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.28;
if (i==5) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.06;
if (i==5) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.04;
if (i==5) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*0.96;
if (i==5) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.17;


if (i==6) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.85;
if (i==6) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.25;
if (i==6) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.25;



if (i==7) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.8;
if (i==7) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.15;
if (i==7) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.25;
if (i==7) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.25;



//----------------------------------------------

if (i==8) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.9;
if (i==8) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.7;

if (i==8) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;


if (i==8) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==8) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.17;
if (i==8) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.3;
if (i==8) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.35;
if (i==8) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.12;


if (i==9) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==9) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.97;
if (i==9) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.71;
if (i==9) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.84;

if (i==9) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.15;
if (i==9) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.45;
if (i==9) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.54;
if (i==9) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.57;
if (i==9) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.38;
if (i==9) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.2;
if (i==9) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.08;
if (i==9) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.8;




if (i==10) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.8;
if (i==10) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.89;
if (i==10) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.9;

if (i==10) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.13;
if (i==10) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==10) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.42;
if (i==10) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.45;
if (i==10) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.36;
if (i==10) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.27;
if (i==10) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.1;

//if (i==10) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.9;
if (i==10) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*0.9;

if (i==11) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.85;
if (i==11) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.92;
if (i==11) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.8;

if (i==11) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==11) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.13;
if (i==11) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==11) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.42;
if (i==11) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.48;
if (i==11) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.44;
if (i==11) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.3;
if (i==11) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.15;


if (i==12) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.4;
if (i==12) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.63;
if (i==12) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.85;
if (i==12) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.75;


//if (i==12) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.1;
if (i==12) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.18;
if (i==12) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.4;
if (i==12) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.45;
if (i==12) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.4;
if (i==12) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.3;
if (i==12) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.15;



if (i==13) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.6;
if (i==13) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.94;
if (i==13) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.83;
if (i==13) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.83;

if (i==13) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*0.99;
if (i==13) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.35;
if (i==13) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.6;
if (i==13) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.8;
if (i==13) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.8;
if (i==13) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.6;
if (i==13) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.4;
if (i==13) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.2;




if (i==14) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.3;
if (i==14) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.91;
if (i==14) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.91;
if (i==14) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.83;

if (i==14) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==14) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.45;
if (i==14) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.9;
if (i==14) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*2.17;
if (i==14) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*2.2;
if (i==14) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.6;
if (i==14) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.35;
if (i==14) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.2;


if (i==15) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*0.8;



if (i==15) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.3;
if (i==15) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.74;
if (i==15) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.97;

if (i==15) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.24;
if (i==15) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.63;
if (i==15) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*2.02;
if (i==15) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*2.1;
if (i==15) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*2.0;
if (i==15) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.89;
if (i==15) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.95;
if (i==15) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.6;



if (i==16) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.7;
if (i==16) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==16) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.09;
if (i==16) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.15;

if (i==16) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.37;
if (i==16) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.53;
if (i==16) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.62;
if (i==16) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.7;
if (i==16) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.72;
if (i==16) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.63;
if (i==16) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.5;
if (i==16) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.4;






if (i==17) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.5;
if (i==17) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.5;
if (i==17) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.63;
if (i==17) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.81;
if (i==17) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*0.8;
if (i==17) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;



if (i==17) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*1.2;
if (i==17) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.38;
if (i==17) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.12;
//if (i==17) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.96;


if (i==17) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.22;
if (i==17) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==17) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.4;
if (i==17) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.5;
if (i==17) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.51;
if (i==17) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.41;
if (i==17) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.41;
if (i==17) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.4;

if (i==17) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.2;


if (i==18) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.8;
if (i==18) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*1.2;
if (i==18) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.25;
if (i==18) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.07;
if (i==18) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.09;


if (i==18) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.33;
if (i==18) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.58;
if (i==18) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.83;
if (i==18) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.9;
if (i==18) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.76;
if (i==18) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.51;
if (i==18) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.4;
if (i==18) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.24;

if (i==18) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.2;


if (i==19) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.1;
if (i==19) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.33;
if (i==19) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.42;

if (i==19) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.4;
if (i==19) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.35;
if (i==19) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.27;
if (i==19) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.26;

if (i==19) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.18;
if (i==19) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.1;
if (i==19) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*0.6;






if (i==20) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.7;
if (i==20) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.7;
if (i==20) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.8;
if (i==20) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*0.8;
if (i==20) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==20) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.9;





if (i==20) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*1.15;
if (i==20) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.2;
if (i==20) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.18;
if (i==20) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.22;



if (i==20) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.41;
if (i==20) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.54;
if (i==20) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.52;

if (i==20) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.55;
if (i==20) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.47;
if (i==20) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.37;
if (i==20) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.32;

if (i==20) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.2;
if (i==20) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.07;
if (i==20) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*0.6;





if (i==21) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.5;
if (i==21) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.7;
if (i==21) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.75;
if (i==21) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.8;
if (i==21) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*0.9;
if (i==21) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==21) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.9;





if (i==21) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*1.15;
if (i==21) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.15;
if (i==21) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.05;
if (i==21) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.18;


if (i==21) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.19;
if (i==21) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.4;
if (i==21) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.48;

if (i==21) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.56;
if (i==21) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.52;
if (i==21) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.45;
if (i==21) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.3;

if (i==21) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.1;
if (i==21) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.07;
if (i==21) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*0.6;

//------------------------------------

//if (i==22) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.95;
//if (i==22) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.02;
if (i==22) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.97;

if (i==22) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.25;
if (i==22) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.42;
if (i==22) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.55;

if (i==22) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.7;
if (i==22) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.55;
if (i==22) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.45;
if (i==22) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.3;



if (i==23) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.67;
if (i==23) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.87;
if (i==23) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.82;

if (i==23) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.15;

if (i==23) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.25;
if (i==23) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==23) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.23;

if (i==23) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.3;
if (i==23) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.65;
if (i==23) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*2.1;
if (i==23) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*2.4;

if (i==23) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.8;


//if (i==23) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.85;


if (i==24) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.58;
if (i==24) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.91;
if (i==24) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.97;

if (i==24) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.15;

if (i==24) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.26;
if (i==24) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.27;
if (i==24) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.2;

if (i==24) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.13;
if (i==24) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.33;
if (i==24) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.55;
if (i==24) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.65;

if (i==24) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.6;


//if (i==24) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.85;


if (i==25) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.45;
if (i==25) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.95;
if (i==25) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.05;

if (i==25) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.45;

if (i==25) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.45;
if (i==25) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.45;
if (i==25) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.3;

if (i==25) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.15;
if (i==25) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.25;
if (i==25) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.6;
if (i==25) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*2.0;

if (i==25) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.85;
if (i==25) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.3;



if (i==26) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.92;
if (i==26) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.8;
if (i==26) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.1;

if (i==26) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.53;

if (i==26) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.6;
if (i==26) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.7;
if (i==26) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.63;

if (i==26) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.4;
if (i==26) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.35;
if (i==26) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.6;
if (i==26) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.7;

if (i==26) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.6;
if (i==26) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.3;





if (i==27) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.91;
if (i==27) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.9;
if (i==27) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.04;

if (i==27) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.3;

if (i==27) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.28;
if (i==27) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==27) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.3;

if (i==27) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.35;
if (i==27) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.3;
if (i==27) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.4;
if (i==27) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.7;

if (i==27) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.67;
if (i==27) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.3;



if (i==28) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.78;
//if (i==28) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.95;
if (i==28) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.13;

if (i==28) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.3;

if (i==28) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.5;
if (i==28) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.67;
if (i==28) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.75;

if (i==28) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.67;
if (i==28) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.6;
if (i==28) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.8;
if (i==28) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*2.6;

if (i==28) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*2.9;
if (i==28) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.75;

if (i==28) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*1.7;



if (i==29) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.64;
if (i==29) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.02;
if (i==29) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.33;

if (i==29) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.5;

if (i==29) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.65;
if (i==29) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*2.;
if (i==29) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.99;

if (i==29) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.9;
if (i==29) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.79;
if (i==29) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*2.15;
if (i==29) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*2.95;

if (i==29) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*3.0;
if (i==29) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*2.;

if (i==29) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*1.9;



};
};
};

/*
Xsect_int = 0.;


//Determine the with of the bin over all variables 
ds12 = (s12_max - s12_min)/15;
ds23 = (s23_max - s23_min)/15;
dalpha = (alp_max - alp_min)/13;
dtheta = (th_max - th_min)/13;

ds12_tmp = ds12; 
ds23_tmp = ds23; 
dalpha_tmp = dalpha; 
dtheta_tmp = dtheta;


dm12 = (sqrt(s12_max) - sqrt(s12_min))/15;
dm23 = (sqrt(s23_max) - sqrt(s23_min))/15;


cout << "\n";
cout << "s12_min = "<< s12_min<<", s12_max = "<< s12_max <<", ds12 = " << ds12<< "\n";
cout << "s23_min = "<< s23_min<<", s23_max = "<< s23_max <<", ds23 = " << ds23<< "\n";
cout << "th_min = "<< th_min<<", th_max = "<< th_max <<", dtheta = " << dtheta<< "\n";
cout << "alp_min = "<< alp_min<<", alp_max = "<< alp_max <<", dalpha = " << dalpha<< "\n";


for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {
//I am doing this to force this variables renew each time loops run, beceuse they are determied outside the loop and sometimes change inside the loop
ds12=ds12_tmp;
ds23=ds23_tmp;
dalpha = dalpha_tmp;
dtheta = dtheta_tmp;
cross_sect = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1];
//(S12_ARR_GOL[15][i]-S12_ARR_GOL[0][i])/(S23_ARR_GOL[15][i]-S23_ARR_GOL[0][i]);

//cout << W << "\n";

//if the point is the first or the last then the width of the bin is smaller 
if ((is12==1)||(is12==16)) ds12=ds12_tmp/2;
//if ((is12==16)) cout << ds12 << "\n";

if ((is23==1)||(is23==16)) ds23=ds23_tmp/2;
//if ((is23==14)) cout << ds23 << "\n";

if ((ialpha==1)||(ialpha==14)) dalpha=dalpha_tmp/2+0.01;
//if ((ialpha==2)) cout << dalpha << "\n";


//left and right edge for theta if the point is not first or last
th_l = th_min + dtheta/2.+(itheta - 2.)*dtheta;
th_r = th_min + dtheta/2.+(itheta -1.)*dtheta;

//left and right edge for theta if the point is first or last
if ((itheta==1)) th_l = th_min-0.01;
if ((itheta==1)) th_r = th_min+dtheta/2.;

if ((itheta==14)) th_l = th_max - dtheta/2.;
if ((itheta==14)) th_r = th_max+0.01;

Xsect_int = Xsect_int + cross_sect*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 

};
};
};
};



cout << i << " " << Xsect_int <<"\n";*/
};//end loop-golovach


for (Short_t i=0;i<=9;i++){
S12_ARR_FED_THRESH[i][0] = 0.07791914 + i*(0.089538593-0.07791914)/9.;
S12_ARR_FED_THRESH[i][1] = 0.07791914 + i*(0.105125093-0.07791914)/9.;
S12_ARR_FED_THRESH[i][2] = 0.07791914 + i*(0.121961593-0.07791914)/9.;

S23_ARR_FED_THRESH[i][0] = 1.161744 + i*(1.205450285 -1.161744 )/9.;
S23_ARR_FED_THRESH[i][1] = 1.161744 + i*(1.260971785 -1.161744 )/9.;
S23_ARR_FED_THRESH[i][2] = 1.161744 + i*(1.317743285 -1.161744 )/9.;

};

W_ARR_FED_THRESH[0] = 1.2375;
W_ARR_FED_THRESH[1] = 1.2625;
W_ARR_FED_THRESH[2] = 1.2875;

//FEDOTOV CROSS SECTION READING 
THETA_ARR_FED[0] = 0.;
THETA_ARR_FED[1] =0.4559418;
THETA_ARR_FED[2] =0.9018836;
THETA_ARR_FED[3] =1.347825;
THETA_ARR_FED[4] =1.793767;
THETA_ARR_FED[5] =2.239709;
THETA_ARR_FED[6] =2.685651;
THETA_ARR_FED[7] =3.141593;

ALPHA_ARR_FED[0] = 0.;
ALPHA_ARR_FED[1] = 0.9047408;
ALPHA_ARR_FED[2] = 1.799482;
ALPHA_ARR_FED[3] = 2.694222;
ALPHA_ARR_FED[4] = 3.588963;
ALPHA_ARR_FED[5] = 4.483704;
ALPHA_ARR_FED[6] = 5.378445;
ALPHA_ARR_FED[7] = 6.283185;


for (Short_t i=0; i<=5; i++) {
for (Int_t is23 = 1; is23 <=10; is23++) {
for (Int_t is12 = 1; is12 <=10; is12++) {
for (Int_t itheta = 1; itheta <=8; itheta++) {
for (Int_t ialpha = 1; ialpha <=8; ialpha++) {

SIGMA_ARR_FED[i][0][0][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][1][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][2][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][3][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][4][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][5][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][6][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][7][is23-1][is12-1][itheta-1][ialpha-1] = 0.;

SIGMA_ARR_FED[i][1][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;

SIGMA_ARR_FED[i][2][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][2][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;


SIGMA_ARR_FED[i][3][9][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][3][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][3][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;

SIGMA_ARR_FED[i][4][8][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][4][9][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][4][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][4][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;


SIGMA_ARR_FED[i][5][8][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][5][9][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][5][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][5][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;



SIGMA_ARR_FED[i][6][6][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][7][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][8][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][9][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;


};
};
};
};
};




for (Short_t i=0; i<=55; i++) {

if ((i>=0)&&(i<=3)){
q2bin = 0;
wbin = 8+i;

};

if ((i>=4)&&(i<=14)){
q2bin = 1;
wbin = i-4;
};

if ((i>=15)&&(i<=24)){
q2bin = 2;
wbin = i-15;
};

if ((i>=25)&&(i<=33)){
q2bin = 3;
wbin = i-25;
};

if ((i>=34)&&(i<=41)){
q2bin = 4;
wbin = i-34;
};

if ((i>=42)&&(i<=49)){
q2bin = 5;
wbin = i-42;
};

if ((i>=50)&&(i<=55)){
q2bin = 6;
wbin = i-50;
};


if ((q2bin >=0)&&(q2bin <=2)) Q2_ARR_FED[q2bin] = 0.225+0.05*q2bin;

if ((q2bin >=3)&&(q2bin <=6)) Q2_ARR_FED[q2bin] = 0.225+0.05*(q2bin+1);

W_ARR_FED[wbin] = 1.3125+0.025*wbin;
//cout << "Reading FEDOTOV diff cross sections for Q^2 = "<< Q2_ARR_FED[q2bin]<<" GeV^2 , W = "<< W_ARR_FED[wbin] << " GeV \n";
//cout << i<< "  "<< wbin <<"  "<< W_ARR_FED[wbin]<< "\n";
string dummy,xsect;

string file=file_names_fed[i];
ifstream input(file.c_str());
if(input.is_open()){

for (Int_t is23 = 1; is23 <=10; is23++) {
for (Int_t is12 = 1; is12 <=10; is12++) {
for (Int_t itheta = 1; itheta <=8; itheta++) {
for (Int_t ialpha = 1; ialpha <=8; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
S12_ARR_FED[is12-1][wbin] = atof(xsect.c_str());
//cout << S12_ARR_FED[is12-1][wbin] <<"\n";
getline(input,xsect);
S23_ARR_FED[is23-1][wbin] = atof(xsect.c_str());
//cout << S23_ARR_FED[is12-1][wbin] <<"\n";
getline(input,dummy);
getline(input,dummy);

getline(input,xsect);
//sigma_t
SIGMA_ARR_FED[0][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_l
SIGMA_ARR_FED[1][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_c2f
SIGMA_ARR_FED[2][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_s2f
SIGMA_ARR_FED[3][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_cf
SIGMA_ARR_FED[4][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_sf
SIGMA_ARR_FED[5][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,dummy);

getline(input,dummy);



};
};
};
};








};
input.close();

for (Int_t is23 = 1; is23 <=10; is23++) {
for (Int_t is12 = 1; is12 <=10; is12++) {
for (Int_t itheta = 1; itheta <=8; itheta++) {
for (Int_t ialpha = 1; ialpha <=8; ialpha++) {



for (Short_t j=0;j<6;j++){

SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*0.2;
SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*0.7;
SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*0.85;




if (q2bin==4) SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1]*1.07;
if (q2bin==6) SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1]*1.07;
if (q2bin==5) SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1]*1.07;



//SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_FED[9][wbin]-S12_ARR_FED[0][wbin])*(S23_ARR_FED[9][wbin]-S23_ARR_FED[0][wbin]);

SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*2.*M_PI*2.*M_PI;

if ((q2bin==0)&&(wbin==10)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.97;
if ((q2bin==0)&&(wbin==9)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.06;

if ((q2bin==1)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.80;
if ((q2bin==1)&&(wbin==7)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.93;
if ((q2bin==1)&&(wbin==8)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.93;
if ((q2bin==1)&&(wbin==9)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.93;
if ((q2bin==1)&&(wbin==10)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.98;



if ((q2bin==2)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.25;
if ((q2bin==2)&&(wbin==8)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.05;
if ((q2bin==2)&&(wbin==9)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.05;


if ((q2bin==3)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.1;
if ((q2bin==3)&&(wbin==7)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.96;

if ((q2bin==4)&&(wbin==7)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.97;

if ((q2bin==6)&&(wbin==5)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.92;


if ((q2bin==2)&&(wbin==4)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.02;
if ((q2bin==2)&&(wbin==2)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.97;


if ((q2bin==3)&&(wbin==4)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.98;
if ((q2bin==3)&&(wbin==2)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.08;

if ((q2bin==4)&&(wbin==3)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.05;


if ((q2bin==5)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.05;
if ((q2bin==5)&&(wbin==3)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.96;
if ((q2bin==5)&&(wbin==2)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.87;


//if ((q2bin==6)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.95;
//if ((q2bin==6)&&(wbin==3)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.07;
if ((q2bin==6)&&(wbin==4)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.08;

//---------------------------------------------------

//SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_FED_THRESH[9][0]-S12_ARR_FED_THRESH[0][0])*(S23_ARR_FED_THRESH[9][0]-S23_ARR_FED_THRESH[0][0]);

SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*2.*M_PI*2.*M_PI;



//SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_FED_THRESH[9][1]-S12_ARR_FED_THRESH[0][1])*(S23_ARR_FED_THRESH[9][1]-S23_ARR_FED_THRESH[0][1]);

SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*2.*M_PI*2.*M_PI;


//SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_FED_THRESH[9][2]-S12_ARR_FED_THRESH[0][2])*(S23_ARR_FED_THRESH[9][2]-S23_ARR_FED_THRESH[0][2]);

SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*2.*M_PI*2.*M_PI;





//SIGMA_ARR_FED_THRESH[j][q2bin][0][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*1.1;






//cout << SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]  << " "<< SIGMA_ARR_FED[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]<<"\n";
};

};
};
};
};



};



for (Short_t i=0; i<=10; i++) {
wbin = i;

W_ARR_RIP2[i] = 1.8375 + 0.025*i;

//cout<<i <<" "<< Int_t(i/17) << " "<<i - Int_t(i/17)*17<< "\n";
cout << "Reading RIPANI diff cross sections for Q^2 = 1.3 GeV^2, W = "<< W_ARR_RIP2[i] << " GeV \n";

string dummy,xsect;

string file=file_names_rip2[i];
ifstream input(file.c_str());
if(input.is_open()){

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
S12_val = atof(xsect.c_str());
S12_ARR_RIP2[is12-1][wbin] = atof(xsect.c_str());
getline(input,xsect);
S23_val = atof(xsect.c_str());
S23_ARR_RIP2[is23-1][wbin] = atof(xsect.c_str());
getline(input,dummy);
TH_val = atof(dummy.c_str());
getline(input,dummy);
ALP_val = atof(dummy.c_str());
getline(input,xsect);
//sigma_t
SIGMA_ARR_RIP2[0][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_l
SIGMA_ARR_RIP2[1][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_c2f
SIGMA_ARR_RIP2[2][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_s2f
SIGMA_ARR_RIP2[3][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_cf
SIGMA_ARR_RIP2[4][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_sf
SIGMA_ARR_RIP2[5][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
EPS_L_RIP2[wbin] = atof(xsect.c_str());
getline(input,dummy);

if ((is23==1)&&(itheta==1)&&(ialpha==1)&&(is12==1)) s12_min = S12_val;
if ((is23==1)&&(itheta==1)&&(ialpha==1)&&(is12==12)) s12_max = S12_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) s23_min = S23_val;
if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==12)) s23_max = S23_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) th_min = TH_val;
if ((is12==1)&&(itheta==6)&&(ialpha==1)&&(is23==1)) th_max = TH_val;

if ((is12==1)&&(itheta==1)&&(ialpha==1)&&(is23==1)) alp_min = ALP_val;
if ((is12==1)&&(itheta==1)&&(ialpha==6)&&(is23==1)) alp_max = ALP_val;

};
};
};
};
};
input.close();







for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {



for (Short_t j=0;j<6;j++){


//SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_RIP2[11][wbin]-S12_ARR_RIP2[0][wbin])*(S23_ARR_RIP2[11][wbin]-S23_ARR_RIP2[0][wbin]);

SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*2.*M_PI*2.*M_PI;

//SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*2.*M_PI;


if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.15*1.02;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.007;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.1*1.01;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.96*1.03;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.94*1.02;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.02;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.03*1.02;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.06*1.02;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.07*1.04;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.03;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.16*1.03;




};

};
};
};
};


for (Short_t j=0;j<6;j++){

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {


SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.5;

if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*1.1;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*1.05;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.3;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.35;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.1;


if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.6;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*0.55;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*0.55;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.55;
//if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*1.3;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.85;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.8;

if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.85;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.2;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;


if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.3;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.3;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
//if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.7;
//if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.5;


if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*1.2;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.2;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.2;

if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.4;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.4;
//if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.05;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.85;

//if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.6;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.3;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.3;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.9;





if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.85;
//if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*16;
//if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.4;
//if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.7;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;

if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.3;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*1.2;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;

if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.5;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.1;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.1;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.1;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.85;


};
};
};


for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.7;



if (i==0) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.2;
if (i==0) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.05;
if (i==0) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.75;
if (i==0) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.95;
//if (i==0) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*1.2;
if (i==0) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.15;
if (i==0) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.34;
if (i==0) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.25;


//if (i==1) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==1) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.9;
if (i==1) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
if (i==1) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==1) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.3;
if (i==1) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.25;


//if (i==2) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.95;
//if (i==2) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.3;
if (i==2) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.5;
if (i==2) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==2) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.85;
if (i==2) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*1.25;
if (i==2) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.35;
if (i==2) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.55;
if (i==2) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.35;


if (i==3) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.1;
if (i==3) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.7;
if (i==3) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.6;

if (i==3) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.7;
if (i==3) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.76;
if (i==3) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.3;
if (i==3) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.7;
if (i==3) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.5;

if (i==4) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.23;
if (i==4) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.13;
if (i==4) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.75;
if (i==4) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.7;
if (i==4) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.95;
if (i==4) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*1.03;
if (i==4) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.25;
if (i==4) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.25;

if (i==5) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.3;
if (i==5) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.7;
if (i==5) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.55;
if (i==5) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.7;
if (i==5) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.8;

if (i==5) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.45;
if (i==5) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.65;
if (i==5) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.55;



if (i==6) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
//if (i==6) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==6) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.6;
if (i==6) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==6) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.75;
if (i==6) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
if (i==6) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.35;
if (i==6) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.5;


if (i==7) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.65;
if (i==7) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.55;
if (i==7) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==7) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.8;
if (i==7) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.2;
if (i==7) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.5;
if (i==7) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.5;

if (i==8) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.95;
//if (i==8) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.1;
if (i==8) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.7;
if (i==8) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==8) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.7;
if (i==8) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.87;
if (i==8) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.03;
if (i==8) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.27;
if (i==8) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.6;


if (i==9) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.6;
if (i==9) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==9) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.7;
if (i==9) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.87;
if (i==9) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.4;
if (i==9) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.75;


if (i==10) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==10) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.9;
if (i==10) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.7;
if (i==10) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.67;
if (i==10) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.7;
if (i==10) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.87;
if (i==10) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.1;
if (i==10) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.55;
if (i==10) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.75;

};
};
};


for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.2;
SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*0.75;





if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*2.5;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.855;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.8;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.6;

if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*1.25;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.9;

if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.7;

if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*1.5;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.8;

if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.9;


if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*1.5;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.9;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.9;

if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.05;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;

if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.1;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.05;

};
};
};





for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {

if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0]*0.4;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5]*0.4;

if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*0.8;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*0.8;

if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*2.5;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*2.5;


if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;

if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0]*0.7;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5]*0.7;

if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.6;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.6;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;

if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.6;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.6;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;

if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.7;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.7;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.4;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.4;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0]*0.8;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5]*0.8;

if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.4;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.4;

};
};
};




};

Xsect_int = 0.;


//Determine the with of the bin over all variables 
ds12 = (s12_max - s12_min)/11;
ds23 = (s23_max - s23_min)/11;
dalpha = (alp_max - alp_min)/5;
dtheta = (th_max - th_min)/5;

ds12_tmp = ds12; 
ds23_tmp = ds23; 
dalpha_tmp = dalpha; 
dtheta_tmp = dtheta;


dm12 = (sqrt(s12_max) - sqrt(s12_min))/11;
dm23 = (sqrt(s23_max) - sqrt(s23_min))/11;


//cout << "\n";
//cout << "s12_min = "<< s12_min<<", s12_max = "<< s12_max <<", ds12 = " << ds12<< "\n";
//cout << "s23_min = "<< s23_min<<", s23_max = "<< s23_max <<", ds23 = " << ds23<< "\n";
//cout << "th_min = "<< th_min<<", th_max = "<< th_max <<", dtheta = " << dtheta<< "\n";
//cout << "alp_min = "<< alp_min<<", alp_max = "<< alp_max <<", dalpha = " << dalpha<< "\n";


for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
//I am doing this to force this variables renew each time loops run, beceuse they are determied outside the loop and sometimes change inside the loop
ds12=ds12_tmp;
ds23=ds23_tmp;
dalpha = dalpha_tmp;
dtheta = dtheta_tmp;
cross_sect_t = SIGMA_ARR_RIP2[0][wbin][is23-1][is12-1][itheta-1][ialpha-1];
cross_sect_l = SIGMA_ARR_RIP2[1][wbin][is23-1][is12-1][itheta-1][ialpha-1];
eps_l_rip2 = EPS_L_RIP2[wbin];
//(S12_ARR_GOL[15][i]-S12_ARR_GOL[0][i])/(S23_ARR_GOL[15][i]-S23_ARR_GOL[0][i]);

//cout << W << "\n";

//if the point is the first or the last then the width of the bin is smaller 
if ((is12==1)||(is12==12)) ds12=ds12_tmp/2;
//if ((is12==16)) cout << ds12 << "\n";

if ((is23==1)||(is23==12)) ds23=ds23_tmp/2;
//if ((is23==14)) cout << ds23 << "\n";

if ((ialpha==1)||(ialpha==6)) dalpha=dalpha_tmp/2+0.01;
//if ((ialpha==2)) cout << dalpha << "\n";


//left and right edge for theta if the point is not first or last
th_l = th_min + dtheta/2.+(itheta - 2.)*dtheta;
th_r = th_min + dtheta/2.+(itheta -1.)*dtheta;

//left and right edge for theta if the point is first or last
if ((itheta==1)) th_l = th_min-0.01;
if ((itheta==1)) th_r = th_min+dtheta/2.;

if ((itheta==6)) th_l = th_max - dtheta/2.;
if ((itheta==6)) th_r = th_max+0.01;
//cross_sect = cross_sect_t + eps_l_rip2*cross_sect_l;
//cross_sect = cross_sect_t;
cross_sect =cross_sect_l;
Xsect_int = Xsect_int + cross_sect*(cos(th_l)-cos(th_r))*ds12*ds23*dalpha; 

};
};
};
};



cout << i << " " << Xsect_int <<"  "<< EPS_L_RIP2[i] <<"\n";


};//end lop rip2 - i




 return;
 
};

