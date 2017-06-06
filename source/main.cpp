#define Pi M_PI

#include <Eigen>
#include <HierarchyCalculator.hpp>
#include <iostream>
#include <Himalaya_interface.hpp>
#include <HierarchyObject.hpp>
#include <vector>

himalaya::Parameters setup_SPS1a(){
   himalaya::Parameters pars;

   pars.scale = 4.67491329E+02;
   pars.mu = 3.52600579E+02;
   pars.g3 = 1.09949966E+00;
   pars.vd = 2.49832484E+01;
   pars.vu = 2.43650538E+02;
   pars.mq2 << 2.99793928E+05, 0, 0,
               0, 2.99792102E+05, 0,
               0, 0, 2.49327504E+05;
   pars.md2 << 2.78275669E+05, 0, 0,
               0, 2.78273780E+05, 0,
               0, 0, 2.74928741E+05;
   pars.mu2 << 2.80477426E+05, 0, 0,
               0, 2.80475621E+05, 0,
               0, 0, 1.80478484E+05;
   pars.Ab = -798.8142296644842;
   pars.At = -506.4162662374052;

   pars.MA = 3.92960954E+02;
   pars.MG = 5.88220143E+02;
   pars.MW = 8.04136643E+01;
   pars.MZ = 9.06817306E+01;
   pars.Mt = 1.52117491E+02;
   pars.Mb = 2.42010269E+00;
   pars.MSt << 3.83255677E+02, 5.70240743E+02;
   pars.MSb << 4.98792475E+02, 5.28677648E+02;
   pars.s2t = sin(2*asin(5.57720315E-01));
   pars.s2b = sin(2*asin(-9.33860207E-01));

   return pars;
}

himalaya::Parameters setup_SPS2(){
   himalaya::Parameters pars;

   pars.scale = 1.11090135E+03;
   pars.mu = 3.73337018E+02;
   pars.g3 = 1.06187116E+00;
   pars.vd = 2.51008404E+01;
   pars.vu = 2.41869332E+02;
   pars.mq2 << 2.36646981E+06, 0, 0,
               0, 2.36644973E+06, 0,
               0, 0, 1.63230152E+06;
   pars.md2 << 2.35612778E+06, 0, 0,
               0, 2.35610884E+06, 0,
               0, 0, 2.31917415E+06;
   pars.mu2 << 2.35685097E+06, 0, 0,
               0, 2.35682945E+06, 0,
               0, 0, 9.05923409E+05;
   pars.Ab = -784.3356416708631;
   pars.At = -527.8746242245387;
   
   pars.MA = 1.48446235E+03;
   pars.MG = 6.69045022E+02;
   pars.MW = 8.04001915E+01;
   pars.MZ = 8.97608307E+01;
   pars.Mt = 1.47685846E+02;
   pars.Mb = 2.38918959E+00;
   pars.MSt << 9.57566721E+02, 1.28878643E+03;
   pars.MSb << 1.27884964E+03, 1.52314587E+03;
   pars.s2t = sin(2*asin(1.13197339E-01));
   pars.s2b = sin(2*asin(-9.99883015E-01));

   return pars;
}

himalaya::Parameters setup_CMSSM_large_m0(){
   himalaya::Parameters pars;

   pars.scale = 4.88927977E+03;
   pars.mu = 3.25904948E+03;
   pars.g3 = 9.82503492E-01;
   pars.vd = 2.54061565E+01;
   pars.vu = 2.41061437E+02;
   pars.mq2 << 3.64876936E+07, 0, 0,
               0, 3.64874478E+07, 0,
               0, 0, 2.92101141E+07;
   pars.md2 << 3.36597459E+07, 0, 0,
               0, 3.36595040E+07, 0,
               0, 0, 3.32395009E+07;
   pars.mu2 << 3.39785780E+07, 0, 0,
               0, 3.39783229E+07, 0,
               0, 0, 1.95557229E+07;
   pars.Ab = -8205.625751354333;
   pars.At = -5328.586025475935;

   pars.MA = 4.76241507E+03;
   pars.MG = 5.95227200E+03;
   pars.MW = 8.03956029E+01;
   pars.MZ = 8.86740277E+01;
   pars.Mt = 1.37624174E+02;
   pars.Mb = 2.19689361E+00;
   pars.MSt << 4.41708777E+03, 5.41195051E+03;
   pars.MSb << 5.40476544E+03, 5.76558587E+03;
   pars.s2t = sin(2*asin(8.00864578E-02));
   pars.s2b = sin(2*asin(-9.99772435E-01));

   return pars;
}

himalaya::Parameters setup_HSSUSY_minmix(){
   himalaya::Parameters pars;

   pars.scale = 1.00000069E+05;
   pars.mu = 1.00012235E+05;
   pars.g3 = 8.78633324E-01;
   pars.vd = 5.13995860E+01;
   pars.vu = 2.36159115E+02;
   pars.mq2 << 1.00025551E+10, 0, 0,
               0, 1.00025551E+10, 0,
               0, 0, 1.00024663E+10;
   pars.md2 << 9.99992530E+09, 0, 0,
               0, 9.99992530E+09, 0,
               0, 0, 9.99993108E+09;
   pars.mu2 << 9.99990230E+09, 0, 0,
               0, 9.99990230E+09, 0,
               0, 0, 9.99961778E+09;
   pars.Ab = 50002.75311060441;
   pars.At = 19962.33330614816;

   pars.MA = 9.99320898E+04;
   pars.MG = 9.99989249E+04;
   pars.MW = 8.03334639E+01;
   pars.MZ = 8.70879717E+01;
   pars.Mt = 1.26274987E+02;
   pars.Mb = 2.06476960E+00;
   pars.MSt << 9.99980718E+04, 1.00012490E+05;
   pars.MSb << 9.99996447E+04, 1.00012359E+05;
   pars.s2t = sin(2*asin(7.92936690E-02));
   pars.s2b = sin(2*asin(-3.29138831E-02));

   return pars;
}

himalaya::Parameters checkH3m(){
   himalaya::Parameters pars;

   pars.scale = 173.3;
   pars.mu = 352.573;
   pars.g3 = 1.1218;
   pars.vd = 24.7443;
   pars.vu = 244.974;
   pars.mq2 << 330852, 0, 0,
               0, 330850, 0,
               0, 0, 274503;
   pars.md2 << 309155, 0, 0,
               0, 309153, 0,
               0, 0, 305087;
   pars.mu2 << 311332, 0, 0,
               0, 311330, 0,
               0, 0, 200009;
   pars.Ab = -529.657;
   pars.At = -529.657;

   pars.MA = 397.338;
   pars.MG = 615.162;
   pars.MW = 80.399;
   pars.MZ = 91.2;
   pars.Mt = 153.608;
   pars.Mb = 4.25;
   pars.MSt << 405.778, 594.116;
   pars.MSb << 405.778, 594.116;
   pars.s2t = sin(2*0.983974);
   pars.s2b = 1.;

   return pars;
}

himalaya::Parameters checkH3m2(){

   himalaya::Parameters pars;

   pars.scale = 1056.26;
   pars.mu = 134.903;
   pars.g3 = 1.07022;
   pars.vd = 25.01481631;
   pars.vu = 241.4436584;
   pars.mq2 << pow(1468.7,2), 0, 0,
               0, pow(1468.7,2), 0,
               0, 0, pow(1214.45,2);
   pars.md2 << pow(1470.63,2), 0, 0,
               0, pow(1470.62,2), 0,
               0, 0, pow(1458.97,2);
   pars.mu2 << pow(1470.37,2), 0, 0,
               0, pow(1470.36,2), 0,
               0, 0, pow(903.144,2);
   pars.Ab = -528.861551;
   pars.At = -361.862050;

   pars.MA = 1437.52;
   pars.MG = 450.708;
   pars.MW = 78.3505;
   pars.MZ = 89.8558;
   pars.Mt = 145.353;
   pars.Mb = 4.25;
   pars.MSt << 911.608, 1223.86;
   pars.MSb << 1216.13, 1459.28;
   pars.s2t = 0.163855;
   pars.s2b = 0.127886;

   return pars;
}

himalaya::Parameters Xt29(){

   himalaya::Parameters pars;

   pars.scale = 1981.84;
   pars.mu = 1999.88;
   pars.g3 = 1.02887;
   pars.vd = 49.723;
   pars.vu = 237.098;
   pars.mq2 <<  4.00295e+06, 0, 0,
               0, 4.00295e+06, 0,
               0, 0, 3.99948e+06;
   pars.md2 << 4.00249e+06, 0, 0,
               0, 4.00249e+06, 0,
               0, 0, 4.00239e+06;
   pars.mu2 << 4.00251e+06, 0, 0,
               0, 4.00251e+06 , 0,
               0, 0, 3.99551e+06;
   pars.Ab = 9997.88;
   pars.At = 6195.23;

   pars.MA = 1995.22;
   pars.MG = 2000.66;
   pars.MW = 77.0967;
   pars.MZ = 88.7893;
   pars.Mt = 145.817;
   pars.Mb = 2.2495;
   pars.MSt << 1781.78, 2204.35;
   pars.MSb << 2000.42, 2000.95;
   pars.s2t = -0.999999;
   pars.s2b = -0.973654;

   return pars;
}

himalaya::Parameters Xt3(){

   himalaya::Parameters pars;

   pars.scale = 1979.99;
   pars.mu = 1999.87;
   pars.g3 = 1.02891;
   pars.vd = 49.6913;
   pars.vu = 236.886;
   pars.mq2 << 4.00326e+06 , 0, 0,
               0, 4.00326e+06, 0,
               0, 0, 3.99918e+06;
   pars.md2 << 4.00275e+06, 0, 0,
               0, 4.00275e+06, 0,
               0, 0, 4.00263e+06;
   pars.mu2 << 4.00277e+06, 0, 0,
               0,4.00277e+06, 0,
               0, 0, 3.99455e+06;
   pars.Ab = 9997.64;
   pars.At = 6394.6;

   pars.MA = 1994.56;
   pars.MG = 2000.73;
   pars.MW = 77.0275;
   pars.MZ = 88.7096;
   pars.Mt = 146.148;
   pars.Mb = 2.24524;
   pars.MSt << 1772.9, 2211.26;
   pars.MSb << 2000.39, 2000.97;
   pars.s2t = -0.999998;
   pars.s2b = -0.895782;

   return pars;
}

himalaya::Parameters Xt31(){

   himalaya::Parameters pars;

   pars.scale = 1978.03;
   pars.mu = 1999.86;
   pars.g3 =  1.02896;
   pars.vd = 49.6565;
   pars.vu = 236.653;
   pars.mq2 <<  4.00358e+06 , 0, 0,
               0, 4.00358e+06, 0,
               0, 0, 3.99882e+06;
   pars.md2 << 4.00302e+06, 0, 0,
               0, 4.00302e+06, 0,
               0, 0, 4.00289e+06;
   pars.mu2 << 4.00304e+06 , 0, 0,
               0, 4.00304e+06, 0,
               0, 0, 3.99344e+06;
   pars.Ab = 9997.39;
   pars.At = 6593.92;

   pars.MA = 1993.83;
   pars.MG = 2000.8;
   pars.MW = 76.9519;
   pars.MZ = 88.6226;
   pars.Mt = 146.502;
   pars.Mb = 2.24085;
   pars.MSt << 1763.87 , 2218.18;
   pars.MSb <<  2000.33, 2001;
   pars.s2t = -0.999997;
   pars.s2b =-0.783043;

   return pars;
}

himalaya::Parameters Xt33(){

   himalaya::Parameters pars;

   pars.scale = 1973.75;
   pars.mu = 1999.82;
   pars.g3 =  1.02907;
   pars.vd = 49.5751;
   pars.vu = 236.115;
   pars.mq2 <<  4.00428e+06 , 0, 0,
               0, 4.00428e+06, 0,
               0, 0, 3.99786e+06;
   pars.md2 << 4.00361e+06, 0, 0,
               0, 4.00361e+06, 0,
               0, 0, 4.00346e+06;
   pars.mu2 << 4.00363e+06 , 0, 0,
               0, 4.00363e+06, 0,
               0, 0, 3.99067e+06;
   pars.Ab = 9996.81;
   pars.At = 6992.34;

   pars.MA = 1992.14;
   pars.MG = 2000.96;
   pars.MW = 76.7777;
   pars.MZ = 88.4219;
   pars.Mt = 147.295;
   pars.Mb = 2.23149;
   pars.MSt << 1745.3 , 2232.1;
   pars.MSb <<  2000.14, 2001.09;
   pars.s2t = -0.999995;
   pars.s2b =-0.550527;

   return pars;
}

himalaya::Parameters MS350(){

   himalaya::Parameters pars;

   pars.scale = 379.219;
   pars.mu = 350.21;
   pars.g3 =  1.1223;
   pars.vd = 49.3791;
   pars.vu = 241.173;
   pars.mq2 <<  121599, 0, 0,
               0, 121599, 0,
               0, 0, 121716;
   pars.md2 << 121725, 0, 0,
               0, 121725, 0,
               0, 0, 121762;
   pars.mu2 << 121722 , 0, 0,
               0, 121722, 0,
               0, 0, 121926;
   pars.Ab = 1752.9;
   pars.At = 73.137;

   pars.MA = 351.585;
   pars.MG = 348.783;
   pars.MW = 79.4527;
   pars.MZ = 90.9831;
   pars.Mt = 154.569;
   pars.Mb = 2.56808;
   pars.MSt << 378.123 , 380.319;
   pars.MSb <<  349.809, 353.456;
   pars.s2t =  -0.265985;
   pars.s2b = -0.0849635;

   return pars;
}

himalaya::Parameters MS400(){

   himalaya::Parameters pars;

 pars.scale   = 425.41;
  pars.mu  = 400.179;
  pars.g3  = 1.11453;
  pars.vd  = 49.4249;
  pars.vu  = 241.011;
  pars.mq2 << 159105 ,     0,      0,
     0 ,159105   ,   0,
     0    ,  0, 159222;
  pars.md2 << 159231  ,    0 ,     0,
     0, 159231     , 0,
     0 ,     0, 159267;
  pars.mu2 << 159228  ,    0 ,     0,
     0, 159228 ,     0,
     0     , 0 ,159432;
  pars.At  = 82.7215;
  pars.Ab  = 2002.52;
  pars.MG  = 398.945;
  pars.MW  = 79.3087;
  pars.MZ  = 90.858;
  pars.Mt  = 153.446;
  pars.Mb  = 2.54742;
  pars.MA  = 401.38;
  pars.MSt << 424.466, 426.356;
  pars.MSb << 399.836, 403.023;
  pars.s2t = -0.125094;
  pars.s2b = -0.101822;

   return pars;
}

himalaya::Parameters MS480(){

  himalaya::Parameters pars;

  pars.scale   = 500.92;
  pars.mu  = 480.143;
  pars.g3  = 1.104;
  pars.vd  = 49.4878;
  pars.vu  = 240.782;
  pars.mq2 << 229519    ,  0 ,     0,
     0 ,229519    ,  0,
     0  ,    0, 229635;
  pars.md2 << 229645,      0,      0,
     0, 229645     , 0,
     0 ,     0, 229680;
  pars.mu2 << 229641 ,     0     , 0,
     0, 229641      ,0,
     0 ,     0, 229843;
  pars.At  = 98.2263;
  pars.Ab  = 2402.06;
  pars.MG  = 479.138;
  pars.MW  = 79.1151;
  pars.MZ  = 90.6904;
  pars.Mt  = 151.895;
  pars.Mb  = 2.51902;
  pars.MA  = 481.136;
  pars.MSt << 500.129, 501.711;
  pars.MSb << 479.872, 482.523;
  pars.s2t = 0.0876104;
  pars.s2b = -0.130184;

   return pars;
}

himalaya::Parameters bug(){

  himalaya::Parameters pars;
   pars.scale   = 1158.13;
    pars.mu  = 1150.17;
    pars.g3  = 1.05606;
    pars.vd  = 49.7917;
    pars.vu  = 239.757;
    pars.mq2 << 1.32202e+06 ,          0 ,          0,
            0, 1.32202e+06    ,       0,
            0 ,          0, 1.32212e+06;
    pars.md2 << 1.32214e+06,           0 ,          0,
            0, 1.32214e+06     ,      0,
            0 ,          0, 1.32217e+06;
    pars.mu2 << 1.32213e+06,           0,           0,
            0, 1.32213e+06    ,       0,
            0 ,          0, 1.32231e+06;
    pars.At  = 230.831;
    pars.Ab  = 5751.4;
    pars.MG  = 1149.82;
    pars.MW  = 78.2546;
    pars.MZ  = 89.962;
    pars.Mt  = 144.595;
    pars.Mb  = 2.38715;
    pars.MA  = 1150.55;
    pars.MSt << 1157.54, 1158.73;
    pars.MSb << 1150.08, 1151.23;
    pars.s2t = 0.841368;
    pars.s2b = -0.381861;
     return pars;
}

himalaya::Parameters bug2(){
  himalaya::Parameters pars;
 pars.scale   = 1158.13;
  pars.mu  = 1150.17;
  pars.g3  = 1.05606;
  pars.vd  = 49.7917;
  pars.vu  = 239.757;
  pars.mq2 << 1.32202e+06,           0,           0,
          0, 1.32202e+06       ,    0,
          0 ,          0, 1.32212e+06;
  pars.md2 << 1.32214e+06,           0 ,          0,
          0, 1.32214e+06    ,       0,
          0 ,          0, 1.32217e+06;
  pars.mu2 << 1.32213e+06 ,          0  ,         0,
          0 ,1.32213e+06   ,        0,
          0  ,         0, 1.32231e+06;
  pars.At  = 230.831;
  pars.Ab  = 5751.4;
  pars.MG  = 1149.82;
  pars.MW  = 78.2546;
  pars.MZ  = 89.962;
  pars.Mt  = 144.595;
  pars.Mb  = 2.38715;
  pars.MA  = 1150.55;
  pars.MSt << 1157.54, 1158.73;
  pars.MSb << 1150.08, 1151.23;
  pars.s2t = 0.841368;
  pars.s2b = -0.381861;
    return pars;

}

int main(int argc, char **argv) {
   try{
      const std::vector<himalaya::Parameters> points = {
	 //setup_SPS1a(),
	 setup_SPS2(),
	 //setup_CMSSM_large_m0(),
	 //setup_HSSUSY_minmix(),
	 //checkH3m(),
	 //checkH3m2(),
	 //Xt29(),
	 //Xt3(),
	 //Xt31(),
	 //Xt33(),
	 //MS350(),
	 //MS400(),
	 //MS480(),
	 //bug(),
	 //bug2()
      }; 
      for (const auto point: points) {
	 // init hierarchy calculator
	 himalaya::HierarchyCalculator hierarchyCalculator(point);

	 // calculate the 3-loop corrections with the suiatble hierarchy
	 //top
	 himalaya::HierarchyObject hoTop = hierarchyCalculator.calculateDMh3L(false);

	 //bottom
	 himalaya::HierarchyObject hoBot = hierarchyCalculator.calculateDMh3L(true);
	 
	 // check terms
	 //hierarchyCalculator.checkTerms();
	 
	 /*std::cout << "hierarchy top: " << hoTop.getSuitableHierarchy() << ", hierarchy bot: " << hoBot.getSuitableHierarchy() << std::endl;
	 std::cout << "error top " << hoTop.getRelDiff2L() << " error bot: " << hoBot.getRelDiff2L() << std::endl;
	 std::cout << "abs err top " << hoTop.getAbsDiff2L() << " abs err bot " << hoBot.getAbsDiff2L() << std::endl; 
	 std::cout << "mdr " << hoTop.getMDRMasses() << std::endl;
	 std::cout << "tree " << hoTop.getDMh(0) << std::endl;
	 std::cout << "1l " << hoTop.getDMh(1) << std::endl;
	 std::cout << "2l " << hoTop.getDMh(2) << std::endl;
	 std::cout << "exp 1 " << hoTop.getExpUncertainty(1) << std::endl;
	 std::cout << "exp 2 " << hoTop.getExpUncertainty(2) << std::endl;
	 std::cout << "exp 3 " << hoTop.getExpUncertainty(3) << std::endl;
	 std::cout << "shift " << hoTop.getDRToMDRShift() << std::endl;
	 std::cout << "3l " << hoTop.getDMh(3) << std::endl;*/
	 std::cout << "----------------------------------" << std::endl;
      }
   }
   catch (std::exception& e){
      std::cout << e.what() << std::endl;
      return EXIT_FAILURE;
   }
   return 0;
}
