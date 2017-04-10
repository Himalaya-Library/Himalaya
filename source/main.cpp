#define Pi M_PI

#include <Eigen>
#include <HierarchyCalculator.hpp>
#include <iostream>
#include <H3m_interface.hpp>
#include <vector>
#include <math.h>

h3m::Parameters setup_SPS1a()
{
   h3m::Parameters pars;

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
   pars.Zt << 5.57720315E-01, 8.30028945E-01, 8.30028945E-01, -5.57720315E-01;
   pars.Zb << -9.33860207E-01, -3.57638243E-01, 3.57638243E-01, -9.33860207E-01;

   return pars;
}

h3m::Parameters setup_SPS2()
{
   h3m::Parameters pars;

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
   pars.Zt << 1.13197339E-01, 9.93572525E-01, 9.93572525E-01, -1.13197339E-01;
   pars.Zb << -9.99883015E-01, -1.52956306E-02, 1.52956306E-02, -9.99883015E-01;

   return pars;
}

h3m::Parameters setup_CMSSM_large_m0()
{
   h3m::Parameters pars;

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
   pars.Zt << 8.00864578E-02, 9.96787921E-01, 9.96787921E-01, -8.00864578E-02;
   pars.Zb << -9.99772435E-01, -2.13325618E-02, 2.13325618E-02, -9.99772435E-01;

   return pars;
}

h3m::Parameters setup_HSSUSY_minmix()
{
   h3m::Parameters pars;

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
   pars.Zt << 7.92936690E-02, 9.96851300E-01, 9.96851300E-01, -7.92936690E-02;
   pars.Zb << -3.29138831E-02, 9.99458191E-01, 9.99458191E-01, 3.29138831E-02;

   return pars;
}

h3m::Parameters checkH3m()
{
   h3m::Parameters pars;

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
   pars.MSb << 405.778, 594.116 ;
   pars.Zt << sin(0.983974), 1., 1., 1.;
   pars.Zb << 1., 1., 1., 1.;

   return pars;
}

int main(int argc, char **argv) {
   try{
      const std::vector<h3m::Parameters> points = {
	 setup_SPS1a(),
	 setup_SPS2(),
	 setup_CMSSM_large_m0(),
	 setup_HSSUSY_minmix(),
	 checkH3m()
      }; 
      for (const auto point: points) {
	 std::cout << "----------------------------------" << std::endl;
	 // init hierarchy calculator
	 h3m::HierarchyCalculator hierarchyCalculator(point);
	 // compare expanded terms at 2-loop level with the exact 2-loop result and choose a suitable hierarchy
	 const int suitableHierarchyTop = hierarchyCalculator.compareHierarchies(false);
	 const int suitableHierarchyBot = hierarchyCalculator.compareHierarchies(true);
	 // calculate the 3-loop corrections with the suiatble hierarchy
	 Eigen::Matrix2d DMh3L = hierarchyCalculator.calculateHierarchy(suitableHierarchyTop, false, 0, 0, 1);
	 Eigen::Matrix2d DMh3Lb = hierarchyCalculator.calculateHierarchy(suitableHierarchyBot, true, 0, 0, 1);
	 std::cout << "hierarchy top: " << suitableHierarchyTop << ", hierarchy bot: " << suitableHierarchyBot << std::endl;
	 std::cout << "DMh3L = " << DMh3L.row(0) << ' ' << DMh3L.row(1) << std::endl;
	 std::cout << "DMh3Lb = " << DMh3Lb.row(0) << ' ' << DMh3Lb.row(1) << std::endl;
      }
   }
   catch (std::exception& e){
      std::cout << e.what() << std::endl;
      return EXIT_FAILURE;
   }
   return 0;
}
