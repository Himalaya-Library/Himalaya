#define BOOST_TEST_MODULE test_FO_expressions

#include <boost/test/unit_test.hpp>
#define private public
#include "HierarchyCalculator.hpp"
#undef private

himalaya::Parameters make_point()
{
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
   pars.s2b = -0.550527;

   return pars;
}


/**
 * Performs a sanity check of the implemented expansion terms by
 * comparing them to their numerical value at a given parameter point.
 */
BOOST_AUTO_TEST_CASE(test_FO_expansions)
{
   using namespace himalaya;

   const auto pars = make_point();

   for (int i = 0; i < NUMBER_OF_HIERARCHIES; i++) {
      himalaya::HierarchyCalculator hc(pars);
      hc.init();

      himalaya::HierarchyObject ho(false);
      ho.setMDRFlag(1);
      ho.setSuitableHierarchy(i);

      Eigen::Matrix2d oloMat  = hc.calculateHierarchy(ho, 1,0,0);
      Eigen::Matrix2d twloMat = hc.calculateHierarchy(ho, 0,1,0);
      Eigen::Matrix2d thloMat = hc.calculateHierarchy(ho, 0,0,1);

      bool ck1LPassed = false, ck2LPassed = false, ck3LPassed = false;

      switch(i){
      case h3:
         ck1LPassed = abs(oloMat(0,0) - (-1033.437882123761)) < 1e-06 && abs(oloMat(1,0) - (-394.3521101999062)) < 1e-06 && abs(oloMat(1,1) - 17633.47392819223) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-13.48340821650015)) < 1e-04 && abs(twloMat(1,0) - 11.12436787252288) < 1e-04 && abs(twloMat(1,1) - 1476.660068002361) < 1e-03;
         ck3LPassed = abs(thloMat(0,0) - 1.096612614742133) < 1e-06 && abs(thloMat(1,0) - 9.986750150481939) < 1e-06 && abs(thloMat(1,1) - 370.2505433664134) < 1e-06;
	 break;
      case h32q2g:
         ck1LPassed = abs(oloMat(0,0) - (-1033.437882123761)) < 1e-06 && abs(oloMat(1,0) - (-394.3521101999062)) < 1e-06 && abs(oloMat(1,1) - 17633.47392819223) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-13.66052379180129)) < 1e-06 && abs(twloMat(1,0) - 11.26755617866339) < 1e-06 && abs(twloMat(1,1) - 1477.465656153518) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 1.113051431370291) < 1e-06 && abs(thloMat(1,0) - 9.903809573970422) < 1e-06 && abs(thloMat(1,1) - 369.7408109643386) < 1e-06;
	 break;
      case h3q22g:
         ck1LPassed = abs(oloMat(0,0) - (-1033.437882123761)) < 1e-06 && abs(oloMat(1,0) - (-394.3521101999062)) < 1e-06 && abs(oloMat(1,1) - 17633.47392819223) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-13.66052379180129)) < 1e-06 && abs(twloMat(1,0) - 11.26755617866339) < 1e-06 && abs(twloMat(1,1) - 1477.465656153518) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 1.058450932536496) < 1e-06 && abs(thloMat(1,0) - 10.0141272838662) < 1e-06 && abs(thloMat(1,1) - 370.3301180635573) < 1e-06;
	 break;
      case h4:
         ck1LPassed = abs(oloMat(0,0) - 0) < 1e-06 && abs(oloMat(1,0) - 0) < 1e-06 && abs(oloMat(1,1) - 6685.123085628641) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - 0) < 1e-06 && abs(twloMat(1,0) - 1183.325484493686) < 1e-06 && abs(twloMat(1,1) - 1458.970501474495) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 162.1379208650191) < 1e-06 && abs(thloMat(1,0) - 326.0219627343553) < 1e-06 && abs(thloMat(1,1) - 431.6926278454841) < 1e-06;
	 break;
      case h5:
         ck1LPassed = abs(oloMat(0,0) - 15921.69462848581) < 1e-06 && abs(oloMat(1,0) - (-388569.2043081555)) < 1e-06 && abs(oloMat(1,1) - 7874.401574063407) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-86.77887344841422)) < 1e-06 && abs(twloMat(1,0) - (-20625.63783863484)) < 1e-06 && abs(twloMat(1,1) - (-42446.62009872038)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 2442.115080578889) < 1e-06 && abs(thloMat(1,0) - (-3859.942907446577)) < 1e-06 && abs(thloMat(1,1) - 60593.055768119) < 1e-06;
	 break;
      case h5g1:
         ck1LPassed = abs(oloMat(0,0) - 15921.69462848581) < 1e-06 && abs(oloMat(1,0) - (-388569.2043081556)) < 1e-06 && abs(oloMat(1,1) - 7874.401574063407) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-114.6037388932203)) < 1e-06 && abs(twloMat(1,0) - (-20341.84471909946)) < 1e-06 && abs(twloMat(1,1) - (-42843.48046642416)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 2415.507513838155) < 1e-06 && abs(thloMat(1,0) - (-3766.750163753644)) < 1e-06 && abs(thloMat(1,1) - 59380.34497121828) < 1e-06;
	 break;
      case h6:
         ck1LPassed = abs(oloMat(0,0) - 9272.477351702315) < 1e-06 && abs(oloMat(1,0) - (-184.7601614832763)) < 1e-06 && abs(oloMat(1,1) - 7581.278122072418) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-1078.578574572312)) < 1e-06 && abs(twloMat(1,0) - 7096.529601647042) < 1e-06 && abs(twloMat(1,1) - (-1927.791631086123)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 245.4412216221288) < 1e-06 && abs(thloMat(1,0) - 573.1296253278389) < 1e-06 && abs(thloMat(1,1) - 8448.4582538127) < 1e-06;
	 break;
      case h6b:
         ck1LPassed = abs(oloMat(0,0) - 9272.477351702311) < 1e-06 && abs(oloMat(1,0) - (-184.7601614832763)) < 1e-06 && abs(oloMat(1,1) - 7581.278122072418) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-1078.578574572312)) < 1e-06 && abs(twloMat(1,0) - 7096.52960164704) < 1e-06 && abs(twloMat(1,1) - (-1900.197036824461)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 283.0253770519464) < 1e-06 && abs(thloMat(1,0) - 566.2182257407396) < 1e-06 && abs(thloMat(1,1) - 10093.33785879814) < 1e-06;
	 break;
      case h6b2qg2:
         ck1LPassed = abs(oloMat(0,0) - 9272.477351702311) < 1e-06 && abs(oloMat(1,0) - (-184.7601614832759)) < 1e-06 && abs(oloMat(1,1) - 7581.278122072418) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-1089.201418061661)) < 1e-06 && abs(twloMat(1,0) - 7145.267026465748) < 1e-06 && abs(twloMat(1,1) - (-2077.345120153528)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 285.3154791763894) < 1e-06 && abs(thloMat(1,0) - 544.3654284413091) < 1e-06 && abs(thloMat(1,1) - 10336.22756889787) < 1e-06;
	 break;
      case h6bq22g:
         ck1LPassed = abs(oloMat(0,0) - 9272.477351702315) < 1e-06 && abs(oloMat(1,0) - (-184.7601614832763)) < 1e-06 && abs(oloMat(1,1) - 7581.278122072418) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-1078.578574572311)) < 1e-06 && abs(twloMat(1,0) - 7096.529601647042) < 1e-06 && abs(twloMat(1,1) - (-1900.197036824461)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 283.0220052455883) < 1e-06 && abs(thloMat(1,0) - 566.2190953470737) < 1e-06 && abs(thloMat(1,1) - 10093.33986048966) < 1e-06;
	 break;
      case h6bq2g2:
         ck1LPassed = abs(oloMat(0,0) - 9272.477351702315) < 1e-06 && abs(oloMat(1,0) - (-184.7601614832759)) < 1e-06 && abs(oloMat(1,1) - 7581.278122072418) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-1089.201418061661)) < 1e-06 && abs(twloMat(1,0) - 7145.267026465748) < 1e-06 && abs(twloMat(1,1) - (-2077.345120153528)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 285.3120881213721) < 1e-06 && abs(thloMat(1,0) - 544.3662758149513) < 1e-06 && abs(thloMat(1,1) - 10336.23012077387) < 1e-06;
	 break;
      case h6g2:
         ck1LPassed = abs(oloMat(0,0) - 9272.477351702315) < 1e-06 && abs(oloMat(1,0) - (-184.7601614832761)) < 1e-06 && abs(oloMat(1,1) - 7581.278122072418) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - (-1089.201418061661)) < 1e-06 && abs(twloMat(1,0) - 7145.267026465748) < 1e-06 && abs(twloMat(1,1) - (-2112.642999123034)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 246.0217489966267) < 1e-06 && abs(thloMat(1,0) - 557.451210096066) < 1e-06 && abs(thloMat(1,1) - 8628.076480526881) < 1e-06;
	 break;
      case h9:
         ck1LPassed = abs(oloMat(0,0) - (-1033.437882123761)) < 1e-06 && abs(oloMat(1,0) - (-394.352110199906)) < 1e-06 && abs(oloMat(1,1) - 17633.47392819223) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - 420.2050380976995) < 1e-06 && abs(twloMat(1,0) - (-554.6021924866435)) < 1e-06 && abs(twloMat(1,1) - (-797.8089039452509)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 132.8584579769461) < 1e-06 && abs(thloMat(1,0) - (-171.9326869339159)) < 1e-06 && abs(thloMat(1,1) - (-800.8408283898472)) < 1e-06;
	 break;
      case h9q2:
         ck1LPassed = abs(oloMat(0,0) - (-1033.437882123761)) < 1e-06 && abs(oloMat(1,0) - (-394.3521101999065)) < 1e-06 && abs(oloMat(1,1) - 17633.47392819223) < 1e-06;
         ck2LPassed = abs(twloMat(0,0) - 420.2050380976993) < 1e-06 && abs(twloMat(1,0) - (-554.6021924866436)) < 1e-06 && abs(twloMat(1,1) - (-797.8089039452487)) < 1e-06;
         ck3LPassed = abs(thloMat(0,0) - 132.6358855624267) < 1e-06 && abs(thloMat(1,0) - (-171.4711818838455)) < 1e-06 && abs(thloMat(1,1) - (-800.9569014303727)) < 1e-06;
	 break;
      default:
         BOOST_FAIL("unknown hierarchy!");
         break;
      }

      BOOST_TEST_MESSAGE("Checking iierarchy " << i << " (" << ho.getH3mHierarchyNotation(i) << ")");
      BOOST_CHECK(ck1LPassed);
      BOOST_CHECK(ck2LPassed);
      BOOST_CHECK(ck3LPassed);
   }
}
