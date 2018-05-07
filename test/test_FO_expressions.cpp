#include "doctest.h"
#define private public
#include "HierarchyCalculator.hpp"
#undef private
#include "Hierarchies.hpp"

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

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
TEST_CASE("test_FO_expansions")
{
   using namespace himalaya;

   const double eps = 1e-13;

   for (int i = 0; i < Hierarchies::NUMBER_OF_HIERARCHIES; i++) {
      himalaya::HierarchyCalculator hc(make_point());
      hc.init();

      himalaya::HierarchyObject ho(false);
      ho.setRenormalizationScheme(RenSchemes::H3m);
      ho.setMDRFlag(1);
      ho.setSuitableHierarchy(i);

      const auto oloMat  = hc.calculateHierarchy(ho, 1, 0, 0);
      const auto twloMat = hc.calculateHierarchy(ho, 0, 1, 0);
      const auto thloMat = hc.calculateHierarchy(ho, 0, 0, 1);
      ho.setMDRFlag(0);
      hc.calculateHierarchy(ho, 0, 0, 1);
      const auto thlomh2 = ho.getZetaHimalaya();
      const auto hier_str = ho.getH3mHierarchyNotation(i);

      INFO("Checking hierarchy " << i << " (" << hier_str << ")");

      switch(i){
      case Hierarchies::h3:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.3521101999062, eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0), -13.48340821650015, 1e-5);
         CHECK_CLOSE(twloMat(1,0),  11.12436787252288, 1e-5);
         CHECK_CLOSE(twloMat(1,1),  1476.660068002361, 1e-5);

         CHECK_CLOSE(thloMat(0,0), 1.096612614742133, 2e-6);
         CHECK_CLOSE(thloMat(1,0), 9.986750150481939, 2e-6);
         CHECK_CLOSE(thloMat(1,1), 370.2505433664134, 2e-6);
	 CHECK_CLOSE(thlomh2, 7795.3852300944, eps);
	 break;
      case Hierarchies::h32q2g:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.3521101999062, eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0), -13.66052379180129, eps);
         CHECK_CLOSE(twloMat(1,0),  11.26755617866339, eps);
         CHECK_CLOSE(twloMat(1,1),  1477.465656153518, eps);

         CHECK_CLOSE(thloMat(0,0), 1.113051431370291, eps);
         CHECK_CLOSE(thloMat(1,0), 9.903809573970422, eps);
         CHECK_CLOSE(thloMat(1,1), 369.7408109643386, eps);
	 CHECK_CLOSE(thlomh2, 7669.5076452224, eps);
         break;
      case Hierarchies::h3q22g:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.3521101999062, eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0), -13.66052379180129, eps);
         CHECK_CLOSE(twloMat(1,0),  11.26755617866339, eps);
         CHECK_CLOSE(twloMat(1,1),  1477.465656153518, eps);

         CHECK_CLOSE(thloMat(0,0), 1.058450932536496, eps);
         CHECK_CLOSE(thloMat(1,0), 10.0141272838662, eps);
         CHECK_CLOSE(thloMat(1,1), 370.3301180635573, eps);
	 CHECK_CLOSE(thlomh2, 7823.8598939601, eps);
         break;
      case Hierarchies::h4:
         CHECK_CLOSE(oloMat(0,0),                 0, eps);
         CHECK_CLOSE(oloMat(1,0),                 0, eps);
         CHECK_CLOSE(oloMat(1,1), 6685.123085628641, eps);

         CHECK_CLOSE(twloMat(0,0),                 0, eps);
         CHECK_CLOSE(twloMat(1,0), 1183.325484493686, eps);
         CHECK_CLOSE(twloMat(1,1), 1458.970501474495, eps);

         CHECK_CLOSE(thloMat(0,0), 162.1379208650191, eps);
         CHECK_CLOSE(thloMat(1,0), 326.0219627343553, eps);
         CHECK_CLOSE(thloMat(1,1), 431.6926278454841, eps);
	 CHECK_CLOSE(thlomh2, -1422.6632313926, eps);
         break;
      case Hierarchies::h5:
         CHECK_CLOSE(oloMat(0,0),  15921.69462848581, eps);
         CHECK_CLOSE(oloMat(1,0), -388569.2043081555, eps);
         CHECK_CLOSE(oloMat(1,1),  7874.401574063407, eps);

         CHECK_CLOSE(twloMat(0,0), -86.77887344841422, eps);
         CHECK_CLOSE(twloMat(1,0), -20625.63783863484, eps);
         CHECK_CLOSE(twloMat(1,1), -42446.62009872038, eps);

         CHECK_CLOSE(thloMat(0,0),  2442.115080578889, eps);
         CHECK_CLOSE(thloMat(1,0), -3859.942907446577, eps);
         CHECK_CLOSE(thloMat(1,1),  60593.055768119  , eps);
	 CHECK_CLOSE(thlomh2, 12162367.7903137561, eps);
         break;
      case Hierarchies::h5g1:
         CHECK_CLOSE(oloMat(0,0),  15921.69462848581, eps);
         CHECK_CLOSE(oloMat(1,0), -388569.2043081556, eps);
         CHECK_CLOSE(oloMat(1,1),  7874.401574063407, eps);

         CHECK_CLOSE(twloMat(0,0), -114.6037388932203, eps);
         CHECK_CLOSE(twloMat(1,0), -20341.84471909946, eps);
         CHECK_CLOSE(twloMat(1,1), -42843.48046642416, eps);

         CHECK_CLOSE(thloMat(0,0),  2415.507513838155, eps);
         CHECK_CLOSE(thloMat(1,0), -3766.750163753644, eps);
         CHECK_CLOSE(thloMat(1,1),  59380.34497121828, eps);
	 CHECK_CLOSE(thlomh2, 11936978.6994632017, eps);
         break;
      case Hierarchies::h6:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702315, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832763, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1078.578574572312, eps);
         CHECK_CLOSE(twloMat(1,0),  7096.529601647042, eps);
         CHECK_CLOSE(twloMat(1,1), -1927.791631086123, eps);

         CHECK_CLOSE(thloMat(0,0), 245.4412216221288, eps);
         CHECK_CLOSE(thloMat(1,0), 573.1296253278389, eps);
         CHECK_CLOSE(thloMat(1,1), 8448.4582538127  , eps);
	 CHECK_CLOSE(thlomh2, 3575581.6748945541, eps);
         break;
      case Hierarchies::h6b:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702311, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832763, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1078.578574572312, eps);
         CHECK_CLOSE(twloMat(1,0),  7096.52960164704 , eps);
         CHECK_CLOSE(twloMat(1,1), -1900.197036824461, eps);

         CHECK_CLOSE(thloMat(0,0), 283.0253770519464, eps);
         CHECK_CLOSE(thloMat(1,0), 566.2182257407396, eps);
         CHECK_CLOSE(thloMat(1,1), 10093.33785879814, eps);
	 CHECK_CLOSE(thlomh2, 2307049.0303978394, eps);
         break;
      case Hierarchies::h6b2qg2:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702311, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832759, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1089.201418061661, eps);
         CHECK_CLOSE(twloMat(1,0),  7145.267026465748, eps);
         CHECK_CLOSE(twloMat(1,1), -2077.345120153528, eps);

         CHECK_CLOSE(thloMat(0,0), 285.3154791763894, eps);
         CHECK_CLOSE(thloMat(1,0), 544.3654284413091, eps);
         CHECK_CLOSE(thloMat(1,1), 10336.22756889787, eps);
	 CHECK_CLOSE(thlomh2, 2335302.18014229, eps);
         break;
      case Hierarchies::h6bq22g:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702315, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832763, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1078.578574572311, eps);
         CHECK_CLOSE(twloMat(1,0),  7096.529601647042, eps);
         CHECK_CLOSE(twloMat(1,1), -1900.197036824461, eps);

         CHECK_CLOSE(thloMat(0,0), 283.0220052455883, eps);
         CHECK_CLOSE(thloMat(1,0), 566.2190953470737, eps);
         CHECK_CLOSE(thloMat(1,1), 10093.33986048966, eps);
	 CHECK_CLOSE(thlomh2, 2307049.804974461, eps);
         break;
      case Hierarchies::h6bq2g2:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702315, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832759, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1089.201418061661, eps);
         CHECK_CLOSE(twloMat(1,0),  7145.267026465748, eps);
         CHECK_CLOSE(twloMat(1,1), -2077.345120153528, eps);

         CHECK_CLOSE(thloMat(0,0), 285.3120881213721, eps);
         CHECK_CLOSE(thloMat(1,0), 544.3662758149513, eps);
         CHECK_CLOSE(thloMat(1,1), 10336.23012077387, eps);
	 CHECK_CLOSE(thlomh2, 2335302.968936102, eps);
         break;
      case Hierarchies::h6g2:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702315, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832761, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1089.201418061661, eps);
         CHECK_CLOSE(twloMat(1,0),  7145.267026465748, eps);
         CHECK_CLOSE(twloMat(1,1), -2112.642999123034, eps);

         CHECK_CLOSE(thloMat(0,0), 246.0217489966267, eps);
         CHECK_CLOSE(thloMat(1,0), 557.451210096066 , eps);
         CHECK_CLOSE(thloMat(1,1), 8628.076480526881, eps);
	 CHECK_CLOSE(thlomh2, 3604063.952735181, eps);
         break;
      case Hierarchies::h9:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.352110199906 , eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0),  420.2050380976995, eps);
         CHECK_CLOSE(twloMat(1,0), -554.6021924866435, eps);
         CHECK_CLOSE(twloMat(1,1), -797.8089039452509, eps);

         CHECK_CLOSE(thloMat(0,0),  132.8584579769461, eps);
         CHECK_CLOSE(thloMat(1,0), -171.9326869339159, eps);
         CHECK_CLOSE(thloMat(1,1), -800.8408283898472, eps);
	 CHECK_CLOSE(thlomh2, -75580.1241347659, eps);
         break;
      case Hierarchies::h9q2:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.352110199906 , eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0),  420.2050380976995, eps);
         CHECK_CLOSE(twloMat(1,0), -554.6021924866435, eps);
         CHECK_CLOSE(twloMat(1,1), -797.8089039452509, eps);

         CHECK_CLOSE(thloMat(0,0),  132.6358855624267, eps);
         CHECK_CLOSE(thloMat(1,0), -171.4711818838455, eps);
         CHECK_CLOSE(thloMat(1,1), -800.9569014303727, eps);
	 CHECK_CLOSE(thlomh2, -74922.1023601382, eps);
         break;
      default:
         REQUIRE_FALSE_MESSAGE(false, "unknown hierarchy!");
         break;
      }
   }
}
