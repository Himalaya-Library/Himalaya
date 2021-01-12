#include "doctest.h"

#include "himalaya/HierarchyCalculator.hpp"

#include "Enums.hpp"
#include "EFTFlags.hpp"

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
   pars.Ad << 0, 0, 0, 0, 0, 0, 0, 0,  9996.81;
   pars.Au << 0, 0, 0, 0, 0, 0, 0, 0,  6992.34;

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
   pars.theta_t = -0.783817;
   pars.theta_b = -0.291498;

   return pars;
}


/**
 * Performs a sanity check of the implemented expansion terms by
 * comparing them to their numerical value at a given parameter point.
 */
TEST_CASE("test_FO_expansions")
{
   using namespace himalaya;
   using namespace himalaya::mh2_eft;
   using namespace himalaya::hierarchies;

   const double eps = 1e-10;

   for (int i = 0; i < Hierarchies::NUMBER_OF_HIERARCHIES; i++) {
      himalaya::HierarchyCalculator hc(make_point());
      himalaya::HierarchyObject ho(false);
      ho.setRenormalizationScheme(RenSchemes::H3m);
      ho.setMDRFlag(1);
      ho.setSuitableHierarchy(i);

      const auto oloMat  = hc.calculateHierarchy(ho, 1, 0, 0);
      const auto twloMat = hc.calculateHierarchy(ho, 0, 1, 0);
      const auto thloMat = hc.calculateHierarchy(ho, 0, 0, 1);
      ho.setMDRFlag(0);
      hc.calculateHierarchy(ho, 0, 0, 1);
      const auto thlomh2 = ho.getDLambdaH3m();
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

         CHECK_CLOSE(thloMat(0,0), 0.9259290788, 2e-6);
         CHECK_CLOSE(thloMat(1,0), 10.2799042295, 2e-6);
         CHECK_CLOSE(thloMat(1,1), 370.7163915042, 2e-6);
	 CHECK_CLOSE(thlomh2, 7971.492926307, eps);
	 break;
      case Hierarchies::h32q2g:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.3521101999062, eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0), -13.66052379180129, eps);
         CHECK_CLOSE(twloMat(1,0),  11.26755617866339, eps);
         CHECK_CLOSE(twloMat(1,1),  1477.465656153518, eps);

         CHECK_CLOSE(thloMat(0,0), 1.1133090886, eps);
         CHECK_CLOSE(thloMat(1,0), 9.9034479467, eps);
         CHECK_CLOSE(thloMat(1,1), 370.8935498384, eps);
	 CHECK_CLOSE(thlomh2, 7902.0277533642, eps);
         break;
      case Hierarchies::h3q22g:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.3521101999062, eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0), -13.66052379180129, eps);
         CHECK_CLOSE(twloMat(1,0),  11.26755617866339, eps);
         CHECK_CLOSE(twloMat(1,1),  1477.465656153518, eps);

         CHECK_CLOSE(thloMat(0,0), 0.8877658852, eps);
         CHECK_CLOSE(thloMat(1,0), 10.3072252597, eps);
         CHECK_CLOSE(thloMat(1,1), 370.7968728561, eps);
	 CHECK_CLOSE(thlomh2, 8000.0243821962, eps);
         break;
      case Hierarchies::h4:
         CHECK_CLOSE(oloMat(0,0),                 0, eps);
         CHECK_CLOSE(oloMat(1,0),                 0, eps);
         CHECK_CLOSE(oloMat(1,1), 6685.123085628641, eps);

         CHECK_CLOSE(twloMat(0,0),                 0, eps);
         CHECK_CLOSE(twloMat(1,0), 1183.325484493686, eps);
         CHECK_CLOSE(twloMat(1,1), 1458.970501474495, eps);

         CHECK_CLOSE(thloMat(0,0), 162.1379211675, eps);
         CHECK_CLOSE(thloMat(1,0), 326.0240309275, eps);
         CHECK_CLOSE(thloMat(1,1), 431.6924771365, eps);
	 CHECK_CLOSE(thlomh2, -1422.7254877425, eps);
         break;
      case Hierarchies::h5:
         CHECK_CLOSE(oloMat(0,0),  15921.69462848581, eps);
         CHECK_CLOSE(oloMat(1,0), -388569.2043081555, eps);
         CHECK_CLOSE(oloMat(1,1),  7874.401574063407, eps);

         CHECK_CLOSE(twloMat(0,0), -86.77887344841422, eps);
         CHECK_CLOSE(twloMat(1,0), -20625.63783863484, eps);
         CHECK_CLOSE(twloMat(1,1), -42446.62009872038, eps);

         CHECK_CLOSE(thloMat(0,0),  2442.1147633261, eps);
         CHECK_CLOSE(thloMat(1,0), -3859.9417043341, eps);
         CHECK_CLOSE(thloMat(1,1),  60592.9761507541  , eps);
	 CHECK_CLOSE(thlomh2, 11229598.8300867677, eps);
         break;
      case Hierarchies::h5g1:
         CHECK_CLOSE(oloMat(0,0),  15921.69462848581, eps);
         CHECK_CLOSE(oloMat(1,0), -388569.2043081556, eps);
         CHECK_CLOSE(oloMat(1,1),  7874.401574063407, eps);

         CHECK_CLOSE(twloMat(0,0), -114.6037388932203, eps);
         CHECK_CLOSE(twloMat(1,0), -20341.84471909946, eps);
         CHECK_CLOSE(twloMat(1,1), -42843.48046642416, eps);

         CHECK_CLOSE(thloMat(0,0),  2415.5071508754, eps);
         CHECK_CLOSE(thloMat(1,0), -3766.7487013104, eps);
         CHECK_CLOSE(thloMat(1,1),  59380.2606366422, eps);
	 CHECK_CLOSE(thlomh2, 11936958.8562686481, eps);
         break;
      case Hierarchies::h6:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702315, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832763, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1078.578574572312, eps);
         CHECK_CLOSE(twloMat(1,0),  7096.529601647042, eps);
         CHECK_CLOSE(twloMat(1,1), -1927.791631086123, eps);

         CHECK_CLOSE(thloMat(0,0), 245.44331586, eps);
         CHECK_CLOSE(thloMat(1,0), 573.1287426057, eps);
         CHECK_CLOSE(thloMat(1,1), 8448.4603822038  , eps);
	 CHECK_CLOSE(thlomh2, 3575575.4157333178, eps);
         break;
      case Hierarchies::h6b:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702311, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832763, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1078.578574572312, eps);
         CHECK_CLOSE(twloMat(1,0),  7096.52960164704 , eps);
         CHECK_CLOSE(twloMat(1,1), -1900.197036824461, eps);

         CHECK_CLOSE(thloMat(0,0), 283.0253419222, eps);
         CHECK_CLOSE(thloMat(1,0), 566.2183368301, eps);
         CHECK_CLOSE(thloMat(1,1), 10093.3364147391, eps);
	 CHECK_CLOSE(thlomh2, 2307049.0303978394, eps);
         break;
      case Hierarchies::h6b2qg2:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702311, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832759, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1089.201418061661, eps);
         CHECK_CLOSE(twloMat(1,0),  7145.267026465748, eps);
         CHECK_CLOSE(twloMat(1,1), -2077.345120153528, eps);

         CHECK_CLOSE(thloMat(0,0), 285.3154430301, eps);
         CHECK_CLOSE(thloMat(1,0), 544.3655499119, eps);
         CHECK_CLOSE(thloMat(1,1), 10336.2260013211, eps);
	 CHECK_CLOSE(thlomh2, 2335302.1799157225, eps);
         break;
      case Hierarchies::h6bq22g:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702315, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832763, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1078.578574572311, eps);
         CHECK_CLOSE(twloMat(1,0),  7096.529601647042, eps);
         CHECK_CLOSE(twloMat(1,1), -1900.197036824461, eps);

         CHECK_CLOSE(thloMat(0,0), 283.021970571, eps);
         CHECK_CLOSE(thloMat(1,0), 566.2192063187, eps);
         CHECK_CLOSE(thloMat(1,1), 10093.3384161623, eps);
	 CHECK_CLOSE(thlomh2, 2307049.8046336095, eps);
         break;
      case Hierarchies::h6bq2g2:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702315, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832759, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1089.201418061661, eps);
         CHECK_CLOSE(twloMat(1,0),  7145.267026465748, eps);
         CHECK_CLOSE(twloMat(1,1), -2077.345120153528, eps);

         CHECK_CLOSE(thloMat(0,0), 285.3120524328, eps);
         CHECK_CLOSE(thloMat(1,0), 544.3663971707, eps);
         CHECK_CLOSE(thloMat(1,1), 10336.2285528529, eps);
	 CHECK_CLOSE(thlomh2, 2335302.9685890833, eps);
         break;
      case Hierarchies::h6g2:
         CHECK_CLOSE(oloMat(0,0),  9272.477351702315, eps);
         CHECK_CLOSE(oloMat(1,0), -184.7601614832761, eps);
         CHECK_CLOSE(oloMat(1,1),  7581.278122072418, eps);

         CHECK_CLOSE(twloMat(0,0), -1089.201418061661, eps);
         CHECK_CLOSE(twloMat(1,0),  7145.267026465748, eps);
         CHECK_CLOSE(twloMat(1,1), -2112.642999123034, eps);

         CHECK_CLOSE(thloMat(0,0), 246.0239126332, eps);
         CHECK_CLOSE(thloMat(1,0), 557.450061436 , eps);
         CHECK_CLOSE(thloMat(1,1), 8628.0794608788, eps);
	 CHECK_CLOSE(thlomh2, 3604058.4466279238, eps);
         break;
      case Hierarchies::h9:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.352110199906 , eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0),  420.2050380976995, eps);
         CHECK_CLOSE(twloMat(1,0), -554.6021924866435, eps);
         CHECK_CLOSE(twloMat(1,1), -797.8089039452509, eps);

         CHECK_CLOSE(thloMat(0,0),  132.8596207671, eps);
         CHECK_CLOSE(thloMat(1,0), -171.9677397533, eps);
         CHECK_CLOSE(thloMat(1,1), -799.518525656, eps);
	 CHECK_CLOSE(thlomh2, -75348.7029377545, eps);
         break;
      case Hierarchies::h9q2:
         CHECK_CLOSE(oloMat(0,0), -1033.437882123761, eps);
         CHECK_CLOSE(oloMat(1,0), -394.352110199906 , eps);
         CHECK_CLOSE(oloMat(1,1),  17633.47392819223, eps);

         CHECK_CLOSE(twloMat(0,0),  420.2050380976995, eps);
         CHECK_CLOSE(twloMat(1,0), -554.6021924866435, eps);
         CHECK_CLOSE(twloMat(1,1), -797.8089039452509, eps);

         CHECK_CLOSE(thloMat(0,0),  131.9045772644, eps);
         CHECK_CLOSE(thloMat(1,0), -170.414374683, eps);
         CHECK_CLOSE(thloMat(1,1), -799.1630905316, eps);
	 CHECK_CLOSE(thlomh2, -73583.5048616582, eps);
         break;
      default:
         REQUIRE_FALSE_MESSAGE(false, "unknown hierarchy!");
         break;
      }
   }
}
