#pragma once

#include "Himalaya_interface.hpp"
#include "dilog.h"
#include <cmath>
#include <limits>

namespace himalaya{
namespace mh2_eft{
   
   /**
    * The Mh2 EFT calculator class
    */
   class Mh2EFTCalculator{
   public:
      /**
       * 	@param at the top yukawa coupling at = yt^2 Sin[beta]^2/(4 Pi).
       * 	@param mt running dr-bar top mass.
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param MR2 the squared renormalization scale.
       * 	@return the squared CP even light Higgs mass at 1-loop level (at) EFT
       */
      double Mh2_EFT_1loop(double at, double mt, double mQ32, double mU32,
			   double Xt, double MR2);
      
      /**
       * 	@param at the top yukawa coupling at = yt^2 Sin[beta]^2/(4 Pi).
       * 	@param mt running dr-bar top mass.
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param MR2 the squared renormalization scale.
       * 	@param g3 the strong coupling constant.
       * 	@param m3 the gluino mass.
       * 	@return the squared CP even light Higgs mass at 2-loop level (at*as) EFT
       */
      double Mh2_EFT_2loop(double at, double mt, double mQ32, double mU32,
			   double Xt, double MR2, double g3, double m3);
      
      /**
       * 	@param at the top yukawa coupling at = yt^2 Sin[beta]^2/(4 Pi).
       * 	@param mt running dr-bar top mass.
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param MR2 the squared renormalization scale.
       * 	@param g3 the strong coupling constant.
       * 	@param m3 the gluino mass.
       * 	@param msq2 the squared common squark mass of the first two generation.
       * 	@return the squared CP even light Higgs mass at 3-loop level (at*as^2) EFT
       */
      double Mh2_EFT_3loop(double at, double mt, double mQ32, double mU32,
			   double Xt, double MR2, double g3, double m3, double msq2);
      
      /**
       * 	@param mQ3 the left-handed soft-breaking mass term of the third generation.
       * 	@param mU3 the right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix.
       * 	@param m3 the gluino mass.
       * 	@param msq the common squark mass of the first two generations.
       * 	@return the 3-loop coefficient O(at*as^2*log^0) with susy logs Log(MR^2/mQ3^2)^0, Log(MR^2/mt^2) -> 0
       */
      double coeff_as2_susy_log0(double mQ3, double mU3, double Xt, double m3, double msq);
      
      /**
       * 	@param mQ3 the left-handed soft-breaking mass term of the third generation.
       * 	@param mU3 the right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix.
       * 	@param m3 the gluino mass.
       * 	@param msq the common squark mass of the first two generations.
       * 	@return the 3-loop coefficient O(at*as^2*log^1) with susy logs Log(MR^2/mQ3^2)^1
       */
      double coeff_as2_susy_log1(double mQ3, double mU3, double Xt, double m3, double msq);
      
      /**
       * 	@param mQ3 the left-handed soft-breaking mass term of the third generation.
       * 	@param mU3 the right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix.
       * 	@param m3 the gluino mass.
       * 	@param msq the common squark mass of the first two generations.
       * 	@return the 3-loop coefficient O(at*as^2*log^2) with susy logs Log(MR^2/mQ3^2)^2
       */
      double coeff_as2_susy_log2(double mQ3, double mU3, double Xt, double m3, double msq);
      
      /**
       * 	@return the 3-loop coefficient O(at*as^2*log^3) with susy logs Log(MR^2/mQ3^2)^3
       */
      double coeff_as2_susy_log3();
      double coeff_as2_susy_log02(double mQ3, double mU3, double m3, double msq, double Mst1, double Xt);
      double coeff_as2_susy_log12(double mQ3, double mU3, double m3, double msq, double Mst1, double Xt);
      double coeff_as2_susy_log22(double mQ3, double mU3, double m3, double msq, double Mst1, double Xt);
      double coeff_as2_susy_log32();
      /**
       *   Returns the matching relation of zeta_lambda^(2) for the degenerated mass case
       *   @param scale the renormalization scale
       *   @param mst1 the mass of the light stop quark
       *   @return zeta_lambda^(2)
       */
     double getZetaDegenerated(double scale, double mst1, double Xt) const;
   private:
      /**
       * 	fin[] function from arXiv:hep-ph/0507139 .
       *
       * 	@param m12 squared mass \f$m_1^2\f$
       * 	@param m22 squared mass \f$m_2^2\f$
       * 	@param MR2 squared renormalization scale
       *
       * 	@return fin(m12, m22)
       */
      double fin(double m12, double m22, double MR2);
      
      /**
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param MR2 the squared renormalization scale.
       * 	@returns the coefficient of order as^0 Log(MR^2/mt^2)^0
       */
      double coeff_as_0_log_0(double mQ32, double mU32, double Xt, double MR2);
      
      /**
       * 	@returns the coefficient of order as^0 Log(MR^2/mt^2)^1
       */
      double coeff_as_0_log_1();
      
      /**
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param m3 the gluino mass.
       * 	@param MR2 the squared renormalization scale.
       * 	@returns the coefficient of order as^1 Log(MR^2/mt^2)^0
       */
      double coeff_as_1_log_0(double mQ32, double mU32, double Xt, double m3, double MR2);
      
      /**
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param m3 the gluino mass.
       * 	@param MR2 the squared renormalization scale.
       * 	@returns the coefficient of order as^1 Log(MR^2/mt^2)^1
       */
      double coeff_as_1_log_1(double mQ32, double mU32, double Xt, double m3, double MR2);
      
      /**
       * 	@returns the coefficient of order as^1 Log(MR^2/mt^2)^2
       */
      double coeff_as_1_log_2();
      
      /**
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param m3 the gluino mass.
       * 	@param msq2 the squared common squark mass of the first two generation.
       * 	@param MR2 the squared renormalization scale.
       * 	@returns the coefficient of order as^2 Log(MR^2/mt^2)^0
       */
      double coeff_as_2_log_0(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2);
      
      /**
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param m3 the gluino mass.
       * 	@param msq2 the squared common squark mass of the first two generation.
       * 	@param MR2 the squared renormalization scale.
       * 	@returns the coefficient of order as^2 Log(MR^2/mt^2)^1
       */
      double coeff_as_2_log_1(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2);
      
      /**
       * 	@param mQ32 the squared left-handed soft-breaking mass term of the third generation.
       * 	@param mU32 the squared right-handed soft-breaking mass term of the third generation.
       * 	@param Xt the Xt parameter of the stop mixing matrix (At - mu*cotb).
       * 	@param m3 the gluino mass.
       * 	@param msq2 the squared common squark mass of the first two generation.
       * 	@param MR2 the squared renormalization scale.
       * 	@returns the coefficient of order as^2 Log(MR^2/mt^2)^2
       */
      double coeff_as_2_log_2(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2);
      
      /**
       * 	@returns the coefficient of order as^2 Log(MR^2/mt^2)^3
       */
      double coeff_as_2_log_3();
      std::string tf(const bool tf);
   };
}	// mh2_eft
}	// himalaya
