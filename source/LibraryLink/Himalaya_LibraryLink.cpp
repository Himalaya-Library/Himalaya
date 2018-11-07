// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "HierarchyCalculator.hpp"
#include "Mh2EFTCalculator.hpp"
#include "MSSM_mass_eigenstates.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>

#include <mathlink.h>
#include <WolframLibrary.h>

namespace himalaya {
namespace {

/********************* put types *********************/

inline void MLPut(MLINK link, const std::string& s)
{
   MLPutSymbol(link, s.c_str());
}

inline void MLPut(MLINK link, int c)
{
   MLPutInteger(link, c);
}

inline void MLPut(MLINK link, double c)
{
   MLPutReal(link, c);
}

inline void MLPut(MLINK link, std::complex<double> c)
{
   if (std::imag(c) == 0.) {
      MLPutReal(link, std::real(c));
   } else {
      MLPutFunction(link, "Complex", 2);
      MLPutReal(link, std::real(c));
      MLPutReal(link, std::imag(c));
   }
}

template <int M>
void MLPut(MLINK link, const Eigen::Array<double,M,1>& a)
{
   double v[M];
   for (int i = 0; i < M; i++)
      v[i] = a(i);
   MLPutRealList(link, v, M);
}

template <int M>
void MLPut(MLINK link, const Eigen::Matrix<double,M,1>& m)
{
   const Eigen::Array<double,M,1> a(m.array());
   MLPut(link, a);
}

template <int M, int N>
void MLPut(MLINK link, const Eigen::Matrix<double,M,N>& m)
{
   double mat[M][N];
   for (int i = 0; i < M; i++)
      for (int k = 0; k < N; k++)
         mat[i][k] = m(i, k);

   long dims[] = { M, N };
   MLPutDoubleArray(link, (double*)mat, dims, NULL, 2);
}

template <int M>
void MLPut(MLINK link, const Eigen::Array<std::complex<double>,M,1>& a)
{
   MLPutFunction(link, "List", M);
   for (int i = 0; i < M; i++)
      MLPut(link, a(i));
}

template <int M>
void MLPut(MLINK link, const Eigen::Matrix<std::complex<double>,M,1>& m)
{
   const Eigen::Array<std::complex<double>,M,1> a(m.array());
   MLPut(link, a);
}

template <int M, int N>
void MLPut(MLINK link, const Eigen::Matrix<std::complex<double>,M,N>& m)
{
   MLPutFunction(link, "List", M);
   for (int i = 0; i < M; i++) {
      MLPutFunction(link, "List", N);
      for (int k = 0; k < N; k++)
         MLPut(link, m(i,k));
   }
}

template <class T>
void MLPut(MLINK link, const std::vector<T>& v)
{
   MLPutFunction(link, "List", v.size());
   for (std::size_t i = 0; i < v.size(); i++)
      MLPut(link, v[i]);
}

/********************* put rules to types *********************/

void MLPutRule(MLINK link, const std::string& name)
{
   MLPutFunction(link, "Rule", 2);
   MLPutUTF8Symbol(link, reinterpret_cast<const unsigned char*>(name.c_str()), name.size());
}

template <class T1, class T2>
void MLPutRuleTo(MLINK link, T1 t, const T2& name)
{
   MLPutRule(link, name);
   MLPut(link, t);
}

/******************************************************************/

void put_message(MLINK link,
                 const std::string& message_function,
                 const std::string& message_str)
{
   MLPutFunction(link, "CompoundExpression", 2);
   MLPutFunction(link, message_function.c_str(), 1);
   MLPutUTF8String(link, reinterpret_cast<const unsigned char*>(message_str.c_str()), message_str.size());
}

/******************************************************************/

class Redirect_output {
public:
   explicit Redirect_output(MLINK link_)
      : link(link_)
      , buffer()
      , old_cout(std::cout.rdbuf(buffer.rdbuf()))
      , old_cerr(std::cerr.rdbuf(buffer.rdbuf()))
      {}

   ~Redirect_output() {
      std::cout.rdbuf(old_cout);
      std::cerr.rdbuf(old_cerr);
      flush();
   }

   void flush() {
      std::string line;
      while (std::getline(buffer, line)) {
         put_message(link, "HimalayaInfoMessage", line);
      }
   }

private:
   MLINK link;               ///< redirect to this link
   std::stringstream buffer; ///< buffer caching stdout
   std::streambuf* old_cout; ///< original stdout buffer
   std::streambuf* old_cerr; ///< original stderr buffer
};

/******************************************************************/

long number_of_args(MLINK link, const std::string& head)
{
   long argc;

   if (!MLCheckFunction(link, head.c_str(), &argc))
      std::cerr << "Error: argument is not a " << head << std::endl;

   return argc;
}

/******************************************************************/

bool check_number_of_args(MLINK link, long number_of_arguments,
                          const std::string& function_name)
{
   const auto n_given = number_of_args(link, "List");
   const bool ok = n_given == number_of_arguments;

   if (!ok) {
      std::cerr << "Error: " << function_name << " expects "
                << number_of_arguments << " argument ("
                << n_given << " given)." << std::endl;
   }

   return ok;
}

/******************************************************************/

std::vector<double> read_list(MLINK link)
{
   int N = 0;

   if (!MLTestHead(link, "List", &N)) {
      throw std::runtime_error("HimalayaCalculateDMh3L expects a list"
                               " as the only argument!");
   }

   std::vector<double> vec(N, 0.);

   for (int i = 0; i < N; i++) {
      double val = 0.;
      if (!MLGetReal64(link, &val)) {
         throw std::runtime_error("Cannot read " + std::to_string(i)
                                  + "'th value of parameter list!");
      }
      vec[i] = val;
   }

   if (!MLNewPacket(link))
      throw std::runtime_error("Cannot create new packet!");

   return vec;
}

/******************************************************************/

struct Data {
   Data(const himalaya::Parameters& pars_, bool bottom_, bool verbose_)
      : pars(pars_)
      , bottom(bottom_)
      , verbose(verbose_)
      {}

   himalaya::Parameters pars{};
   bool bottom{false};
   bool verbose{true};
};

/******************************************************************/

Data make_data(const std::vector<double>& parsvec)
{
   const int N_input_parameters = 123; // number of Himalaya input parameters

   if (parsvec.size() != N_input_parameters) {
      throw std::runtime_error("HimalayaCalculateDMh3L expects "
                               + std::to_string(N_input_parameters) + ", but "
                               + std::to_string(parsvec.size()) + " given!");
   }

   int c = 0; // counter

   const bool bottom = parsvec.at(c++);
   const bool verbose = parsvec.at(c++);

   himalaya::Parameters pars;
   pars.scale = parsvec.at(c++);
   pars.mu = parsvec.at(c++);
   pars.g1 = parsvec.at(c++);
   pars.g2 = parsvec.at(c++);
   pars.g3 = parsvec.at(c++);
   pars.vd = parsvec.at(c++);
   pars.vu = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.mq2(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.md2(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.mu2(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.ml2(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.me2(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.Au(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.Ad(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.Ae(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.Yu(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.Yd(i,k) = parsvec.at(c++);
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         pars.Ye(i,k) = parsvec.at(c++);
   pars.MA = parsvec.at(c++);
   pars.M1 = parsvec.at(c++);
   pars.M2 = parsvec.at(c++);
   pars.MG = parsvec.at(c++);

   const double MW = parsvec.at(c++);
   const double MZ = parsvec.at(c++);
   const double Mt = parsvec.at(c++);
   const double Mb = parsvec.at(c++);
   const double Mtau = parsvec.at(c++);

   if (MW > 0) pars.MW = MW;
   if (MZ > 0) pars.MZ = MZ;
   if (Mt > 0) pars.Mt = Mt;
   if (Mb > 0) pars.Mb = Mb;
   if (Mtau > 0) pars.Mtau = Mtau;

   Eigen::Vector2d MSt, MSb;
   MSt(0) = parsvec.at(c++);
   MSt(1) = parsvec.at(c++);
   MSb(0) = parsvec.at(c++);
   MSb(1) = parsvec.at(c++);
   const double s2t = parsvec.at(c++);
   const double s2b = parsvec.at(c++);

   if (MSt.minCoeff() > 0. && std::abs(s2t) <= 1.) {
      pars.MSt = MSt;
      pars.s2t = s2t;
   }

   if (MSb.minCoeff() > 0. && std::abs(s2b) <= 1.) {
      pars.MSb = MSb;
      pars.s2b = s2b;
   }

   pars.validate(verbose);

   if (c != N_input_parameters) {
      throw std::runtime_error(
         "Bug: Expecting to read " + std::to_string(N_input_parameters) +
         " input parameters from input vector of size " +
         std::to_string(parsvec.size()) + ", but only " + std::to_string(c) +
         " parameters have been read.");
   }

   return Data(pars, bottom, verbose);
}

/******************************************************************/

struct Results {
   using Loop_corrections = std::tuple<double,double,double,double>;

   himalaya::HierarchyObject ho{false};
   Loop_corrections eft; ///< fixed-order corrections for v^2 << MS^2
   Loop_corrections fo;  ///< fixed-order corrections
};

/******************************************************************/

himalaya::Parameters make_gaugeless(const himalaya::Parameters& pars)
{
   auto gl = pars;
   gl.g1 = 0.;
   gl.g2 = 0.;
   gl.MW = himalaya::NaN;
   gl.MZ = himalaya::NaN;
   gl.validate(false);
   return gl;
}

/******************************************************************/

Results calculate_results(const Data& data)
{
   Results res;

   himalaya::HierarchyCalculator hc(data.pars, data.verbose);
   res.ho = hc.calculateDMh3L(data.bottom);

   // create point with g1 = g2 = 0
   const auto pars_gl = make_gaugeless(data.pars);

   // calculate fixed-order corrections for v^2 << MS^2
   himalaya::mh2_eft::Mh2EFTCalculator meft(pars_gl);
   const auto dmh2_eft_0l = meft.getDeltaMh2EFT0Loop();
   const auto dmh2_eft_1l = meft.getDeltaMh2EFT1Loop(1,1);
   const auto dmh2_eft_2l = meft.getDeltaMh2EFT2Loop(1,1);

   res.eft = std::make_tuple(dmh2_eft_0l, dmh2_eft_1l,
                             dmh2_eft_2l, res.ho.getDMh2EFT(3));

   // calculate fixed-order corrections
   himalaya::mh2_fo::MSSM_mass_eigenstates mfo(pars_gl);
   const auto dmh_fo    = mfo.calculate_Mh2(); // 0L, 1L, 2L
   const auto dmh_fo_3l = res.ho.getDMh2(3);   // 3L

   res.fo = std::tuple_cat(dmh_fo, std::tie(dmh_fo_3l));

   return res;
}

/******************************************************************/

void put_result(const Results& res, MLINK link)
{
   MLPutFunction(link, "List", 12);

   const auto& ho = res.ho;
   const auto& eft = res.eft;
   const auto& fo = res.fo;

   const auto hierarchy = ho.getSuitableHierarchy();
   const std::string msf = ho.getIsAlphab() ? "MsbottomMDRPrime" : "MstopMDRPrime";

   Eigen::Vector4d expansion_uncertainty;
   expansion_uncertainty << 0., ho.getDMhExpUncertainty(1),
      ho.getDMhExpUncertainty(2), ho.getDMhExpUncertainty(3);

   Eigen::Vector4d lambda;
   lambda << ho.getDLambda(0), ho.getDLambda(1),
             ho.getDLambda(2), ho.getDLambda(3);

   Eigen::Vector4d lambda_uncertainty;
   lambda_uncertainty << 0., 0., 0., ho.getDLambdaUncertainty(3);

   Eigen::Vector4d lambda_shift_DRp_to_MS;
   lambda_shift_DRp_to_MS <<
      ho.getDLambdaDRbarPrimeToMSbarShift(0),
      ho.getDLambdaDRbarPrimeToMSbarShift(1),
      ho.getDLambdaDRbarPrimeToMSbarShift(2),
      ho.getDLambdaDRbarPrimeToMSbarShift(3);

   std::vector<Eigen::Matrix2d> Mh2;
   for (int i = 0; i < 4; i++)
      Mh2.push_back(ho.getDMh(i));

   std::vector<Eigen::Matrix2d> Mh2_shift_DRp_to_MDRp;
   Mh2_shift_DRp_to_MDRp.push_back(Eigen::Matrix2d::Zero());
   Mh2_shift_DRp_to_MDRp.push_back(Eigen::Matrix2d::Zero());
   Mh2_shift_DRp_to_MDRp.push_back(Eigen::Matrix2d::Zero());
   Mh2_shift_DRp_to_MDRp.push_back(ho.getDMhDRbarPrimeToMDRbarPrimeShift());

   std::vector<Eigen::Matrix2d> Mh2_shift_DRp_to_H3m;
   Mh2_shift_DRp_to_H3m.push_back(Eigen::Matrix2d::Zero());
   Mh2_shift_DRp_to_H3m.push_back(Eigen::Matrix2d::Zero());
   Mh2_shift_DRp_to_H3m.push_back(Eigen::Matrix2d::Zero());
   Mh2_shift_DRp_to_H3m.push_back(ho.getDMhDRbarPrimeToH3mShift());

   Eigen::Vector4d Mh2_eft;
   Mh2_eft << std::get<0>(eft), std::get<1>(eft),
              std::get<2>(eft), std::get<3>(eft);

   Eigen::Vector4d Mh2_fo;
   Mh2_fo << std::get<0>(fo), std::get<1>(fo),
             std::get<2>(fo), std::get<3>(fo);

   MLPutRuleTo(link, hierarchy, "hierarchyID");
   MLPutRuleTo(link, ho.getH3mHierarchyNotation(hierarchy), "hierarchyName");
   MLPutRuleTo(link, ho.getMDRMasses(), msf);
   MLPutRuleTo(link, Mh2, "Mh2");
   MLPutRuleTo(link, Mh2_shift_DRp_to_MDRp, "Mh2ShiftDRbarPrimeToMDRPrime");
   MLPutRuleTo(link, Mh2_shift_DRp_to_H3m, "Mh2ShiftDRbarPrimeToH3m");
   MLPutRuleTo(link, expansion_uncertainty, "expansionUncertainty");
   MLPutRuleTo(link, Mh2_eft, "Mh2EFT");
   MLPutRuleTo(link, Mh2_fo, "Mh2FO");
   MLPutRuleTo(link, lambda, "lambda");
   MLPutRuleTo(link, lambda_uncertainty, "lambdaUncertainty");
   MLPutRuleTo(link, lambda_shift_DRp_to_MS, "lambdaShiftDRbarPrimeToMSbar");

   MLEndPacket(link);
}

} // anonymous namespace
} // namespace himalaya

extern "C" {

/******************************************************************/

DLLEXPORT int HimalayaCalculateDMh3L(
   WolframLibraryData /* libData */, MLINK link)
{
   using namespace himalaya;

   if (!check_number_of_args(link, 1, "HimalayaCalculateDMh3L"))
      return LIBRARY_TYPE_ERROR;

   try {
      Redirect_output rd(link);

      const auto res = calculate_results(make_data(read_list(link)));

      rd.flush();

      put_result(res, link);
   } catch (const std::exception& e) {
      put_message(link, "HimalayaErrorMessage", e.what());
      MLPutSymbol(link, "$Failed");
   } catch (...) {
      put_message(link, "HimalayaErrorMessage", "An unknown exception has been thrown.");
      MLPutSymbol(link, "$Failed");
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT mint WolframLibrary_getVersion()
{
   return WolframLibraryVersion;
}

/******************************************************************/

DLLEXPORT int WolframLibrary_initialize(WolframLibraryData /* libData */)
{
   return LIBRARY_NO_ERROR;
}

} // extern "C"
