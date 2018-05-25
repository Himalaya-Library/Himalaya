#include "HierarchyCalculator.hpp"

#include <iostream>
#include <string>
#include <sstream>

#include <mathlink.h>
#include <WolframLibrary.h>

namespace himalaya {

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

/********************* put single heads *********************/

void MLPutHeads(MLINK link, const std::vector<std::string>& heads)
{
   for (const auto& h: heads)
      MLPutFunction(link, h.c_str(), 1);
}

/********************* put rules to types *********************/

void MLPutRule(MLINK link, const std::string& name, const std::vector<std::string>& heads = {})
{
   MLPutFunction(link, "Rule", 2);
   MLPutHeads(link, heads);
   MLPutUTF8Symbol(link, reinterpret_cast<const unsigned char*>(name.c_str()), name.size());
}

void MLPutRule(MLINK link, int number, const std::vector<std::string>& heads = {})
{
   MLPutFunction(link, "Rule", 2);
   MLPutHeads(link, heads);
   MLPutInteger(link, number);
}

void MLPutRule(MLINK link, long number, const std::vector<std::string>& heads = {})
{
   MLPutFunction(link, "Rule", 2);
   MLPutHeads(link, heads);
   MLPutLongInteger(link, number);
}

template <class T1, class T2>
void MLPutRuleTo(MLINK link, T1 t, const T2& name, const std::vector<std::string>& heads = {})
{
   MLPutRule(link, name, heads);
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

void put_result(const himalaya::HierarchyObject& ho, MLINK link)
{
   MLPutFunction(link, "List", 15);

   const auto hierarchy = ho.getSuitableHierarchy();
   const auto ren_scheme = ho.getRenormalizationScheme();
   const auto ren_scheme_str =
      ren_scheme == himalaya::RenSchemes::H3m ||
      ren_scheme == himalaya::RenSchemes::H3mMDRBAR ? "H3m" : "DRbarPrime";

   Eigen::Matrix<double,4,1> expansion_uncertainty;
   expansion_uncertainty << 0., ho.getExpUncertainty(1),
      ho.getExpUncertainty(2), ho.getExpUncertainty(3);

   MLPutRuleTo(link, ren_scheme_str, "renormalizationScheme");
   MLPutRuleTo(link, hierarchy, "hierarchyID");
   MLPutRuleTo(link, ho.getH3mHierarchyNotation(hierarchy), "hierarchyName");
   MLPutRuleTo(link, ho.getMDRMasses(), "Mstop");
   MLPutRuleTo(link, ho.getDMh(0), "Mh2Tree");
   MLPutRuleTo(link, ho.getDMh(1), "Mh21Loop");
   MLPutRuleTo(link, ho.getDMh(2), "Mh22Loop");
   MLPutRuleTo(link, ho.getDMh(3), "Mh23Loop");
   MLPutRuleTo(link, expansion_uncertainty, "expansionUncertainty");
   MLPutRuleTo(link, ho.getDeltaLambdaHimalaya(), "deltaLambda3LoopHimalayaDRbarPrime");
   MLPutRuleTo(link, ho.getDRbarPrimeToMSbarShiftHimalaya(), "deltaLambda3LoopHimalayaShiftDRbarPrimeToMSbar");
   MLPutRuleTo(link, ho.getDeltaLambdaUncertaintyHimalaya(), "deltaLambda3LoopHimalayaUncertainty");
   MLPutRuleTo(link, ho.getDeltaLambdaEFT(), "deltaLambda3LoopEFTDRbarPrime");
   MLPutRuleTo(link, ho.getDRbarPrimeToMSbarShiftEFT(), "deltaLambda3LoopEFTShiftDRbarPrimeToMSbar");
   MLPutRuleTo(link, ho.getDeltaLambdaUncertaintyEFT(), "deltaLambda3LoopEFTUncertainty");

   MLEndPacket(link);
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

Data make_data(const std::vector<double>& parsvec)
{
   const int N_input_parameters = 48; // number of Himalaya input parameters

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
   pars.Ab = parsvec.at(c++);
   pars.At = parsvec.at(c++);
   pars.MA = parsvec.at(c++);
   pars.MG = parsvec.at(c++);
   pars.MW = parsvec.at(c++);
   pars.MZ = parsvec.at(c++);
   pars.Mt = parsvec.at(c++);
   pars.Mb = parsvec.at(c++);

   Eigen::Matrix<double,2,1> MSt, MSb;
   double s2t, s2b;
   MSt(0) = parsvec.at(c++);
   MSt(1) = parsvec.at(c++);
   MSb(0) = parsvec.at(c++);
   MSb(1) = parsvec.at(c++);
   s2t = parsvec.at(c++);
   s2b = parsvec.at(c++);

   if (MSt.minCoeff() > 0. && std::abs(s2t) <= 1.) {
      pars.MSt = MSt;
      pars.s2t = s2t;
   }

   if (MSb.minCoeff() > 0. && std::abs(s2b) <= 1.) {
      pars.MSb = MSb;
      pars.s2b = s2b;
   }

   if (c != N_input_parameters) {
      throw std::runtime_error(
         "Bug: Expecting to read " + std::to_string(N_input_parameters) +
         " input parameters from input vector of size " +
         std::to_string(parsvec.size()) + ", but only " + std::to_string(c) +
         " parameters have been read.");
   }

   return Data(pars, bottom, verbose);
}

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

      const auto data = make_data(read_list(link));

      himalaya::HierarchyCalculator hc(data.pars, data.verbose);
      const auto ho = hc.calculateDMh3L(data.bottom);

      rd.flush();

      put_result(ho, link);
   } catch (const std::exception& e) {
      put_message(link, "HimalayaErrorMessage", e.what());
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
