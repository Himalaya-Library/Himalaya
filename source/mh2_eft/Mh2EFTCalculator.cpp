// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "Mh2EFTCalculator.hpp"
#include "ThresholdCalculator.hpp"
#include "EFTFlags.hpp"
#include "Logger.hpp"
#include "Powers.hpp"
#include <cmath>
#include <iostream>
#include <string>

/**
 * @file Mh2EFTCalculator.cpp
 * @brief Implementation of EFT Higgs mass calculation class.
 */

#define CALC_IF(cond,expr) ((cond) ? (expr) : 0)

namespace himalaya
{
namespace mh2_eft
{

namespace
{
const double zt3 = 1.2020569031595942853997381615114; // zeta(3)
const double Pi  = 3.1415926535897932384626433832795;

template <typename T>
bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept {
    return std::fabs(a) < prec;
}

template <typename T>
bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept {
    return is_zero(a - b, prec);
}

template <typename T>
bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept {
    if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
        return true;

    if (std::abs(a) < std::numeric_limits<T>::epsilon() ||
    std::abs(b) < std::numeric_limits<T>::epsilon())
        return false;

    return std::abs((a - b) / a) < prec;
}

/// Returns var if not NaN, 0 otherwise
double isNaN(double var, const std::string& msg = "")
{
    if (std::isnan(var)) {
        WARNING_MSG("NaN appeared in calculation of threshold correction"
                    << (msg.empty() ? "" : " " + msg) << "!");
        return 0.;
    }
    return var;
}

/// Re(B0(s,x,x,q2)), Eq.(2.4) from [hep-ph/0701051]
double fB(double s, double x, double q2)
{
    using std::asin;
    using std::log;
    using std::sqrt;

    if (is_zero(s) && is_zero(x))
        return 0.0;

    if (is_zero(s))
        return -log(x / q2);

    if (is_zero(x))
        return 2.0 - log(s / q2);

    if (is_equal(s, x))
        return 2.0 - 1.813799364234218 - log(x / q2);

    if (s <= 4.0 * x)
        return 2.0 - log(x / q2) - 2.0 * sqrt(4.0 * x / s - 1.0) * asin(sqrt(s / (4.0 * x)));

    const double sq = sqrt(1.0 - 4.0 * x / s);

    // s > 4*x
    return 2.0 - log(x / q2) + sq * log(s * (1.0 - sq) / (2 * x) - 1.0);
}

} // anonymous namespace

/**
 * Constructor
 * @param p_ a HimalayaInterface struct
 * @param msq2_ the averaged squark mass of the first two generations squared
 * @param verbose a bool enable the output of the parameter validation. Enabled by default
 */
Mh2EFTCalculator::Mh2EFTCalculator(
    const himalaya::Parameters& p_, double msq2_, bool verbose)
    : p(p_), msq2(msq2_)
{
    p.validate(verbose);

    if (!std::isfinite(msq2_))
        msq2 = p.calculateMsq2();

    // fill orders
    orders.fill(1);

    const double eps = 1e-10;

    //TODO check consistency
    if (std::abs(p.g1) < eps) {
        setCorrectionFlag(EFTOrders::G12G22, 0);
        setCorrectionFlag(EFTOrders::G12YB2, 0);
        setCorrectionFlag(EFTOrders::G14, 0);
        setCorrectionFlag(EFTOrders::G12YB2, 0);
        setCorrectionFlag(EFTOrders::G12YTAU2, 0);
        setCorrectionFlag(EFTOrders::G12YT2, 0);
    }

    if (std::abs(p.g2) < eps) {
        setCorrectionFlag(EFTOrders::G24, 0);
        setCorrectionFlag(EFTOrders::G22YB2, 0);
        setCorrectionFlag(EFTOrders::G22YTAU2, 0);
        setCorrectionFlag(EFTOrders::G22YT2, 0);
    }

    if (std::abs(p.Mt) < eps) {
        setCorrectionFlag(EFTOrders::G12YT2  , 0);
        setCorrectionFlag(EFTOrders::G22YT2  , 0);
        setCorrectionFlag(EFTOrders::YT4     , 0);
        setCorrectionFlag(EFTOrders::G32YT4  , 0);
        setCorrectionFlag(EFTOrders::YT6     , 0);
        setCorrectionFlag(EFTOrders::YTAU2YT4, 0);
        setCorrectionFlag(EFTOrders::YT2YB4  , 0);
        setCorrectionFlag(EFTOrders::YB2YT4  , 0);
    }

    if (std::abs(p.Mb) < eps) {
        setCorrectionFlag(EFTOrders::G12YB2  , 0);
        setCorrectionFlag(EFTOrders::G22YB2  , 0);
        setCorrectionFlag(EFTOrders::YB4     , 0);
        setCorrectionFlag(EFTOrders::G32YB4  , 0);
        setCorrectionFlag(EFTOrders::YB6     , 0);
        setCorrectionFlag(EFTOrders::YTAU2YB4, 0);
        setCorrectionFlag(EFTOrders::YT2YB4  , 0);
        setCorrectionFlag(EFTOrders::YB2YT4  , 0);
        setCorrectionFlag(EFTOrders::YTAU4YB2, 0);
        setCorrectionFlag(EFTOrders::YT6     , 0);
    }

    if (std::abs(p.Mtau) < eps) {
        setCorrectionFlag(EFTOrders::G12YTAU2, 1);
        setCorrectionFlag(EFTOrders::G22YTAU2, 1);
        setCorrectionFlag(EFTOrders::YTAU4   , 0);
        setCorrectionFlag(EFTOrders::YTAU2YB4, 0);
        setCorrectionFlag(EFTOrders::YTAU2YT4, 0);
        setCorrectionFlag(EFTOrders::YTAU6   , 0);
        setCorrectionFlag(EFTOrders::YTAU4YB2, 0);
    }

    // For now, disable all 1L corrections, except 1L O(at)
    setCorrectionFlag(EFTOrders::G14     , 0);
    setCorrectionFlag(EFTOrders::G24     , 0);
    setCorrectionFlag(EFTOrders::G12G22  , 0);
    setCorrectionFlag(EFTOrders::G12YT2  , 0);
    setCorrectionFlag(EFTOrders::G22YT2  , 0);
    setCorrectionFlag(EFTOrders::G12YB2  , 0);
    setCorrectionFlag(EFTOrders::G22YB2  , 0);
    setCorrectionFlag(EFTOrders::G12YTAU2, 0);
    setCorrectionFlag(EFTOrders::G22YTAU2, 0);
    setCorrectionFlag(EFTOrders::YB4, 0);
    setCorrectionFlag(EFTOrders::YTAU4, 0);

    // For now, disable all 2L corrections, except 2L O(at*as + at^2)
    setCorrectionFlag(EFTOrders::G32YB4  , 0);
    setCorrectionFlag(EFTOrders::YB6     , 0);
    setCorrectionFlag(EFTOrders::YTAU2YB4, 0);
    setCorrectionFlag(EFTOrders::YT2YB4  , 0);
    setCorrectionFlag(EFTOrders::YB2YT4  , 0);
    setCorrectionFlag(EFTOrders::YTAU4YB2, 0);
    setCorrectionFlag(EFTOrders::YTAU6   , 0);
}

void Mh2EFTCalculator::setCorrectionFlag(int order, int flag)
{
    if (flag < 0 || flag > 1)
        ERROR_MSG("You can only enable (1) or disable (0) corrections!");

    if (order < EFTOrders::EFTOrders::FIRST ||
            order >= EFTOrders::EFTOrders::NUMBER_OF_EFT_ORDERS)
        ERROR_MSG("Your order is not defined in the EFTOrders enum!");

    orders.at(order) = flag;
}

/**
 * Returns the tree-level EFT contribution to the light CP-even Higgs mass
 */
double Mh2EFTCalculator::getDeltaMh2EFT0Loop() const
{
    return pow2(p.MZ * std::cos(2 * std::atan(p.vu / p.vd)));
}

/**
 * Returns the 1-loop EFT contribution to the light CP-even Higgs mass
 *
 * @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 * @return 1-loop EFT contribution to the light CP-even Higgs mass
 */
double Mh2EFTCalculator::getDeltaMh2EFT1Loop(int omitSMLogs, int omitMSSMLogs) const
{
    ThresholdCalculator thresholdCalculator(p, msq2);

    using std::log;
    const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));

    const double v2 = pow2(p.vu) + pow2(p.vd);
    const double gt = sqrt(2) * p.Mt / std::sqrt(v2);

    // 1-Loop prefactor at
    const double pref_at = 1. / pow2(4 * Pi) * pow2(p.Mt * gt);

    const double q2 = pow2(p.Mt);
    const double beta = atan(p.vu / p.vd);
    const double cbeta = cos(beta);
    const double c2beta = cos(2 * beta);
    const double sbeta = sin(beta);
    const double mhtree = std::abs(c2beta * p.MZ);
    const double yt = sqrt(2.) * p.Mt / p.vu;
    const double yb = sqrt(2.) * p.Mb / p.vd;
    const double ytau = sqrt(2.) * p.Mtau / p.vd;
    const int Xi = 1;        // gauge parameter

    // Threshold corrections
    const double dlambdayb2g12 = CALC_IF(orders.at(EFTOrders::G12YB2),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_YB2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdag14 = CALC_IF(orders.at(EFTOrders::G14),
                                      thresholdCalculator.getThresholdCorrection(
                                          ThresholdVariables::LAMBDA_G14, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdaregg14 = CALC_IF(orders.at(EFTOrders::G14),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_REG_G14, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdachig14 = CALC_IF(orders.at(EFTOrders::G14),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_CHI_G14, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dg1g1 = CALC_IF(orders.at(EFTOrders::G14),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::G1_G1, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdachig24 = CALC_IF(orders.at(EFTOrders::G24),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_CHI_G24, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdag24 = CALC_IF(orders.at(EFTOrders::G24),
                                      thresholdCalculator.getThresholdCorrection(
                                          ThresholdVariables::LAMBDA_G24, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dg2g2 = CALC_IF(orders.at(EFTOrders::G24),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::G2_G2, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdaregg24 = CALC_IF(orders.at(EFTOrders::G24),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_REG_G24, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdag12g22 = CALC_IF(orders.at(EFTOrders::G12G22),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_G12_G22, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdaregg12g22 = CALC_IF(orders.at(EFTOrders::G12G22),
                                            thresholdCalculator.getThresholdCorrection(
                                                    ThresholdVariables::LAMBDA_REG_G12_G22, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdachig12g22 = CALC_IF(orders.at(EFTOrders::G12G22),
                                            thresholdCalculator.getThresholdCorrection(
                                                    ThresholdVariables::LAMBDA_CHI_G12_G22, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayb2g22 = CALC_IF(orders.at(EFTOrders::G22YB2),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_YB2_G22, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayb4 = CALC_IF(orders.at(EFTOrders::YB4),
                                      thresholdCalculator.getThresholdCorrection(
                                          ThresholdVariables::LAMBDA_YB4, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayt2g12 = CALC_IF(orders.at(EFTOrders::G12YT2),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_YT2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayt2g22 = CALC_IF(orders.at(EFTOrders::G22YT2),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_YT2_G22, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdaytau2g12 = CALC_IF(orders.at(EFTOrders::G12YTAU2),
                                           thresholdCalculator.getThresholdCorrection(
                                                   ThresholdVariables::LAMBDA_YTAU2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdaytau2g22 = CALC_IF(orders.at(EFTOrders::G22YTAU2),
                                           thresholdCalculator.getThresholdCorrection(
                                                   ThresholdVariables::LAMBDA_YTAU2_G22, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdaytau4 = CALC_IF(orders.at(EFTOrders::YTAU4),
                                        thresholdCalculator.getThresholdCorrection(
                                            ThresholdVariables::LAMBDA_YTAU4, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dvg12 = 0.;/*
        CALC_IF(orders.at(EFTOrders::G12G22) || orders.at(EFTOrders::G14),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::VEV_G12, RenSchemes::DRBARPRIME, omitMSSMLogs));*/
    const double dvg22 = 0.;/*
        CALC_IF(orders.at(EFTOrders::G12G22) || orders.at(EFTOrders::G24),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::VEV_G22, RenSchemes::DRBARPRIME, omitMSSMLogs));*/
    const double dvyt2 =
        CALC_IF(orders.at(EFTOrders::G12YT2) || orders.at(EFTOrders::G22YT2),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::VEV_YT2, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dvyb2 =
        CALC_IF(orders.at(EFTOrders::G12YB2) || orders.at(EFTOrders::G22YB2),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::VEV_YB2, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dvytau2 =
        CALC_IF(orders.at(EFTOrders::G12YTAU2) || orders.at(EFTOrders::G22YTAU2),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::VEV_YTAU2, RenSchemes::DRBARPRIME, omitMSSMLogs));

    const double bbhDR = fB(pow2(mhtree), pow2(mhtree), q2);
    const double bbwDR = fB(pow2(mhtree), pow2(p.MW), q2);
    const double bbzDR = fB(pow2(mhtree), pow2(p.MZ), q2);
    const double B00DR = fB(pow2(mhtree), 0., q2);

    // corrections to Mh2
    const double dmh2g12g22 = isNaN(orders.at(EFTOrders::G12G22) * (
                                        (v2 * (24 + 40 * dlambdachig12g22 + 40 * dlambdag12g22 + 40 *
                                                dlambdaregg12g22 - 36 * lmMt + 40 * dvg12 * pow2(c2beta) + 24 * dvg22 * pow2(
                                                        c2beta) + 36 * lmMt * pow2(c2beta) - 12 * lmMt * Xi * pow2(c2beta) - 6 * bbwDR * (-2
                                                                + pow2(c2beta)) * pow2(c2beta) - 6 * pow4(c2beta) - 27 * bbhDR * pow4(c2beta) -
                                                36 * lmMt * pow4(c2beta) - 3 * bbzDR * (12 - 4 * pow2(c2beta) + pow4(c2beta)))) /
                                        80.
                                    ),
                                    "dmh2g12g22");

    const double dmh2g14 = isNaN(orders.at(EFTOrders::G14) * (
                                     -(v2 * (-4 * (18 + 100 * dlambdachig14 + 100 * dlambdag14 + 100 *
                                             dlambdaregg14 - 27 * lmMt) + 6 * (40 * dg1g1 - 40 * dvg12 + 3 * lmMt * (-3 + Xi)) *
                                             pow2(c2beta) + 9 * (9 * bbhDR + 2 * (1 + bbwDR + 6 * lmMt)) * pow4(c2beta) + 9 *
                                             bbzDR * (12 - 4 * pow2(c2beta) + pow4(c2beta)))) / 800.
                                 ),
                                 "dmh2g14");
    const double dmh2g24 = isNaN(orders.at(EFTOrders::G24) * (
                                     -(v2 * (-24 - 16 * dlambdachig24 - 16 * dlambdag24 - 16 * dlambdaregg24 +
                                             36 * lmMt + 16 * dg2g2 * pow2(c2beta) - 16 * dvg22 * pow2(c2beta) - 18 * lmMt * pow2(
                                                     c2beta) + 6 * lmMt * Xi * pow2(c2beta) + 2 * pow4(c2beta) + 9 * bbhDR * pow4(
                                                             c2beta) + 12 * lmMt * pow4(c2beta) + 2 * bbwDR * (12 - 4 * pow2(c2beta) + pow4(
                                                                     c2beta)) + bbzDR * (12 - 4 * pow2(c2beta) + pow4(c2beta)))) / 32.
                                 ),
                                 "dmh2g24");
    const double dmh2g12yb2 = isNaN(orders.at(EFTOrders::G12YB2) * (((
                                        10 * dlambdayb2g12 - 3 * (3 * B00DR - 2 * dvyb2 + 3 * lmMt) * pow2(c2beta))
                                    * pow2(cbeta) * v2) / 20.),
                                    "dmh2g12yb2");
    const double dmh2g22yb2 = isNaN(orders.at(EFTOrders::G22YB2) * (((
                                        2 * dlambdayb2g22 + (-3 * B00DR + 2 * dvyb2 - 3 * lmMt) * pow2(c2beta)) * pow2(
                                        cbeta) * v2) / 4.),
                                    "dmh2g22yb2");
    const double dmh2yb4 = isNaN(orders.at(EFTOrders::YB4) * (((12 * B00DR
                                 + dlambdayb4 + 12 * lmMt) * v2 * pow4(cbeta)) / 2.),
                                 "dmh2yb4");
    const double dmh2g12yt2 = isNaN(orders.at(EFTOrders::G12YT2) * (((10
                                    * dlambdayt2g12 + (6 + 6 * dvyt2 - 9 * lmMt) * pow2(c2beta)) * pow2(sbeta) * v2)
                                    / 20.),
                                    "dmh2g12yt2");
    const double dmh2g22yt2 = isNaN(orders.at(EFTOrders::G22YT2) * (((2
                                    * dlambdayt2g22 + (2 + 2 * dvyt2 - 3 * lmMt) * pow2(c2beta)) * pow2(sbeta)
                                    * v2) / 4.),
                                    "dmh2g22yt2");
    const double dmh2yt4 = isNaN(orders.at(EFTOrders::YT4) * (pref_at * (12 * lmMt +
                                 thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT,
                                         RenSchemes::DRBARPRIME, omitMSSMLogs))),
                                 "dmh2yt4");
    const double dmh2g12ytau2 = isNaN(orders.at(EFTOrders::G12YTAU2) * (-((-10
                                      * dlambdaytau2g12 + 3 * B00DR * pow2(c2beta) + 3 * (-2 * dvytau2 + lmMt) *
                                      pow2(c2beta)) * pow2(cbeta) * v2) / 20.),
                                      "dmh2g12ytau2");
    const double dmh2g22ytau2 = isNaN(orders.at(EFTOrders::G22YTAU2) * (-((-2
                                      * dlambdaytau2g22 + B00DR * pow2(c2beta) + (-2 * dvytau2 + lmMt) * pow2(
                                          c2beta)) * pow2(cbeta) * v2) / 4.),
                                      "dmh2g22ytau2");
    const double dmh2ytau4 = isNaN(orders.at(EFTOrders::YTAU4) * (((4 * B00DR
                                   + dlambdaytau4 + 4 * lmMt) * v2 * pow4(cbeta)) / 2.),
                                   "dmh2ytau4");

    // Loop factor
    const double k = 1 / pow2(4.*Pi);

    return dmh2yt4 + k * (pow2(p.g1 * p.g2) * dmh2g12g22 + pow4(p.g1) * dmh2g14 + pow4(p.g2) *
                          dmh2g24 + pow2(p.g1 * yb) * dmh2g12yb2 + pow2(p.g2 * yb) * dmh2g22yb2 + pow4(yb) *
                          dmh2yb4 + pow2(p.g1 * ytau) * dmh2g12ytau2 + pow2(p.g2 * ytau) * dmh2g22ytau2 +
                          pow4(ytau) * dmh2ytau4 + pow2(p.g1 * yt) * dmh2g12yt2 + pow2(p.g2 * yt) * dmh2g22yt2);
}

/**
 * Returns the 2-loop EFT contribution to the light CP-even Higgs mass
 * @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 * @return 2-loop EFT contribution to the light CP-even Higgs mass
 */
double Mh2EFTCalculator::getDeltaMh2EFT2Loop(int omitSMLogs, int omitMSSMLogs) const
{
    ThresholdCalculator thresholdCalculator(p, msq2);

    using std::log;
    const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
    // couplings
    const double v2 = pow2(p.vu) + pow2(p.vd);
    const double gt = sqrt(2) * p.Mt / std::sqrt(v2);
    const double g32 = pow2(p.g3);
    const double yt2 = pow2(sqrt(2.) * p.Mt / p.vu);
    const double yb2 = pow2(sqrt(2.) * p.Mb / p.vd);
    const double ytau2 = pow2(sqrt(2.) * p.Mtau / p.vd);
    const double yt4 = pow2(yt2);
    const double yb4 = pow2(yb2);
    const double yt6 = pow3(yt2);
    const double yb6 = pow3(yb2);
    const double ytau4 = pow2(ytau2);
    const double ytau6 = pow3(ytau2);
    const double beta = atan(p.vu / p.vd);
    const double cbeta = cos(beta);
    const double sbeta = sin(beta);
    const double lmbMt = log(pow2(p.Mb / p.Mt));

    // 2-Loop prefactor at*as
    const double pref = 1. / pow4(4 * Pi) * pow2(p.Mt * gt * p.g3);
    const double B00DR = 0.;

    // Threshold corrections
    const double dytas = CALC_IF(orders.at(EFTOrders::G32YT4),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::YT_AS, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayb4g32 = CALC_IF(orders.at(EFTOrders::G32YB4),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_YB4_G32, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayb4 =
        CALC_IF(orders.at(EFTOrders::YB6) || orders.at(EFTOrders::YTAU2YB4) ||
                orders.at(EFTOrders::YT2YB4) || orders.at(EFTOrders::G32YB4),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::LAMBDA_YB4, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayb6 = CALC_IF(orders.at(EFTOrders::YB6),
                                      thresholdCalculator.getThresholdCorrection(
                                          ThresholdVariables::LAMBDA_YB6, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayt4 =
        CALC_IF(orders.at(EFTOrders::YT6) || orders.at(EFTOrders::YB2YT4),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::LAMBDA_AT, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dytyt = CALC_IF(orders.at(EFTOrders::YT6),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::YT_YT, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayt6 = CALC_IF(orders.at(EFTOrders::YT6),
                                      thresholdCalculator.getThresholdCorrection(
                                          ThresholdVariables::LAMBDA_YT6, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dvyt2 = CALC_IF(orders.at(EFTOrders::YT6) || orders.at(EFTOrders::YT2YB4),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::VEV_YT2, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dytauytau = CALC_IF(orders.at(EFTOrders::YTAU6),
                                     thresholdCalculator.getThresholdCorrection(
                                         ThresholdVariables::YTAU_YTAU, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdaytau4 =
        CALC_IF(orders.at(EFTOrders::YTAU4YB2) || orders.at(EFTOrders::YTAU6),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::LAMBDA_YTAU4, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdaytau6 = CALC_IF(orders.at(EFTOrders::YTAU6),
                                        thresholdCalculator.getThresholdCorrection(
                                            ThresholdVariables::LAMBDA_YTAU6, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayt2yb4 = CALC_IF(orders.at(EFTOrders::YT2YB4),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_YT2_YB4, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayt4yb2 = CALC_IF(orders.at(EFTOrders::YB2YT4),
                                         thresholdCalculator.getThresholdCorrection(
                                                 ThresholdVariables::LAMBDA_YT4_YB2, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dytyb = CALC_IF(orders.at(EFTOrders::YB2YT4),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::YT_YB, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dvytau2 =
        CALC_IF(orders.at(EFTOrders::YTAU2YB4) || orders.at(EFTOrders::YTAU6),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::VEV_YTAU2, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dvyb2 =
        CALC_IF(orders.at(EFTOrders::YB6) || orders.at(EFTOrders::YTAU4YB2) ||
                orders.at(EFTOrders::YB2YT4),
                thresholdCalculator.getThresholdCorrection(
                    ThresholdVariables::VEV_YB2, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dytauyb = CALC_IF(orders.at(EFTOrders::YTAU4YB2),
                                   thresholdCalculator.getThresholdCorrection(
                                       ThresholdVariables::YTAU_YB, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayb4ytau2 = CALC_IF(orders.at(EFTOrders::YTAU2YB4),
                                           thresholdCalculator.getThresholdCorrection(
                                                   ThresholdVariables::LAMBDA_YTAU2_YB4, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dlambdayb2ytau4 = CALC_IF(orders.at(EFTOrders::YTAU4YB2),
                                           thresholdCalculator.getThresholdCorrection(
                                                   ThresholdVariables::LAMBDA_YTAU4_YB2, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dybyt = CALC_IF(orders.at(EFTOrders::YT2YB4),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::YB_YT, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dybas = CALC_IF(orders.at(EFTOrders::G32YB4),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::YB_AS, RenSchemes::DRBARPRIME, omitMSSMLogs));
    const double dybyb = CALC_IF(orders.at(EFTOrders::YB6),
                                 thresholdCalculator.getThresholdCorrection(
                                     ThresholdVariables::YB_YB, RenSchemes::DRBARPRIME, omitMSSMLogs));

    // Corrections to Mh
    const double dmh2yt4g32 = isNaN(orders.at(EFTOrders::G32YT4) * (pref * (96 * pow2(lmMt)
                                    + (-32 + 48 * dytas) * lmMt - 24 * dytas
                                    + thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT_AS,
                                            RenSchemes::DRBARPRIME, omitMSSMLogs))),
                                    "dmh2yt4g32");
    const double dmh2yb4g32 =
        isNaN(orders.at(EFTOrders::G32YB4) * ((dlambdayb4g32
                + 16 * (3 * B00DR * (dybas + 16 * lmMt) + lmMt * (-2 + 3 * dybas
                        + 24 * lmMt))) * v2 * pow4(cbeta)) / 2.,
              "dmh2yb4g32");
    const double dmh2yb6 =
        isNaN(orders.at(EFTOrders::YB6) * (-((-144 * dybyb * (B00DR
                                           + lmMt) + pow2(cbeta) * (-49 - 3 * dlambdayb6 + 72 * dvyb2
                                                   - 72 * B00DR * dvyb2 - 60 * lmbMt + 234 * lmMt + 1296 * B00DR * lmMt - 72 * dvyb2
                                                   * lmMt + dlambdayb4 * (-9 + 9 * B00DR - 6 * dvyb2 + 9 * lmMt) + 36 * pow2(lmbMt)
                                                   + 648 * pow2(lmMt) - 6 * pow2(Pi))) * v2 * pow4(cbeta)) / 6.),
              "dmh2yb6");
    const double dmh2yt6 = isNaN(orders.at(EFTOrders::YT6) * ((v2 * (4 * dytyt * (
                                     -6 + dlambdayt4 + 12 * lmMt) + (-12 + dlambdayt6 + dlambdayt4 * (2 + 2 * dvyt2
                                             - 3 * lmMt) + 24 * dvyt2 * (-1 + lmMt) - 18 * lmMt * (1 + 3 * lmMt) - 2 * pow2(Pi))
                                 * pow2(sbeta)) * pow4(sbeta)) / 2.),
                                 "dmh2yt6");
    const double dmh2yb4ytau2 = isNaN(orders.at(EFTOrders::YTAU2YB4) * (((
                                          dlambdayb4ytau2 - 96 * B00DR * lmMt - dlambdayb4 * (-1 + B00DR + lmMt) + 2
                                          * dvytau2 * (-12 + 12 * B00DR + dlambdayb4 + 12 * lmMt) - 48 * pow2(lmMt)) * v2
                                      * pow6(cbeta)) / 2.),
                                      "dmh2yb4ytau2");
    const double dmh2yb2ytau4 = isNaN(orders.at(EFTOrders::YTAU4YB2) * (((4 * dytauyb
                                      * (4 * B00DR + dlambdaytau4 + 4 * lmMt) + pow2(cbeta) * (dlambdayb2ytau4 - 8 * dvyb2
                                              + 8 * B00DR * dvyb2 + dlambdaytau4 * (3 - 3 * B00DR + 2 * dvyb2 - 3 * lmMt) - 6 * lmMt
                                              - 24 * B00DR * lmMt + 8 * dvyb2 * lmMt - 12 * pow2(lmMt))) * v2 * pow4(cbeta)) / 2.),
                                      "dmh2yb2ytau4");
    const double dmh2ytau6 = isNaN(orders.at(EFTOrders::YTAU6) * (-((-12 * dytauytau
                                   * (4 * B00DR + dlambdaytau4 + 4 * lmMt) + pow2(cbeta) * (-3 * dlambdaytau6 + 3
                                           * dlambdaytau4 * (-1 + B00DR - 2 * dvytau2 + lmMt) + 2 * (6 + 30 * (1 + B00DR) * lmMt
                                                   - 12 * dvytau2 * (-1 + B00DR + lmMt) + 15 * pow2(lmMt) + pow2(Pi)))) * v2
                                   * pow4(cbeta)) / 6.),
                                   "dmh2ytau6");
    const double dmh2yt2yb4 =
        isNaN(orders.at(EFTOrders::YT2YB4) * (((48 * dybyt * lmMt
                + (dlambdayt2yb4 + dlambdayb4 * (2 + 2 * dvyt2 - 3 * lmMt) + 3 * (-15 + 8 * lmbMt
                        + 8 * dvyt2 * (-1 + lmMt) + 18 * lmMt - 24 * pow2(lmMt) - 2 * pow2(Pi))) * pow2(sbeta)
                + 12 * B00DR * (4 * dybyt + (2 * dvyt2 - 9 * lmMt) * pow2(sbeta))) * v2 * pow4(cbeta)) / 2.),
              "dmh2yt2yb4");
    const double dmh2yt4yb2 = isNaN(orders.at(EFTOrders::YB2YT4) * (((4 * dytyb * (-6
                                    + dlambdayt4 + 12 * lmMt) + pow2(cbeta) * (dlambdayt4yb2 + dlambdayt4 * (3 - 3 * B00DR
                                            + 2 * dvyb2 - 3 * lmMt) - 6 * (4 * dvyb2 + lmMt + 6 * B00DR * lmMt - 4 * dvyb2 * lmMt
                                                    + 3 * pow2(lmMt) - pow2(Pi)))) * v2 * pow4(sbeta)) / 2.),
                                    "dmh2yt4yb2");

    // Loop factor
    const double k2 = 1 / pow4(4.*Pi);

    return k2 * (g32 * yb4 * dmh2yb4g32 + yb6 * dmh2yb6 + yt6 * dmh2yt6 + yb4 * ytau2
                 * dmh2yb4ytau2 + yb2 * ytau4 * dmh2yb2ytau4 + ytau6 * dmh2ytau6 + yt2 * yb4
                 * dmh2yt2yb4 + yt4 * yb2 * dmh2yt4yb2) + dmh2yt4g32;
}

/**
 * Returns the 3-loop EFT contribution to the light CP-even Higgs mass.
 *
 * @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 * @param omitDeltaLambda3L an integer flag to disable the MSSM contribution to delta_lambda_3L
 * @return 3-loop contribution to the light CP-even Higgs mass
 */
double Mh2EFTCalculator::getDeltaMh2EFT3Loop(
    int omitSMLogs, int omitMSSMLogs, int omitDeltaLambda3L) const
{
    ThresholdCalculator thresholdCalculator(p, msq2);

    using std::log;
    using std::sqrt;

    const double catas2 = 248.1215180432007;
    const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
    // threshold correction of yt_as DRbar'
    const double dytas = thresholdCalculator.getThresholdCorrection(
                             ThresholdVariables::YT_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);
    // threshold correction of yt_as2 DRbar'
    const double dytas2 = thresholdCalculator.getThresholdCorrection(
                              ThresholdVariables::YT_AS2, RenSchemes::DRBARPRIME, omitMSSMLogs);
    // threshold correction of g3_as DRbar'
    const double dg3as = thresholdCalculator.getThresholdCorrection(
                             ThresholdVariables::G3_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);

    const double gt = sqrt(2.) * p.Mt / sqrt(pow2(p.vu) + pow2(p.vd));

    // 3-Loop prefactor at*as^2
    const double pref = 1. / pow6(4 * Pi) * pow2(p.Mt * gt * pow2(p.g3));

    return pref * (736 * pow3(lmMt) + (160 + 192 * dg3as + 384 * dytas) * pow2(lmMt)
                   + (-128 * zt3 - 2056 / 3. + -64 * dg3as - 512 * dytas + 72 * pow2(dytas)
                      + 48 * dytas2) * lmMt + 64 * dytas - 84 * pow2(dytas) - 24 * dytas2
                   + catas2
                   + omitDeltaLambda3L * thresholdCalculator.getThresholdCorrection(
                       ThresholdVariables::LAMBDA_AT_AS2,
                       RenSchemes::DRBARPRIME, omitMSSMLogs));
}

/**
 * Returns the matching relation of delta_Lambda 3L for the degenerate
 * mass case.
 *
 * @param scale the renormalization scale
 * @param mst1 the mass of the light stop quark
 * @param Xt stop mixing parameter
 * @param omitlogs factor which multiplies the logs
 * @return delta_Lambda 3L
 */
double Mh2EFTCalculator::getDeltaLambdaDegenerate(
    double scale, double mst1, double Xt, int omitlogs) const
{
    using std::log;

    const double LS = omitlogs * log(pow2(scale / p.MSt(0))) + log(pow2(p.MSt(0) / mst1));

    const double gt = sqrt(2) * p.Mt / std::sqrt(pow2(p.vu) + pow2(p.vd));

    // 3-Loop prefactor
    const double pref = 1. / pow6(4 * Pi) * pow2(p.Mt * gt * pow2(p.g3));

    // to obtain delta_lambda one has to divide the difference of the two calculations by v^2
    const double v2 = pow2(p.vd) + pow2(p.vu);

    const double xt = Xt / mst1;
    const double catas2 = 248.1215180432007;

    const double deltaLambda3L = 2 / 27.*pref * (6082 - 27832 * LS + 14856 * pow2(LS)
                                 - 4032 * pow3(LS) - 15408 * zt3 + 1728 * zt3 * LS - 27 * catas2 / 2.
                                 + xt * (7616 * LS - 11712 * pow2(LS) + 32 * (-940 + 477 * zt3))
                                 + pow2(xt) * (28848 - 2640 * LS + 1008 * pow2(LS) - 11880 * zt3)
                                 + pow3(xt) * (160 * LS + 864 * pow2(LS) + 8 * (2722 - 2259 * zt3))) / v2;

    return deltaLambda3L;
}

/**
 * Calculates the loop corrections in the approximation v^2 << MS^2
 * and prints the result.
 *
 * @param ostr output stream
 * @param mhc Mh2EFTCalculator object
 *
 * @return output stream
 */
std::ostream& operator<<(std::ostream& ostr, const Mh2EFTCalculator& mhc)
{
    ostr << "=========================\n"
         << "Himalaya Mh2EFTCalculator\n"
         << "=========================\n"
         << mhc.p
         << "msq2  = " << mhc.msq2 << " GeV^2\n";

    const auto dmh2_0l =  mhc.getDeltaMh2EFT0Loop();
    const auto dmh2_1l =  mhc.getDeltaMh2EFT1Loop(1, 1);
    const auto dmh2_2l =  mhc.getDeltaMh2EFT2Loop(1, 1);

    ostr << "Mh^2_EFT_0L  = " << dmh2_0l << " GeV^2 O(g1^2, g2^2)\n";
    ostr << "ΔMh^2_EFT_1L = " << dmh2_1l << " GeV^2 O(full)\n";
    ostr << "ΔMh^2_EFT_2L = " << dmh2_2l << " GeV^2 O((αt+ab)*αs + (αt+αb)^2 + ab*aτ + aτ^2)\n";

    return ostr;
}

} // namespace mh2_eft
} // namespace himalaya
