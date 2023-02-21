#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "hermiteRootsWeights.h"
#include "sphericalBessel.h"

double T(double t) {
    double mult = pow(t, 2.0 / 3.0);
    double term0 = 1.0;
    double term2 = (5.0 / 48.0) * pow(t, -2);
    double term4 = -(5.0 / 36.0) * pow(t, -4);
    double term6 = (77125.0 / 82944.0) * pow(t, -6);
    double term8 = -(108056875.0 / 6967296.0) * pow(t, -8);
    double term10 = (162375596875.0 / 334430208.0) * pow(t, -10);
    double result = mult * (term0 + term2 + term4 + term6 + term8 + term10);
    return result;
}

void GetAiryRoots(std::vector<double> &airyRoots, const std::vector<double> &first11AiryRoots, std::size_t m) {
    airyRoots.resize(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        double x = static_cast<double>(i + 1);
        double t = (3.0 * M_PI / 8.0) * (4.0 * x - 1.0);
        airyRoots[i] = -T(t);
    }
    for (std::size_t i = 0; i < 11; ++i) airyRoots[i] = first11AiryRoots[i];
}

void GetXinit_Gatteschi(std::vector<double> &xinitAiry, const std::vector<double> &airyRoots, double a, double nu) {
    std::size_t size = airyRoots.size();
    std::vector<double> xinit(size, 0.0);

    for (std::size_t k = 0; k < size; ++k) {
        double aR = airyRoots[k];
        double aR2 = aR * aR;
        double aR3 = aR2 * aR;
        double aR4 = aR3 * aR;
        double aR5 = aR4 * aR;
        double xSq1 = nu + pow(2.0, 2.0 / 3.0) * aR * pow(nu, 1.0 / 3.0);
        double xSq2 = (1.0 / 5.0) * pow(2.0, 4.0 / 3.0) * aR2 * pow(nu, -1.0 / 3.0);
        double xSq3 = (11.0 / 35.0 - a * a - 12.0 / 175.0) * aR3 * pow(nu, -1.0);
        double xSq4 = ((16.0 / 1575.0) * aR + (92.0 / 7875.0) * aR4) * pow(2.0, 2.0 / 3.0) * pow(nu, -5.0 / 3.0);
        double xSq5 =
            ((15152.0 / 3031875.0) * aR5 + (1088.0 / 121275.0) * aR2) * pow(2.0, 1.0 / 3.0) * pow(nu, -7.0 / 3.0);
        double sum = xSq1 + xSq2 + xSq3 + xSq4 - xSq5;
        xinit[k] = sqrt(abs(sum));
    }
    std::reverse(xinit.begin(), xinit.end());
    xinitAiry = xinit;
}

void GetXinit_Tricomi(std::vector<double> &xinitSin, std::size_t m, double a, double nu) {
    std::vector<double> Tnk0(m, M_PI / 2);

    std::vector<double> rhs(m, 0.0);
    for (std::size_t j = 0; j < m; ++j) {
        double term1 = 4.0 * m + 3.0;
        double term2 = 4.0 * (j + 1);
        rhs[j] = M_PI * (term1 - term2) / nu;
    }
    // first get roots to x - sinx = rhs
    for (std::size_t k = 1; k <= 7; ++k) {
        for (std::size_t j = 0; j < m; ++j) {
            double tnk0 = Tnk0[j];
            double val = tnk0 - sin(tnk0) - rhs[j];
            double dval = 1.0 - cos(tnk0);
            double dtnk0 = val / dval;
            Tnk0[j] = tnk0 - dtnk0;
        }
    }
    std::vector<double> tnk(m, 0.0);
    for (std::size_t j = 0; j < m; ++j) {
        double t = cos(Tnk0[j] / 2.0);
        tnk[j] = t * t;
    }

    xinitSin.resize(m, 0.0);
    for (std::size_t j = 0; j < m; ++j) {
        double term1 = nu * tnk[j];
        double term2a = (tnk[j] + 1.0 / 4.0) / ((tnk[j] - 1.0) * (tnk[j] - 1.0));
        double term2b = 3.0 * a * a - 1.0;
        double term2 = (term2a + term2b) / (3.0 * nu);
        double term = term1 - term2;
        xinitSin[j] = sqrt(term);
    }
}

bool GetInitialGuessesForHermiteRoots(std::vector<double> &initialGuesses, std::size_t n) {
    // [1] L. Gatteschi, Asymptotics and bounds for the zeros of Laguerre
    // polynomials: a survey, J. Comput. Appl. Math., 144 (2002), pp. 7-27.
    //
    // [2] F. G. Tricomi, Sugli zeri delle funzioni di cui si conosce una
    // rappresentazione asintotica, Ann. Mat. Pura Appl. 26 (1947), pp. 283-300.

    // Error if n < 20 because initial guesses are based on asymptotic expansions:
    // assert n ≥ 20

    initialGuesses.clear();
    if (n <= 20) {
        std::cout << "n too small for guesses based on asymptotics" << std::endl;
        return false;
    }

    std::vector<double> first11AiryRoots(11, 0.0);
    first11AiryRoots[0] = -2.338107410459767;
    first11AiryRoots[1] = -4.08794944413097;
    first11AiryRoots[2] = -5.520559828095551;
    first11AiryRoots[3] = -6.786708090071759;
    first11AiryRoots[4] = -7.944133587120853;
    first11AiryRoots[5] = -9.022650853340981;
    first11AiryRoots[6] = -10.04017434155809;
    first11AiryRoots[7] = -11.00852430373326;
    first11AiryRoots[8] = -11.93601556323626;
    first11AiryRoots[9] = -12.828776752865757;
    first11AiryRoots[10] = -13.69148903521072;

    std::vector<double> initialGuessGatteschi;
    // Gatteschi formula involving airy roots [1].
    // These initial guess are good near x = sqrt(n+1/2)
    std::size_t m = 0;
    double a = 0.0;
    if (n % 2 == 1) {
        m = (n - 1) / 2;
        // bess = (1:m)*π
        a = 0.5;
    } else {
        m = n / 2;
        // bess = ((0:m-1) .+ 0.5)*π
        a = -0.5;
    }
    double nu = 4.0 * m + 2.0 * a + 2.0;
    std::vector<double> airyRoots;
    GetAiryRoots(airyRoots, first11AiryRoots, m);
    std::vector<double> xInitAiry;
    std::vector<double> xinitSin;
    GetXinit_Gatteschi(xInitAiry, airyRoots, a, nu);
    GetXinit_Tricomi(xinitSin, m, a, nu);

    // patch together
    double p = 0.4985 + std::numeric_limits<double>::epsilon();
    std::size_t upperSin = static_cast<std::size_t>(floor(p * n));
    std::cout << "upperSin = " << upperSin << std::endl;

    initialGuesses.resize(m, 0.0);
    for (std::size_t j = 0; j < upperSin; ++j) initialGuesses[j] = xinitSin[j];
    for (std::size_t j = upperSin; j < m; ++j) initialGuesses[j] = xInitAiry[j];

    return true;
}

void GetHermitePolynomialValByRecursion(double &herm, double &derivHerm, std::size_t n, double x) {
    // from the julia package by Alex Townsend Cornell Univ.
    double dn = static_cast<double>(n);
    double w = exp(-x * x / (4.0 * dn));
    std::size_t wc = 0;  // num times a correction is applied for scaling
    double Hold = 1.0;
    if (n == 0) {
        herm = 1.0;
        derivHerm = 0.0;
        return;
    }
    double H = x;
    for (std::size_t k = 1; k < n; ++k) {
        double dk = static_cast<double>(k);
        double tmp1 = x * H / sqrt(dk + 1.0);
        double tmp2 = Hold / sqrt(1.0 + 1.0 / dk);
        Hold = H;
        H = tmp1 - tmp2;
        while ((abs(H) >= 100) && (wc < n)) {
            H *= w;
            Hold *= w;
            wc += 1;
        }
    }
    for (std::size_t j = wc + 1; j <= n; ++j) {
        H *= w;
        Hold *= w;
    }
    herm = H;
    derivHerm = -x * H + sqrt(dn) * Hold;
}

void GetRoot_and_WeightOfHermite(double &root, double &weight, std::size_t n, double guessRoot) {
    // newton's method
    double herm = 0.0;
    double derivHerm = 0.0;
    double x = sqrt(2.0) * guessRoot;
    for (std::size_t j = 0; j < 10; ++j) {
        GetHermitePolynomialValByRecursion(herm, derivHerm, n, x);
        double dx = herm / derivHerm;
        if (isnan(dx)) dx = 0.0;
        x = x - dx;
        if (abs(dx) < std::numeric_limits<double>::epsilon()) break;
    }
    x /= sqrt(2.0);
    root = x;
    weight = 1.0 / (derivHerm * derivHerm);
}

void GetRoots_and_WeightsOfHermite(std::vector<double> &roots, std::vector<double> &weights, std::size_t n,
                                   const std::vector<double> &initialGuesses) {
    roots.clear();
    weights.clear();
    std::size_t size = initialGuesses.size();
    roots.reserve(size);
    weights.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        double guessRoot = initialGuesses[i];
        double root = 0.0;
        double weight = 0.0;
        GetRoot_and_WeightOfHermite(root, weight, n, guessRoot);
        roots.push_back(root);
        weights.push_back(weight);
    }
    std::reverse(roots.begin(), roots.end());
    std::reverse(weights.begin(), weights.end());
    double sum = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) {
        sum += 2.0 * exp(-roots[i] * roots[i]) * weights[i];  // factor of 2 because we don't count neg roots
    }

    for (std::size_t i = 0; i < weights.size(); ++i) weights[i] *= sqrt(M_PI) * exp(-roots[i] * roots[i]) / sum;
}

double GetEwaldCoeff(double THRESH, double cutoff) {
    // find beta so that erfc(beta*cutoff)/cutoff = 10^(-THRESH)
    double tol = pow(10.0, -THRESH);
    std::cout << "tol = " << tol << std::endl;

    double betaLo = 0.0;
    double betaHi = 1.0;
    double fBetaLo = erfc(betaLo * cutoff) / cutoff;
    double fBetaHi = erfc(betaHi * cutoff) / cutoff;
    if (fBetaLo <= tol) return betaLo;
    if (fBetaHi > tol) {
        while (fBetaHi > tol) {
            betaHi *= 2.0;
            fBetaHi = erfc(betaHi * cutoff) / cutoff;
        }
    }
    // note erfc(beta*cutoff)/cutoff is monotone decreasing in beta
    double midBeta = 0.0;
    double diff = betaHi - betaLo;
    while (diff > 1.0e-15) {
        double midBeta = (betaLo + betaHi) / 2.0;
        double fMidBeta = erfc(midBeta * cutoff) / cutoff;
        if (fMidBeta > tol) {
            betaLo = midBeta;
        } else {
            betaHi = midBeta;
        }
        diff = betaHi - betaLo;
    }
    return (betaLo + betaHi) / 2.0;
}

void GetNandNprime(std::size_t &N, std::size_t &Nprime, double ewaldCoeff, double THRESH, double cutoff) {
    // formulas from
    /*
    "Resolutions of the Coulomb Operator: VII.
          Evaluation of Long-Range Coulomb and Exchange Matrices"
    Taweetham Limpanuparb, Josh Milthorpe, Alistair P. Rendell, and Peter M. W. Gill
     JCTC 9, 863-867 Jan 8 2013
    */

    double tol = pow(10.0, -THRESH);
    double r12Max = 4.0 * sqrt(3.0) * cutoff;
    double x = (M_PI / (4.0 * ewaldCoeff)) * tol;
    double omegaR = ewaldCoeff * r12Max;
    double rN = ceil(omegaR * omegaR + 2.0 * omegaR * (sqrt(-log10(x)) - 1.0) + 2.0);
    N = static_cast<std::size_t>(rN);
    double rNprime = ceil((2.0 / M_PI) * sqrt(-(rN + 1.0) * log(x)) - 1.0);
    Nprime = static_cast<std::size_t>(rNprime);
    if (Nprime > N) Nprime = N;
}

bool GetTrimmedGaussHermiteRoots_and_Weights(std::vector<double> &trimmedRoots, std::vector<double> &trimmedWeights,
                                             double &ewaldCoeff, double THRESH, double cutoff) {
    trimmedRoots.clear();
    trimmedWeights.clear();

    ewaldCoeff = GetEwaldCoeff(THRESH, cutoff);

    std::size_t N;
    std::size_t Nprime;

    GetNandNprime(N, Nprime, ewaldCoeff, THRESH, cutoff);
    std::cout << "N,Nprime = " << N << ", " << Nprime << std::endl;
    std::size_t n = 2 * N;

    std::vector<double> initialGuesses;
    if (!GetInitialGuessesForHermiteRoots(initialGuesses, n)) return false;

    std::vector<double> roots;
    std::vector<double> weights;
    GetRoots_and_WeightsOfHermite(roots, weights, n, initialGuesses);
    std::cout << "size of roots = " << roots.size() << std::endl;

    trimmedRoots.resize(Nprime, 0.0);
    trimmedWeights.resize(Nprime, 0.0);
    std::size_t istart = N - Nprime;
    for (std::size_t i = 0; i < Nprime; ++i) {
        trimmedRoots[i] = roots[istart + i];
        trimmedWeights[i] = weights[istart + i];
    }

    return true;
}

unsigned GetDegreeOfSphericalDesignFromHermiteRootsWeights(double THRESH, double directSumCutoff, double ewaldCoeff,
                                                           const std::vector<double> &gaussHermiteRoots,
                                                           const std::vector<double> &gaussHermiteWeights) {
    std::size_t Nprime = gaussHermiteRoots.size();
    double r12Max = 4.0 * sqrt(3.0) * directSumCutoff;  // max distance in 4x4x4 cubes of size cutoff
    double tol = pow(10.0, -THRESH) / Nprime;
    double lhs = 1.0;
    unsigned degree = 0;
    // std::cout << "r12Max, tol, N = " << r12Max << ", " << tol << ", " << N << std::endl;
    while (lhs > tol) {
        double max = 0.0;
        degree += 1;
        for (std::size_t n = 0; n < Nprime; ++n) {
            double lambda_n = 2.0 * ewaldCoeff * gaussHermiteRoots[n];
            double q_n = 4.0 * sqrt(ewaldCoeff * gaussHermiteWeights[n]);
            double sphBess = SphericalBessel(degree, lambda_n * r12Max);
            double compare = (q_n * sphBess) * (q_n * sphBess) / (4.0 * M_PI);
            if (compare > max) max = compare;
        }
        lhs = max;
    }
    unsigned degreeOfSpherePts = degree + 1;
    if (degreeOfSpherePts % 2 == 0) degreeOfSpherePts += 1;
    return degreeOfSpherePts;
}
