#include <cmath>
#include <cstdint>
#include <iostream>

// boilerplate stuff from openmc
double prn(uint64_t* seed)
{
  constexpr uint64_t prn_mult {6364136223846793005ULL}; // multiplication
  constexpr uint64_t prn_add {1442695040888963407ULL};  // additive factor, c
  constexpr uint64_t prn_stride {152917LL}; // stride between particles
  // Advance the LCG
  *seed = (prn_mult * (*seed) + prn_add);

  // Permute the output
  uint64_t word =
    ((*seed >> ((*seed >> 59u) + 5u)) ^ *seed) * 12605985483714917081ull;
  uint64_t result = (word >> 43u) ^ word;

  // Convert output from unsigned integer to double
  return ldexp(result, -64);
}
double normal_variate(double mean, double standard_deviation, uint64_t* seed)
{
  // Sample a normal variate using Marsaglia's polar method
  double x, y, r2;
  do {
    x = 2.0 * prn(seed) - 1.0;
    y = 2.0 * prn(seed) - 1.0;
    r2 = x * x + y * y;
  } while (r2 > 1 || r2 == 0);
  double z = std::sqrt(-2.0 * std::log(r2) / r2);
  return mean + standard_deviation * z * x;
}
double normal_percentile(double p)
{
  constexpr double p_low = 0.02425;
  constexpr double a[6] = {-3.969683028665376e1, 2.209460984245205e2,
    -2.759285104469687e2, 1.383577518672690e2, -3.066479806614716e1,
    2.506628277459239e0};
  constexpr double b[5] = {-5.447609879822406e1, 1.615858368580409e2,
    -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1};
  constexpr double c[6] = {-7.784894002430293e-3, -3.223964580411365e-1,
    -2.400758277161838, -2.549732539343734, 4.374664141464968,
    2.938163982698783};
  constexpr double d[4] = {7.784695709041462e-3, 3.224671290700398e-1,
    2.445134137142996, 3.754408661907416};

  // The rational approximation used here is from an unpublished work at
  // http://home.online.no/~pjacklam/notes/invnorm/

  double z;
  double q;

  if (p < p_low) {
    // Rational approximation for lower region.

    q = std::sqrt(-2.0 * std::log(p));
    z = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);

  } else if (p <= 1.0 - p_low) {
    // Rational approximation for central region
    q = p - 0.5;
    double r = q * q;
    z = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) *
        q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);

  } else {
    // Rational approximation for upper region

    q = std::sqrt(-2.0 * std::log(1.0 - p));
    z = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }

  // Refinement based on Newton's method

  z = z - (0.5 * std::erfc(-z / std::sqrt(2.0)) - p) * std::sqrt(2.0 * M_PI) *
            std::exp(0.5 * z * z);

  return z;
}

/*
 * Samples x from the distribution proportional to:
 *
 * erf(sqrt(x) + a) - erf(sqrt(x) - a)
 *
 * where a>0.
 *
 * There are two ways to do this. The CDF inversion
 * technique simplifies considerably when a>4 because
 * exp(-a^2) ~= 0 and erf(a) ~= 1.
 *
 * For a<4, we use an expansion of the distribution
 * into Erlang random variables. Our paper shows how to
 * obtain that formula.
 *
 * seed - a uint64 reference for random number generation.
 */
double sample_erfs(double a, uint64_t& seed) {
  constexpr double rt_pi = std::sqrt(M_PI);
  double y;
  if (a < 4.0) {

    double prob = 0.0; // cumulative discrete prob
    double gamma = rt_pi; // sequence of half-index gamma fcn
    double p = std::erf(a);

    double xi = prn(&seed); // product of IID uniform r.v.
    double ctr = 0.5;
    double a_pow = a;

    // used for sampling the discrete sum
    const double expaa = std::exp(-a*a);
    const double xi_sum = prn(&seed) * (a / rt_pi * expaa + p * (0.5 + a*a));
    while (prob += p, prob < xi_sum) {
      gamma *= ctr;
      ctr += 1.0;
      p -= a_pow / gamma * expaa;
      a_pow *= a*a;
      xi *= prn(&seed);
    }
    y = -std::log(xi);

  } else {
    // Newton-Halley converges unconditionally and quickly.
    const double xi = prn(&seed);

    double corr = 0.0;
    y = 1e-9; // avoid division by zero
    const double nc = 1.0 + 2.0 * a*a; // normalizing constant
    do {
      const double rt_y = std::sqrt(y);
      const double expt = std::exp(-std::pow(rt_y - a, 2));
      const double erft = std::erf(rt_y - a);

      double f = (0.5 * nc - expt * (a + rt_y) / rt_pi + y + (0.5*nc - y) *
          erft) / nc - xi;
      double fp = (1.0 - erft) / nc;
      double fpp = -expt / (rt_pi * nc) / rt_y;
      corr = f / (fp - 0*f/fp*fpp*0.5);
      y -= corr;
    } while (std::abs(corr) > 1e-7);
  }
  return y;
}

// Samples an inverse Gaussian random variable
// The technique is standard.
double sample_ig(double mu, double lam, uint64_t* seed) {
  double w = mu * std::pow(normal_variate(0.0, 1.0, seed), 2);
  double c = 0.5 * mu / lam;
  double x1 = mu + c * (w - std::sqrt(w*(4*lam+w)));
  double x = x1;
  if (prn(seed) >= mu / (mu + x1)) {
    x = mu * mu / x1;
  }
  return x;
}


/*
 * Samples from the distribution proportional to:
 * exp(-lambda x) (erf(sqrt(x)+a) - erf(sqrt(x)-a))
 *
 * NOTE: this assumes that lambda > 0. It is wrong
 * in the case where lambda=0, which requires a
 * slightly different expression in the normalizing
 * constant. Assumes a>0.
 */
extern "C" double sample_erfs_lam(double a, double lam, uint64_t seed) {
  constexpr double rt_pi = std::sqrt(M_PI);
  double y;

  double prob = 0.0; // cumulative discrete prob
  double gamma = rt_pi; // sequence of half-index gamma fcn
  double p = std::erf(a);

  double xi = prn(&seed); // product of IID uniform r.v.
  double ctr = 0.5;
  double a_pow = a;

  // used for sampling the discrete sum
  const double expaa = std::exp(-a*a);
  const double rt_opl = std::sqrt(1.0 + lam);
  const double xi_sum = prn(&seed) / lam * (erf(a) -
      std::exp(-a*a*lam/(1.0+lam)) * erf(a/rt_opl)/rt_opl);
  double pow_lam = 1.0+lam;
  while (prob += p/pow_lam, prob < xi_sum) {
    pow_lam *= (1.0 + lam);
    gamma *= ctr;
    ctr += 1.0;
    p -= a_pow / gamma * expaa;
    a_pow *= a*a;
    xi *= prn(&seed);
  }
  y = -std::log(xi) / (1.0 + lam);

  return y;
}

/*
 * Samples from the Egelstaff-Schofield law using exp(-lambda alpha) as the
 * proposal distribution, exact sampling from the beta distribution without
 * kinematic boundaries, then rejects on valid kinematic boundary to obtain
 * the correct overall distribution.
 *
 * A - ratio of target mass to neutron mass
 * Ein - incident energy of the neutron in eV
 * T - temperature of the scattering material in kelvin
 * wt - weight of the Egelstaff-Schofield law. See NJOY2016 manual for detail
 * c - Egelstaff-Schofield diffusion coefficient
 * lambda - Debye-Waller coefficient (not to be confused with W)
 * seed - uint64 for random number generation
 *
 * alpha - where to save sampled alpha to
 * beta - where to save sampled beta to
 */
extern "C" void sample_egelstaff_schofield_rejection(double A, double Ein, double T, double wt,
    double c, double lambda, uint64_t seed, double* alpha, double* beta) {

  constexpr double kb = 8.617e-5;
  const double kbT = kb * T;

  bool reject = true;
  double r;
  do {
    // Proposal distribution for alpha
    *alpha = -std::log(prn(&seed)) / lambda;

    // Now sample the conditional distribution of beta over the full range
    const double delta = 2.0 * c * wt;
    const double alpha_ninvg = (*alpha) * std::sqrt(c*c + 0.25);
    const double beta_ninvg = -(*alpha) * 0.5;
    const double z = sample_ig(delta/std::sqrt(alpha_ninvg*alpha_ninvg-beta_ninvg*beta_ninvg),
        delta*delta, &seed);

    // Sample from NIG distribution
    r = std::sqrt(z)*normal_percentile(prn(&seed)) + beta_ninvg * z;

    // Check if the bounds are valid
    const double rt_alphahat = std::sqrt(A * kbT * (*alpha) / Ein);
    const double rmin = A * (1.0 - 2.0/rt_alphahat); // min beta/alpha ratio
    const double rmax = A * (1.0 + 2.0/rt_alphahat); // max beta/alpha ratio
    if (rmin < r && r < rmax) reject = false;
  } while (reject);
  *beta = r * (*alpha);
}

/*
 * Note: this samples beta/alpha in an unbounded fashion for use with the convolved
 * scattering law, which rejection has to be done on.
 */
extern "C" void sample_egelstaff_schofield_given_alpha(double A, double Ein, double T, double wt,
    double c, uint64_t seed, double alpha, double* beta) {

  constexpr double kb = 8.617e-5;
  const double kbT = kb * T;

  // variance and mean of the normal mixture variable come from z, below
  const double delta = 2.0 * c * wt;
  const double alpha_ninvg = alpha * std::sqrt(c*c + 0.25);
  const double beta_ninvg = -alpha * 0.5;
  const double z = sample_ig(delta/std::sqrt(alpha_ninvg*alpha_ninvg-beta_ninvg*beta_ninvg),
      delta*delta, &seed);

  // find the (unbounded) sample of beta/alpha ratio
  const double r = std::sqrt(z)*normal_percentile(prn(&seed)) + beta_ninvg * z;

  *beta = r * (alpha);
}

/*
 * After a value of alpha has been sampled, this samples a value of beta that is guaranteed
 * to lie within the kinematically valid boundaries without resort to any rejection sampling.
 */
extern "C" void sample_egelstaff_schofield_given_alpha_bounded(double A, double Ein, double T, double wt,
    double c, uint64_t seed, double alpha, double beta_s, double* beta) {

  constexpr double kb = 8.617e-5;
  const double kbT = kb * T;

  // variance and mean of the normal mixture variable come from z, below
  const double delta = 2.0 * c * wt;
  const double alpha_ninvg = alpha * std::sqrt(c*c + 0.25);
  const double beta_ninvg = -alpha * 0.5;
  const double z = sample_ig(delta/std::sqrt(alpha_ninvg*alpha_ninvg-beta_ninvg*beta_ninvg),
      delta*delta, &seed);

  // find the sample of beta/alpha ratio
  const double rt_alphahat = std::sqrt(A * kbT * alpha / Ein);
  const double bmin = alpha * A * (1.0 - 2.0/rt_alphahat); // min beta
  const double bmax = alpha * A * (1.0 + 2.0/rt_alphahat); // max beta
  constexpr double rt2 = std::sqrt(2.0);
  const double rtz = std::sqrt(z);
  const double pct_min = 0.5 * (1.0 + std::erf(( (bmin-beta_s)/alpha - beta_ninvg * z)/(rt2*rtz)));
  const double pct_max = 0.5 * (1.0 + std::erf(( (bmax-beta_s)/alpha- beta_ninvg * z)/(rt2*rtz)));
  const double r = rtz*normal_percentile(pct_min + prn(&seed) * (pct_max - pct_min)) + beta_ninvg * z;

  *beta = r * alpha;
}

// Samples only the E&S law when independently on its own, i.e. without any
// solid contributions. This is for the quasielastic term in the scattering
// law, and it forces the sampled value of beta to fall within the valid range.
extern "C" void sample_egelstaff_schofield(double A, double Ein, double T, double wt,
    double c, double lambda, uint64_t seed, double* alpha, double* beta) {

  constexpr double kb = 8.617e-5;
  const double kbT = kb * T;

  // first sample alpha
  const double cc = 2.0 * std::sqrt(Ein * A / kbT);
  const double eta = 2.0 * std::sqrt((1.0 + 0.25/(c*c)) * wt);
  const double scale = std::pow((wt + A)/eta, 2);
  *alpha = sample_erfs_lam(cc/eta, lambda/scale, seed) / scale;

  // variance and mean of the normal mixture variable come from z, below
  const double delta = 2.0 * c * wt;
  const double alpha_ninvg = *alpha * std::sqrt(c*c + 0.25);
  const double beta_ninvg = -*alpha * 0.5;
  const double z = sample_ig(delta/std::sqrt(alpha_ninvg*alpha_ninvg-beta_ninvg*beta_ninvg),
      delta*delta, &seed);

  // Find the percentiles associated with the allowable min and max beta/alpha
  // ratio used to sample a bounded normal inverse gaussian random variable

  // find the sample of beta/alpha ratio
  // BOUNDED VERSION
  const double rt_alphahat = std::sqrt(A * kbT * (*alpha) / Ein);
  const double rmin = A * (1.0 - 2.0/rt_alphahat); // min beta/alpha ratio
  const double rmax = A * (1.0 + 2.0/rt_alphahat); // max beta/alpha ratio
  constexpr double rt2 = std::sqrt(2.0);
  const double rtz = std::sqrt(z);
  const double pct_min = 0.5 * (1.0 + std::erf((rmin - beta_ninvg * z)/(rt2*rtz)));
  const double pct_max = 0.5 * (1.0 + std::erf((rmax - beta_ninvg * z)/(rt2*rtz)));
  const double r = rtz*normal_percentile(pct_min + prn(&seed) * (pct_max - pct_min)) + beta_ninvg * z;

  *beta = r * (*alpha);
}
