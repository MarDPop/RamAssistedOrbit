#ifndef LIB_MATH_H
#define LIB_MATH_H

#include <algorithm>
#include <cmath>
#include <cstint>
#include <limits>

#define _USE_MATH_DEFINES

constexpr double INV_LN2 = 1.0 / M_LN2;

inline sqrt_approx(float x)
{
	union 
	{
		float f;
		uint32_t i;
	} val = {x};
	
	val.i -= (1 << 23);
	val.i >>= 1;
	val.i += (1 << 29);
}

inline double sqrt_approx(double d)
{
	double x;
	asm (
		"movq %1, %%xmm0 ]n"
		"sqrtsd %%xmm0, %%xmm1 \n"
		"movq %%xmm1, %0 \n"
		: "=r" (x)
		: "g" (d)
		: "xmm0", "xmm1", "memory"
	);
	return x;
}

inline float fast_rsqrt(float x)
{
	union {
		float f;
		uint32_t i;
	} conv = {x};
	conv.i = 0x5f3759df - (conv.i >> 1);
	conv.f *= 1.5f - (x*0.5f*conv.f*conf.f);
	return conv.f;
}

inline double fast_arccos(double x)
{
	x = 1.0 - x;
	double y  = 0.007 + x*13.7;
	y = y*0.5 + x/y;
	return y = y*0.5 + x/y;
}

inline double fast_ln(float x)
{
	union {
		float f;
		uint32_t i;
	} conv = {x};
	
	int t = static_cast<int>((conv.i >> 23) - 127);
	conv.i = 1065353216 | (conv.i & 8388607);
	return -1.49278f + (2.11263f + conv.f*(conv.f*0.10969f - 0.729104f))*conv.f + M_LN2*t;
}

template<typename Float>
constexpr Float sqrt_impl(Float x, Float curr, Float prev)
{
	return (curr == prev) ? curr : sqrt_impl(x, 0.5*(curr + x /curr), curr);
}

template<typename Float>
constexpr Float sqrt(Float x)
{
	return x >= 0.0 ? sqrt_impl(x, x, 0.0) : x / 0.0;
}

template<typename Float>
constexpr Float pow(Float x, unsigned n)
{
	return n == 1u ? x : n == 0 ? 1.0 : pow(x, n - (n >> 1)) * pow(x, n >> 1);
}

template<typename Float>
constexpr unsigned get_exp_k(Float x)
{
	return static_cast<unsigned>(x*INV_LN2 + 0.5);
}

template<typename Float>
constexpr Float expm1_impl(Float x, Float curr, Float prev, Float x_n, Float n)
{
	return (curr == prev || n > 17.5) ? curr : expm1_impl(x, curr + x_n, curr, x_n*x/n, n + 1.0);
}

template<typename Float>
constexpr Float exp_n_impl(Float x, unsigned k)
{
	return (1.0 + expm1_impl(k * M_LN2 + x, -1.0, 0.0, 1.0, 1.0))/(1lu << k);
}

template<typename Float>
constexpr Float exp_p_impl(Float x, unsigned k)
{
	return (1lu << k) * (1.0 + expm1_impl(x - k*M_LN2, -1.0, 0.0, 1.0, 1.0));
}

template<typename Float>
constexpr Float exp(Float x)
{
	return (fabs(x) > 44.01484) ? std::numeric_limits<Float>::infinity() : (x >= 0.0 ? exp_p_impl(x, get_exp_k(x)): exp_n_impl(x, get_exp_k(-x)));
}

template<typename Float>
constexpr unsigned get_ln_k(Float x, Float ten, unsigned ex)
{
	return ten > x ? ex : get_ln_k(x, ten*static_cast<Float>(10.0), ex + 1);
}

template<typename Float>
constexpr Float ln_impl(Float x, Float curr, Float prev)
{
	return (curr == prev) ? curr : ln_impl(x, curr + 2.0*(x - exp(curr))/(x + exp(curr)), curr);
}

template<typename Float>
constexpr Float ln_inter(Float x, unsigned k)
{
	return ln_impl(x, pow(10.0, k), 1.0, 0.0) + k * M_LN10;
}

template<typename Float>
constexpr Float ln(Float x)
{
	return x > 0.0 : ln_inter(x, get_ln_k(x, static_cast<Float>(10.0), 0u)): x / 0.0;
}

template<typename Float>
constexpr Float exp_frac_approx(Float x)
{
	const Float x2 = x*x;
	return 1.0 + (x + x)/(2.0 - x + x2/(6.0 + x2*0.1));
}


namespace poly 
{
	template<typename Float>
	struct poly3_x_t 
	{
		Float* x2;
		Float invNumPts;
		Float sumX;
		Float sumX2;
		Float sumX3;
		Float sumX4;
		Float meanX;
		Float meanX2;
		Float S11;
		Float S12;
		Float S22;
		Float denDet;
		
		poly3_x_t(const Float* x, const unsigned n) : x2(new Float[n]), invNumPts(1.0/static_cast<Float>(n))
		{
			for(auto i = 0u; i < n; i++)
			{
				x2[i] = x[i]*x[i];
				sumX += x[i];
				sumX2 += x2[i];
				sumX3 += x2[i]*x[i];
				sumX4 += x2[i]*x2[i];
			}
			
			meanX = sumX * invNumPts;
			meanX2 = sumX2 * invNumPts;
			
			S11 = sumX2 - sumX*meanX;
			S12 = sumX3 - sumX*meanX2;
			S22 = sumX4 - sumX2*meanX2;
			
			denDet = 1.0/(S22*S11 - S12*S12);
		}
		
		~poly3_x_t()
		{
			delete[] x2;
		}
	};
	
	// highest coef first
	template<typename Float>
	void poly3(const Float* x, const poly3_x_t<Float>& x_const, const Float* y, const unsigned n, Float* coef)
	{
		Float sumY = 0.0;
		Float sumXY = 0.0;
		Float sumX2Y = 0.0;
		
		for(auto i = 0u; i < n; i++)
		{
			sumY += y[i];
			sumXY += x[i]*y[i];
			sumX2Y += x_const.x2[i]*y[i];
		}
		
		const Float Sy1 = sumXY - x_const.meanX*sumY;
		const Float Sy2 = sumX2Y - x_const.meanX2*sumY;
		
		coef[0] = (Sy2*x_const.S11 - Sy1*x_const.S12)*x_const.denDet;
		coef[1] = (Sy1*x_const.S22 - Sy2*x_const.S12)*x_const.denDet;
		
		coef[2] = (sumY*x_const.invNumPts - coef[0]*x_const.meanX2 - coef[1]*x_const.meanX);
	}
	
	template<typename Float>
	struct poly2_x_t 
	{
		Float invNumPts;
		Float sumX;
		Float invSumX2;
		
		poly2_x_t(const Float* x, const unsigned n) : invNumPts(1.0/static_cast<Float>(n))
		{
			for(auto i = 0u; i < n; i++)
			{
				meanX += x[i];
				sumX2 += x[i]*x[i];
			}
			invSumX2 = 1.0/invSumX2;
		}
	};
	
	void poly2(const Float* x, const poly2_x_t<Float>& x_const, const Float* y, const unsigned n, Float* coef)
	{
		Float sumY = 0.0;
		Float sumXY = 0.0;
		
		for(auto i = 0u; i < n; i++)
		{
			sumY += y[i];
			sumXY += x[i]*y[i];
		}
		
		coef[1] = sumXY * x_const.invSumX2;
		coef[0] = (sumY - coef[1])*x_const.meanX;
	}

	
}

#endif