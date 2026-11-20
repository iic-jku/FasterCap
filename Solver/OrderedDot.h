#ifndef FASTERCAP_ORDERED_DOT_H
#define FASTERCAP_ORDERED_DOT_H

#include <cmath>
#include <cstddef>

inline void AccumulateOrderedDot(double value, double &sum, double &correction)
{
	double next = sum + value;
	correction += std::fabs(sum) >= std::fabs(value) ?
		(sum-next)+value : (value-next)+sum;
	sum = next;
}

// Fixed lanes allow independent arithmetic without making the reduction depend
// on buffer alignment, SIMD width or worker count. Compensation retains small
// terms during cancellation. Nonfinite terms propagate to a nonfinite result.
inline double OrderedDotProduct(const double *a, const double *b, std::size_t size)
{
	double sums[4] = {}, corrections[4] = {};
	std::size_t i = 0;
	for(; size-i>=4; i+=4) {
		for(std::size_t lane=0; lane<4; ++lane) {
			AccumulateOrderedDot(a[i+lane]*b[i+lane], sums[lane], corrections[lane]);
		}
	}
	double sum = 0.0, correction = 0.0;
	for(std::size_t lane=0; lane<4; ++lane) {
		AccumulateOrderedDot(sums[lane], sum, correction);
		AccumulateOrderedDot(corrections[lane], sum, correction);
	}
	for(; i<size; ++i) AccumulateOrderedDot(a[i]*b[i], sum, correction);
	return sum + correction;
}

#endif
