#include "OrderedDot.h"
#include <cstdio>
#include <limits>

int main()
{
    alignas(128) double a[1024], b[1024];
    // Exact sum is 480. BLAS may regroup these equal vectors differently at
    // different addresses; the observed aligned/unaligned results were 480/360.
    for(int offset=0; offset<32; ++offset) {
        for(int i=0; i<960; ++i) {
            a[i+offset] = i%4==0 ? 1e16 : i%4==2 ? -1e16 : 1.0;
            b[i+offset] = 1.0;
        }
        double value = OrderedDotProduct(a+offset,b+offset,960);
        if(value != 480.0) {
            std::printf("offset %d: %.17g, expected 480\n",offset,value);
            return 1;
        }
    }
    // Exact integer sums exercise empty, short and partial final groups.
    for(int size=0; size<33; ++size) {
        double expected = 0;
        for(int i=0; i<size; ++i) {
            a[i] = i%5-2; b[i] = i%3-1;
            expected += a[i]*b[i];
        }
        if(OrderedDotProduct(a,b,size) != expected) return 2;
    }
    for(int position=0; position<9; ++position) {
        for(int i=0; i<9; ++i) a[i] = b[i] = 1.0;
        a[position] = std::numeric_limits<double>::quiet_NaN();
        if(std::isfinite(OrderedDotProduct(a,b,9))) return 3;
        a[position] = std::numeric_limits<double>::infinity();
        if(std::isfinite(OrderedDotProduct(a,b,9))) return 4;
    }
    return 0;
}
