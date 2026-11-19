#include "stdafx.h"
#include "Solver/Autorefine.h"
#include <cstring>
#include <new>
#include <iostream>

class TestRefiner : public CAutoRefine {
public:
    unsigned long BlockCount() const { return m_ulBlocksNum; }
};

int main()
{
    // Allocator reuse must not determine whether an empty refiner has disk blocks.
    alignas(TestRefiner) unsigned char storage[sizeof(TestRefiner)];
    std::memset(storage, 0xff, sizeof(storage));
    TestRefiner *refiner = new (storage) TestRefiner;
    if(refiner->BlockCount() != 0 || refiner->GetInteractionLevel() != AUTOREFINE_HIER_PRE_0_LEVEL) {
        std::cerr << "Fresh refiner has invalid cleanup state\n";
        return 1;
    }
    refiner->Clean(AUTOREFINE_DEALLMEM_ALL);
    refiner->~TestRefiner();
    return 0;
}
