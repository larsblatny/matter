#ifndef BLOCK_SCAN_GRID_HPP
#define BLOCK_SCAN_GRID_HPP

#include <vector>
#include <atomic>
#include <cstdint>
#include <cassert>
#include <cmath>

class BlockScanGrid {
public:
    struct I3 {int x, y, z;};

    BlockScanGrid(){};
    ~BlockScanGrid(){};

    // Initialize/Reinitialize the sparse grid 
    void init(int B_, I3 bmin_, I3 bmax_, int n_threads_) {
        assert(B_ > 0);
    }
};

#endif