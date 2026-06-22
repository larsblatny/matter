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
        B = B_;
        nodes_per_block = B * B * B;
        n_threads = std::max(1, n_threads_);

        bmin = bmin_;
        bmax = bmax_;

        nbx = bmax.x - bmin.x;
        nby = bmax.y - bmin.y;
        nbz = bmax.z - bmin.z;

        const int n_blocks_dense = nbx * nby * nbz;

        block_active.assign(n_blocks_dense, 0u);
        block_id.assign(n_blocks_dense, -1);
        block_xyz_by_id.clear();
        block_xyz_by_id.shrink_to_fit();
        scan_sum.assign(n_blocks_dense, 0);

        n_active_blocks = 0;
    }

    int get_B() const {
        return B;
    }

    int num_active_nodes() const {
        return n_active_blocks * nodes_per_block;
    }

    // Clear active blocks
    void clear_active() {
        std::fill(block_active.begin(), block_active.end(), uint8_t(0));
    }

    // Mark a block active
    void mark_block_active(int bx, int by, int bz) {
        const int idx = block_flattened_index(bx, by, bz);
        if (idx < 0) {
            return;
        }

        #pragma omp atomic write // TODO: try to replace with normal writing without atomic operation
        block_active[idx] = 1u;
    }

    // After marking active blocks, rebuild sparse grid
    void rebuild_sparse_grid() {
        const int n_blocks_dense = block_active.size();
        if (n_blocks_dense==0) {
            n_active_blocks = 0;
            block_xyz_by_id.clear();
            return;
        }

        // 1. Scan
        exclusive_scan(block_active.data(), scan_sum.data(), n_blocks_dense, n_threads, n_active_blocks);

        // 2. Reinitialize block_xyz_by_id
        block_xyz_by_id.assign(n_active_blocks, I3{0, 0, 0});

        // 3. Fill block_id and block_xyz_by_id
        #pragma omp parallel for num_threads(n_threads)
        for (int i=0; i<n_blocks_dense; ++i) {
            if (block_active[i]) {
                const int id = scan_sum[i];
                block_id[i] = id;

                // i -> (bx, by, bz)
                I3 b = decode_flattened_index(i);
                block_xyz_by_id[id] = b;
            } else {
                block_id[i] = -1;
            }
        }
    }

    // Grid node (i, j, k) -> flattened node id
    int ijk_to_gid(int ix, int iy, int iz) const {
        int bx = floor_div(ix, B);
        int by = floor_div(iy, B);
        int bz = floor_div(iz, B);

        // local
        int lx = ix - bx * B;
        int ly = iy - by * B;
        int lz = iz - bz * B;

        const int lid = lx + B * (ly + B * lz);

        const int bidx = block_flattened_index(bx, by, bz);
        if (bidx < 0) {
            return -1;
        }

        const int bid = block_id[bidx];
        if (bid < 0) {
            return -1;
        }

        return bid * nodes_per_block + lid;
    }

    // Flattened node id -> node (i, j, k)
    I3 gid_to_ijk(int gid) const {
        const int bid = gid / nodes_per_block;
        const int lid = gid - bid * nodes_per_block;

        const int lx = lid % B;
        const int t = lid / B;
        const int ly = t % B;
        const int lz = t / B;

        const I3 b = block_xyz_by_id[bid];
        return I3{b.x*B+lx, b.y*B+ly, b.z*B+lz};
    }

    // Floor division supporting negative integers
    static int floor_div(int a, int b) {
        // for possible negative a
        int q = a / b;
        int r = a % b;
        if (r < 0) q -= 1; // for negative a
        return q;
    }

private:
    int B = 4; // block size
    int nodes_per_block = B * B * B;
    int n_threads = 1;

    I3 bmin{0, 0, 0}; // minimum block index
    I3 bmax{0, 0, 0}; // maximum block index

    int nbx = 0; // number of blocks in x
    int nby = 0; // number of blocks in y
    int nbz = 0; // number of blocks in z

    std::vector<uint8_t> block_active; // binary flag map
    std::vector<int> block_id; // block dense id -> compact id
    std::vector<I3> block_xyz_by_id; // compact id -> block spatial index
    std::vector<int> scan_sum; // prefix sum offset

    int n_active_blocks = 0;

    // (i, j, k) -> flattened block index
    int block_flattened_index(int bx, int by, int bz) const {
        if (!is_inside(bx, by, bz)) {
            return -1;
        }

        const int lx = bx - bmin.x;
        const int ly = by - bmin.y;
        const int lz = bz - bmin.z;
        return (lx * nby + ly) * nbz + lz;
    }

    // Flattened block index -> (i, j, k)
    I3 decode_flattened_index(int idx) const {
        // idx = (lx*nby + ly)*nbz + lz
        int lz = idx % nbz;
        int t = idx / nbz;
        int ly = t % nby;
        int lx = t / nby;
        return I3{lx + bmin.x, ly + bmin.y, lz + bmin.z};
    }

    // Check whether the block is inside the dense domain
    bool is_inside(int bx, int by, int bz) const {
        return (bx >= bmin.x && bx < bmax.x &&
                by >= bmin.y && by < bmax.y &&
                bz >= bmin.z && bz < bmax.z);
    }

    // Parallel exclusive scan
    static void exclusive_scan(const uint8_t* in, int* out, int n_blocks_dense, int n_threads, int& n_active_blocks) {
        if (n_blocks_dense==0) {
            n_active_blocks = 0;
            return;
        }
        n_threads = std::max(1, n_threads);
        std::vector<int> thread_sums(std::size_t(n_threads + 1), 0);

        // Thread-wise local scan
        #pragma omp parallel num_threads(n_threads)
        {
            int tid = omp_get_thread_num();
            int start = (n_blocks_dense * tid) / n_threads;
            int end = (n_blocks_dense * (tid + 1)) / n_threads;

            int thread_value = 0;
            for (int i=start; i<end; ++i) {
                out[i] = thread_value;
                thread_value += int(in[i] != 0);
            }
            thread_sums[tid+1] = thread_value;
        }

        // Prefix sum of thread_sums
        for (int t=1; t<=n_threads; ++t) {
            thread_sums[t] += thread_sums[t-1];
        }

        // Thread-wise offset
        #pragma omp parallel num_threads(n_threads)
        {
            int tid = omp_get_thread_num();
            int start = (n_blocks_dense * tid) / n_threads;
            int end = (n_blocks_dense * (tid + 1)) / n_threads;

            int offset = thread_sums[tid];
            for (int i=start; i<end; ++i) {
                out[i] += offset;
            }
        }

        // Get n_active_blocks from thread_sums
        n_active_blocks = thread_sums[n_threads];
    }
};

#endif