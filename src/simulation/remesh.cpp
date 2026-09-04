// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"

void Simulation::resizeGrid(){
    grid.v.resize(grid_nodes);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(grid_nodes); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.force.resize(grid_nodes); std::fill( grid.force.begin(), grid.force.end(), TV::Zero() );
    grid.mass.resize(grid_nodes); std::fill( grid.mass.begin(), grid.mass.end(), 0.0        );
    if (use_mibf) {
        grid.friction.resize(grid_nodes); 
        std::fill( grid.friction.begin(), grid.friction.end(), 0.0 );
    }
}

void Simulation::resetSparseGridScan() {
    grid.v.resize(scan_sparse_grid.num_active_nodes()); std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(scan_sparse_grid.num_active_nodes()); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.force.resize(scan_sparse_grid.num_active_nodes()); std::fill( grid.force.begin(), grid.force.end(), TV::Zero() );
    grid.mass.resize(scan_sparse_grid.num_active_nodes()); std::fill( grid.mass.begin(), grid.mass.end(), 0.0        );
    if (use_mibf) {
        grid.friction.resize(scan_sparse_grid.num_active_nodes()); std::fill( grid.friction.begin(), grid.friction.end(), 0.0 );
    }
}

// Get particles min and max node index
void Simulation::getParticlesMinMax() {
    const int node_lo_offset = 0;
    const int node_hi_offset = 2;

    std::vector<int> tminx(n_threads, INT32_MAX);
    std::vector<int> tminy(n_threads, INT32_MAX);
    std::vector<int> tminz(n_threads, INT32_MAX);
    std::vector<int> tmaxx(n_threads, -INT32_MAX);
    std::vector<int> tmaxy(n_threads, -INT32_MAX);
    std::vector<int> tmaxz(n_threads, -INT32_MAX);

    #pragma omp parallel num_threads(n_threads)
    {
        int tid = omp_get_thread_num();
        int lminx = INT32_MAX;
        int lminy = INT32_MAX;
        int lminz = INT32_MAX;
        int lmaxx = -INT32_MAX;
        int lmaxy = -INT32_MAX;
        int lmaxz = -INT32_MAX;

        #pragma omp for nowait
        for (int p=0; p<Np; p++) {
            TV xp = particles.x[p];
            int i_base = std::floor(xp(0)*one_over_dx - 0.5); 
            int j_base = std::floor(xp(1)*one_over_dx - 0.5); 
            int k_base = std::floor(xp(2)*one_over_dx - 0.5);

            // get local min/max base
            int ilo = i_base + node_lo_offset;
            int ihi = i_base + node_hi_offset;
            int jlo = j_base + node_lo_offset;
            int jhi = j_base + node_hi_offset;
            int klo = k_base + node_lo_offset;
            int khi = k_base + node_hi_offset;

            lminx = std::min(lminx, ilo);
            lminy = std::min(lminy, jlo);
            lminz = std::min(lminz, klo);
            lmaxx = std::max(lmaxx, ihi);
            lmaxy = std::max(lmaxy, jhi);
            lmaxz = std::max(lmaxz, khi);
        }

        tminx[tid] = lminx;
        tminy[tid] = lminy;
        tminz[tid] = lminz;
        tmaxx[tid] = lmaxx;
        tmaxy[tid] = lmaxy;
        tmaxz[tid] = lmaxz;
    }

    particles.minx_id = *std::min_element(tminx.begin(), tminx.end());
    particles.miny_id = *std::min_element(tminy.begin(), tminy.end());
    particles.minz_id = *std::min_element(tminz.begin(), tminz.end());
    particles.maxx_id = *std::max_element(tmaxx.begin(), tmaxx.end());
    particles.maxy_id = *std::max_element(tmaxy.begin(), tmaxy.end());
    particles.maxz_id = *std::max_element(tmaxz.begin(), tmaxz.end());
}

// Mark active blocks
void Simulation::markActiveBlocksScan() {
    scan_sparse_grid.clear_active();
    const int B = scan_sparse_grid.get_B();

    #pragma omp parallel for num_threads(n_threads)
    for (int p=0; p<Np; p++) {
        TV xp = particles.x[p];
        int i_base = std::floor(xp(0)*one_over_dx - 0.5); 
        int j_base = std::floor(xp(1)*one_over_dx - 0.5); 
        int k_base = std::floor(xp(2)*one_over_dx - 0.5);

        // loop on 3x3x3 stencil for quadratic B-splines
        for (int i=i_base; i<i_base+3; ++i)
        for (int j=j_base; j<j_base+3; ++j)
        for (int k=k_base; k<k_base+3; ++k) {
            int bx = BlockScanGrid::floor_div(i, B);
            int by = BlockScanGrid::floor_div(j, B);
            int bz = BlockScanGrid::floor_div(k, B);
            scan_sparse_grid.mark_block_active(bx, by, bz);
        }
    } // end particles loop
    
    scan_sparse_grid.rebuild_sparse_grid();
}

// A fixed grid - must hard-coded for every simulation
void Simulation::remeshFixed(unsigned int extra_nodes){

    grid.x = arange(-dx*(1+extra_nodes), Lx+(2+extra_nodes)*dx, dx);
    grid.y = arange(-dx,                 Ly+(2+extra_nodes)*dx, dx);
    #ifdef THREEDIM
        grid.z = arange(-dx*(1+extra_nodes), Lz+(2+extra_nodes)*dx, dx);
    #endif

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];
    #ifdef THREEDIM
        grid.zc = grid.z[0];
    #endif

    Nx = grid.x.size();
    Ny = grid.y.size();
    #ifdef THREEDIM
        Nz = grid.z.size();
        grid_nodes = Nx*Ny*Nz;
    #else
        grid_nodes = Nx*Ny;
    #endif
}

void Simulation::remeshFixedInit(unsigned int sfx, unsigned int sfy, unsigned int sfz){

    // ACTUAL min and max position of particles
    auto max_x_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T max_x = (*max_x_it)(0);
    auto max_y_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T max_y = (*max_y_it)(1);
#ifdef THREEDIM
    auto max_z_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T max_z = (*max_z_it)(2);
#endif
    auto min_x_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T min_x = std::floor((*min_x_it)(0)*one_over_dx) * dx; 
    auto min_y_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T min_y = std::floor((*min_y_it)(1)*one_over_dx) * dx; 
#ifdef THREEDIM
    auto min_z_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T min_z = std::floor((*min_z_it)(2)*one_over_dx) * dx; 
#endif
    
    // Save for remeshFixedCont
    max_x_init = max_x;
    min_x_init = min_x;

    max_y_init = max_y;
    min_y_init = min_y;

#ifdef THREEDIM
    max_z_init = max_z;
    min_z_init = min_z;
#endif

    // safety_factor = 2 means we have a grid which has a grid point 2*dx from the boundary particle
    // Assuming a local approach, a grid point 2dx away from a particle will not influence this particle
    unsigned int safety_factor_x = sfx;
    unsigned int safety_factor_y = sfy;
    unsigned int safety_factor_z = sfz;

    low_x_init    = min_x - dx * safety_factor_x;
    high_x_init   = max_x + dx * safety_factor_x;
    low_y_init    = min_y - dx * safety_factor_y;
    high_y_init   = max_y + dx * safety_factor_y;
#ifdef THREEDIM
    low_z_init  = min_z - dx * safety_factor_z;
    high_z_init = max_z + dx * safety_factor_z;
#endif

    grid.x = arange(low_x_init, high_x_init+dx, dx);
    grid.y = arange(low_y_init, high_y_init+dx, dx);
#ifdef THREEDIM
    grid.z = arange(low_z_init, high_z_init+dx, dx);
#endif

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];
#ifdef THREEDIM
    grid.zc = grid.z[0];
#endif

    Nx      = grid.x.size();
    Ny      = grid.y.size();
    Nx_init = Nx;
    Ny_init = Ny;
#ifdef THREEDIM
    Nz      = grid.z.size();
    Nz_init = Nz;
    grid_nodes = Nx*Ny*Nz;
#else
    grid_nodes = Nx*Ny;
#endif

    #ifdef WARNINGS
        #ifdef THREEDIM
        debug("               grid        = (", Nx, ", ", Ny, ", ", Nz, ")"  );
        #else
        debug("               grid        = (", Nx, ", ", Ny, ")"  );
        #endif
        debug("               min_x       = ", min_x);
        debug("               max_x       = ", max_x);
        debug("               min_y       = ", min_y);
        debug("               max_y       = ", max_y);
        #ifdef THREEDIM
        debug("               min_z       = ", min_z);
        debug("               max_z       = ", max_z);
        debug("               high_z_init = ", high_z_init);
        debug("               low_z_init  = ", low_z_init);
        #endif
        debug("               grid.xc     = ", grid.xc);
        debug("               grid.yc     = ", grid.yc);
        debug("               len(grid.x) = ", grid.x.size());
        debug("               len(grid.y) = ", grid.y.size());
        debug("               Nx          = ", Nx);
        debug("               Ny          = ", Ny);
        debug("               dx (obs)    = ", grid.x[1]-grid.x[0]);
    #endif

}

void Simulation::remeshFixedCont(){

    auto max_x_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T max_x = (*max_x_it)(0);
    auto min_x_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T min_x = (*min_x_it)(0);

    auto max_y_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T max_y = (*max_y_it)(1);
    auto min_y_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T min_y = (*min_y_it)(1);
#ifdef THREEDIM
    auto max_z_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T max_z = (*max_z_it)(2);
    auto min_z_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T min_z = (*min_z_it)(2);
#endif


    T high_x;
    if (max_x < max_x_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(max_x_init-max_x)/dx - 1e-8*dx) );
        high_x = high_x_init - reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(max_x-max_x_init)/dx - 1e-8*dx) );
        high_x = high_x_init + expansion_factor * dx;
    }

    T low_x;
    if (min_x > min_x_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(min_x-min_x_init)/dx - 1e-8*dx) );
        low_x = low_x_init + reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(min_x_init-min_x)/dx - 1e-8*dx) );
        low_x = low_x_init - expansion_factor * dx;
    }

    T high_y;
    if (max_y < max_y_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(max_y_init-max_y)/dx - 1e-8*dx) );
        high_y = high_y_init - reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(max_y-max_y_init)/dx - 1e-8*dx) );
        high_y = high_y_init + expansion_factor * dx;
    }

    T low_y;
    if (min_y > min_y_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(min_y-min_y_init)/dx - 1e-8*dx) );
        low_y = low_y_init + reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(min_y_init-min_y)/dx - 1e-8*dx) );
        low_y = low_y_init - expansion_factor * dx;
    }
#ifdef THREEDIM
    T high_z;
    if (max_z < max_z_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(max_z_init-max_z)/dx - 1e-8*dx) );
        high_z = high_z_init - reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(max_z-max_z_init)/dx - 1e-8*dx) );
        high_z = high_z_init + expansion_factor * dx;
    }

    T low_z;
    if (min_z > min_z_init){
        unsigned int reduction_factor = std::floor( std::max(0.0,(min_z-min_z_init)/dx - 1e-8*dx) );
        low_z = low_z_init + reduction_factor * dx;
    } else{
        unsigned int expansion_factor = std::floor( std::max(0.0,(min_z_init-min_z)/dx - 1e-8*dx) );
        low_z = low_z_init - expansion_factor * dx;
    }
#endif

    grid.x = arange(low_x, high_x+dx, dx);
    grid.y = arange(low_y, high_y+dx, dx);
#ifdef THREEDIM
    grid.z = arange(low_z, high_z+dx, dx);
#endif

    Nx      = grid.x.size();
    Ny      = grid.y.size();
#ifdef THREEDIM
    Nz      = grid.z.size();
    grid_nodes = Nx*Ny*Nz;
#else
    grid_nodes = Nx*Ny;
#endif

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];
#ifdef THREEDIM
    grid.zc = grid.z[0];
#endif

    #ifdef WARNINGS
        #ifdef THREEDIM
        debug("               grid        = (", Nx, ", ", Ny, ", ", Nz, ")"  );
        #else
        debug("               grid        = (", Nx, ", ", Ny, ")"  );
        #endif
    #endif

}
