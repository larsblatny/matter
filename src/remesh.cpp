#include "simulation.hpp"

void Simulation::remeshFixed(){
    // A fixed grid hard-coded for ever simulation

    int extra_nodes = 0;
    grid.x = arange(-dx*(1+extra_nodes), dx+(2+extra_nodes)*dx, dx);
    grid.y = arange(-dx*(1+extra_nodes), Ly+(2+extra_nodes)*dx, dx);

    grid.xc = grid.x[0];
    grid.yc = grid.y[0];

    Nx = grid.x.size();
    Ny = grid.y.size();
    grid_nodes = Nx*Ny;

    grid.v.resize(grid_nodes);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(grid_nodes); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(grid_nodes); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );


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
    T min_x = (*min_x_it)(0);
    auto min_y_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T min_y = (*min_y_it)(1);
#ifdef THREEDIM
    auto min_z_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV &x1, const TV &x2 )
                                             {
                                                 return x1(2) < x2(2);
                                             } );
    T min_z = (*min_z_it)(2);
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
    unsigned int safety_factor_x = sfx; // std::max(sfx, nonlocal_support);
    unsigned int safety_factor_y = sfy; // std::max(sfy, nonlocal_support);
    unsigned int safety_factor_z = sfz; // std::max(sfz, nonlocal_support);

    low_x_init    = min_x - dx * safety_factor_x;
    high_x_init   = max_x + dx * safety_factor_x;
    low_y_init    = min_y - dx * safety_factor_y;
    high_y_init   = max_y + dx * safety_factor_y;
#ifdef THREEDIM
    T low_z_init  = min_z - dx * safety_factor_z;
    T high_z_init = max_z + dx * safety_factor_z;
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

    grid.v.resize(grid_nodes);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(grid_nodes); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(grid_nodes); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );

    #ifdef WARNINGS
        #ifdef THREEDIM
            debug("               grid   = (", Nx, ", ", Ny, ", ", Ny, ")"  );
        #else
            debug("               grid   = (", Nx, ", ", Ny, ")"  );
        #endif

        debug("               min_x   = ", min_x);
        debug("               max_x   = ", max_x);
        debug("               min_y   = ", min_y);
        debug("               max_y   = ", max_y);
        debug("               grid.xc = ", grid.xc);
        debug("               grid.yc = ", grid.yc);
        debug("               len(grid.x) = ", grid.x.size());
        debug("               len(grid.y) = ", grid.y.size());
        debug("               Nx      = ", Nx);
        debug("               Ny      = ", Ny);
        debug("               dx_meas = ", grid.x[1]-grid.x[0]);
        // std::vector<double> tmp = linspace(0.0, 1.0, 5);
        // for (double i: tmp)
        //     std::cout << i << ' ';
        // std::cout << std::endl;
        // std::vector<double> tmp = arange(0.0, 1.0, 0.2);
        // for (double i: tmp)
        //     std::cout << i << ' ';
        // std::cout << std::endl;
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

    grid.v.resize(grid_nodes);    std::fill( grid.v.begin(),    grid.v.end(),    TV::Zero() );
    grid.flip.resize(grid_nodes); std::fill( grid.flip.begin(), grid.flip.end(), TV::Zero() );
    grid.mass.resize(grid_nodes); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );

    #ifdef WARNINGS
        #ifdef THREEDIM
            debug("               grid   = (", Nx, ", ", Ny, ", ", Ny, ")"  );
        #else
            debug("               grid   = (", Nx, ", ", Ny, ")"  );
        #endif
    #endif

}
