#ifndef SWM_MINI_APP_UTILS_H_
#define SWM_MINI_APP_UTILS_H_

#include <AMReX.H>

void parse_input(int & nx, int & ny,
                 amrex::Real & dx, amrex::Real & dy,
                 int & max_chunk_size,
                 int & n_time_steps, amrex::Real & dt,
                 int & plot_interval);


void initialize_geometry(const int nx, const int ny,
                         const amrex::Real dx, const amrex::Real dy,
                         amrex::Geometry & geom);


// Linear mapping of a value (x) from one interval [x_min, x_max] to another [xi_min, xi_max].
// Precondition: x_min <= x <= x_max
amrex::Real linear_map_coordinates(const amrex::Real x, 
                                   const amrex::Real x_min,  const amrex::Real x_max,
                                   const amrex::Real xi_min, const amrex::Real xi_max);

void initialize_variables(amrex::MultiFab & psi,
                          amrex::MultiFab & p,
                          amrex::MultiFab & u,
                          amrex::MultiFab & v,
                         const amrex::Geometry geom);

void write_output( amrex::MultiFab & output_values,
                   const amrex::MultiFab & psi,
                   const amrex::MultiFab & p,
                   const amrex::MultiFab & u,
                   const amrex::MultiFab & v,
                   const amrex::Geometry geom,
                   const amrex::Real time,
                   const int time_step);



#endif // SWM_MINI_APP_UTILS_H_