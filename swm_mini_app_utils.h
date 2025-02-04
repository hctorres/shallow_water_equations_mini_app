#ifndef SWM_MINI_APP_UTILS_H_
#define SWM_MINI_APP_UTILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

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

void define_cell_centered_MultiFab(const int nx, const int ny,
                                   const int max_chunk_size,
                                   amrex::MultiFab & cell_centered_MultiFab);

void define_x_face_MultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                            amrex::MultiFab & x_face_MultiFab);

void define_y_face_MultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                            amrex::MultiFab & y_face_MultiFab);

void define_nodal_MultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                           amrex::MultiFab & nodal_MultiFab);

void initialize_variables(const amrex::Geometry & geom,
                          amrex::MultiFab & psi,
                          amrex::MultiFab & p,
                          amrex::MultiFab & u,
                          amrex::MultiFab & v);

void write_output(const amrex::MultiFab & psi,
                  const amrex::MultiFab & p,
                  const amrex::MultiFab & u,
                  const amrex::MultiFab & v,
                  const amrex::Geometry & geom,
                  const amrex::Real time,
                  const int time_step,
                  amrex::MultiFab & output_values);

amrex::MultiFab createMultiFab(const amrex::MultiFab & mf);

void Copy(const amrex::MultiFab & src, amrex::MultiFab & dest);

void computeIntermediateVariables(amrex::Real fsdx, amrex::Real fsdy, const amrex::Geometry& geom,
                                 const amrex::MultiFab& p, const amrex::MultiFab& u, const amrex::MultiFab& v,
                                 amrex::MultiFab& cu, amrex::MultiFab& cv, amrex::MultiFab& h, amrex::MultiFab& z);

void updateNewVariables(const double tdtsdx, const double tdtsdy, const double tdts8, const amrex::Geometry& geom,
                        const amrex::MultiFab& p_old, const amrex::MultiFab& u_old, const amrex::MultiFab& v_old,
                        const amrex::MultiFab& cu, const amrex::MultiFab& cv, const amrex::MultiFab& h, const amrex::MultiFab& z,
                        amrex::MultiFab& p_new, amrex::MultiFab& u_new, amrex::MultiFab& v_new);

#endif // SWM_MINI_APP_UTILS_H_