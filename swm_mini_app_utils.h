#ifndef SWM_MINI_APP_UTILS_H_
#define SWM_MINI_APP_UTILS_H_

#include <AMReX.H>

void write_output( amrex::MultiFab & output_values,
                   const amrex::MultiFab & psi,
                   const amrex::MultiFab & p,
                   const amrex::MultiFab & u,
                   const amrex::MultiFab & v,
                   const amrex::Geometry geom,
                   const amrex::Real time,
                   const int time_step);

#endif // SWM_MINI_APP_UTILS_H_