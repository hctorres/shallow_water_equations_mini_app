
#include <string>

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>

#include "swm_mini_app_utils.h"

void write_output( amrex::MultiFab & output_values,
                   const amrex::MultiFab & psi,
                   const amrex::MultiFab & p,
                   const amrex::MultiFab & u,
                   const amrex::MultiFab & v,
                   const amrex::Geometry geom,
                   const amrex::Real time,
                   const int time_step) 
{

    // Interpolate all values to cell centers
    for (amrex::MFIter mfi(output_values); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real const>& phi_old_array = psi.const_array(mfi);
        const amrex::Array4<amrex::Real const>& p_array = p.const_array(mfi);
        const amrex::Array4<amrex::Real const>& u_array = u.const_array(mfi);
        const amrex::Array4<amrex::Real const>& v_array = v.const_array(mfi);

        const amrex::Array4<amrex::Real>& output_values_array = output_values.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            output_values_array(i,j,k,0) = phi_old_array(i,j,k);
            output_values_array(i,j,k,1) = (p_array(i,j,k) + p_array(i+1,j,k) + p_array(i,j+1,k) + p_array(i+1,j+1,k))/4.0;
            output_values_array(i,j,k,2) = (u_array(i,j,k) + u_array(i,j+1,k))/2.0;
            output_values_array(i,j,k,3) = (v_array(i,j,k) + v_array(i+1,j,k))/2.0;
        });
    }

    // Write output file
    int min_digits = 5;
    const std::string& pltfile = amrex::Concatenate("plt", time_step, min_digits);
    amrex::WriteSingleLevelPlotfile(pltfile, output_values, {"psi", "p", "u", "v"}, geom, time, time_step);

    return;
}
             