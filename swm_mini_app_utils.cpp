#include <string>
#include <cmath>

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>

#include "swm_mini_app_utils.h"

amrex::Real linear_map_coordinates(const amrex::Real x, 
                                   const amrex::Real x_min, const amrex::Real x_max,
                                   const amrex::Real xi_min, const amrex::Real xi_max)
{
    return x_min + ((xi_max-xi_min)/(x_max-x_min))*x;
}

void initialize_psi(amrex::MultiFab & psi,
                    const amrex::Geometry geom)
{
    // coefficient for initialization of stream function
    const amrex::Real a = 1000000;

    const amrex::Real x_min = geom.ProbLo(0);
    const amrex::Real x_max = geom.ProbHi(0);
    const amrex::Real y_min = geom.ProbLo(1);
    const amrex::Real y_max = geom.ProbHi(1);

    // loop over cell centers
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& phi_array = psi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            // Three ways to get the cell center:

            // 1.) Using CellSizeArray and computing by hand         
            //amrex::GpuArray<amrex::Real,2> dx = geom.CellSizeArray();
            //const amrex::Real x_cell_center = (i+0.5) * dx[0];
            //const amrex::Real y_cell_center = (j+0.5) * dx[1];

            // 2.) Using CellSize and computing by hand         
            //const amrex::Real dx = geom.CellSize(0);
            //const amrex::Real dy = geom.CellSize(1);
            //const amrex::Real x_cell_center = (i+0.5) * dx;
            //const amrex::Real y_cell_center = (j+0.5) * dy;

            // 3.) Using built in CellCenter function
            const amrex::Real x_cell_center = geom.CellCenter(i, 0);
            const amrex::Real y_cell_center = geom.CellCenter(j, 1);

            const amrex::Real x_transformed = linear_map_coordinates(x_cell_center, x_min, x_max, 0.0, 2*std::numbers::pi);
            const amrex::Real y_transformed = linear_map_coordinates(y_cell_center, y_min, y_max, 0.0, 2*std::numbers::pi);

            phi_array(i,j,k) = a*std::sin(x_transformed)*std::sin(y_transformed);
        });
    }

    psi.FillBoundary(geom.periodicity());

    return;
}

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
