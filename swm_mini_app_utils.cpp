#include <string>
#include <cmath>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "swm_mini_app_utils.h"


void parse_input(int & nx, int & ny,
                 amrex::Real & dx, amrex::Real & dy,
                 int & max_chunk_size,
                 int & n_time_steps, amrex::Real & dt,
                 int & plot_interval)
{
    // ParmParse is way of reading inputs from the inputs file
    // pp.get means we require the inputs file to have it
    // pp.query means we optionally need the inputs file to have it - but we must supply a default here
    amrex::ParmParse pp;

    pp.get("nx",nx);
    pp.get("ny",ny);

    pp.get("dx",dx);
    pp.get("dy",dy);

    pp.get("max_chunk_size",max_chunk_size);

    pp.get("n_time_steps",n_time_steps);

    pp.get("dt",dt);

    // Default plot_interval to -1, allow us to set it to something else in the inputs file
    //  If plot_interval < 0 then no plot files will be written
    plot_interval = -1;
    pp.query("plot_interval",plot_interval);

}

void define_cell_centered_MultiFab(const int nx, const int ny,
                                   const int max_chunk_size,
                                   amrex::MultiFab & cell_centered_MultiFab)
{
    // lower and upper indices of domain
    const amrex::IntVect domain_low_index(0,0);
    const amrex::IntVect domain_high_index(nx-1, ny-1);
    
    // create box of indicies for cells
    const amrex::Box cell_centered_box(domain_low_index, domain_high_index);

    // initialize the boxarray "cell_box_array" from the single box "cell_centered_box"
    amrex::BoxArray cell_box_array(cell_centered_box);
    //cell_box_array.define(cell_centered_box);

    // break up boxarray "cell_box_array" into chunks no larger than "max_chunk_size" along a direction
    cell_box_array.maxSize(max_chunk_size);

    // assigns processor to each box in the box array
    amrex::DistributionMapping distribution_mapping(cell_box_array);

    // number of components for each array
    int Ncomp = 1;

    // number of ghost cells for each array
    int Nghost = 1;

    cell_centered_MultiFab.define(cell_box_array, distribution_mapping, Ncomp, Nghost);
}

void initialize_geometry(const int nx, const int ny,
                         const amrex::Real dx, const amrex::Real dy,
                         amrex::Geometry & geom)
{
  // lower and upper indices of domain
  const amrex::IntVect domain_low_index(0,0);
  const amrex::IntVect domain_high_index(nx-1, ny-1);

  // create box of indicies for cells
  const amrex::Box cell_centered_box(domain_low_index, domain_high_index);

  // physical min and max boundaries of cells
  const amrex::RealBox real_box({0, 0},
                                {nx*dx, ny*dy});

  // This, a value of 0, says we are using Cartesian coordinates
  int coordinate_system = 0;

  // This sets the boundary conditions in each direction to periodic
  amrex::Array<int,AMREX_SPACEDIM> is_periodic {1,1};

  // This defines a Geometry object
  geom.define(cell_centered_box, real_box, coordinate_system, is_periodic);
 // geom.define(cell_centered_box, real_box, amrex::CoordSys::cartesian, is_periodic);

  //amrex::Print() << "geom " << geom << std::endl;

  return;
}

amrex::Real linear_map_coordinates(const amrex::Real x, 
                                   const amrex::Real x_min, const amrex::Real x_max,
                                   const amrex::Real xi_min, const amrex::Real xi_max)
{
    return x_min + ((xi_max-xi_min)/(x_max-x_min))*x;
}

void initialize_variables(const amrex::Geometry & geom,
                          amrex::MultiFab & psi,
                          amrex::MultiFab & p,
                          amrex::MultiFab & u,
                          amrex::MultiFab & v)
{

    const amrex::Real x_min = geom.ProbLo(0);
    const amrex::Real x_max = geom.ProbHi(0);
    const amrex::Real y_min = geom.ProbLo(1);
    const amrex::Real y_max = geom.ProbHi(1);

    const amrex::Real dx = geom.CellSize(0);
    const amrex::Real dy = geom.CellSize(1);

    ////////////////////////////////////////////////////////////////////////// 
    // Initialization of stream function (psi)
    ////////////////////////////////////////////////////////////////////////// 

    // coefficient for initialization psi
    const amrex::Real a = 1000000;

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& phi_array = psi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const amrex::Real x_cell_center = (i+0.5) * dx;
            const amrex::Real y_cell_center = (j+0.5) * dy;

            const amrex::Real x_transformed = linear_map_coordinates(x_cell_center, x_min, x_max, 0.0, 2*std::numbers::pi);
            const amrex::Real y_transformed = linear_map_coordinates(y_cell_center, y_min, y_max, 0.0, 2*std::numbers::pi);

            phi_array(i,j,k) = a*std::sin(x_transformed)*std::sin(y_transformed);
        });
    }
    
    psi.FillBoundary(geom.periodicity());

    ////////////////////////////////////////////////////////////////////////// 
    // Initialization of pressure (p)
    ////////////////////////////////////////////////////////////////////////// 

    // coefficient for pressure
    double el = geom.ProbLength(0);
    amrex::Real pcf = (std::numbers::pi * std::numbers::pi * a * a)/(el * el);

    for (amrex::MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& p_array = p.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            amrex::Real x_node = i * dx;
            amrex::Real y_node = j * dy;
            
            const amrex::Real x_transformed = linear_map_coordinates(x_node, x_min, x_max, 0.0, 2*std::numbers::pi);
            const amrex::Real y_transformed = linear_map_coordinates(y_node, y_min, y_max, 0.0, 2*std::numbers::pi);

            p_array(i,j,k) = pcf * (std::cos(2*x_transformed) + std::cos(2*y_transformed)) + 5000;
        });
    }

    p.FillBoundary(geom.periodicity());

    ////////////////////////////////////////////////////////////////////////// 
    // Initialization of x velocity (u)
    ////////////////////////////////////////////////////////////////////////// 

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& u_array = u.array(mfi);
        const amrex::Array4<amrex::Real const>& phi_old_array = psi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            u_array(i,j,k) = -(phi_old_array(i,j,k)-phi_old_array(i,j-1,k))/dy;
        });
    }

    u.FillBoundary(geom.periodicity());

    ////////////////////////////////////////////////////////////////////////// 
    // Initialization of y velocity (v)
    ////////////////////////////////////////////////////////////////////////// 

    for (amrex::MFIter mfi(v); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& v_array = v.array(mfi);
        const amrex::Array4<amrex::Real const>& phi_old_array = psi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            v_array(i,j,k) = (phi_old_array(i,j,k)-phi_old_array(i-1,j,k))/dx;
        });
    }

    v.FillBoundary(geom.periodicity());

    return;
}


void write_output(const amrex::MultiFab & psi,
                  const amrex::MultiFab & p,
                  const amrex::MultiFab & u,
                  const amrex::MultiFab & v,
                  const amrex::Geometry & geom,
                  const amrex::Real time,
                  const int time_step,
                  amrex::MultiFab & output_values)
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
