/*
 *
 */

#include <numbers>
#include <cmath>

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include "swm_mini_app_utils.h"


int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {

    // ***********************************************************************
    // Simulation Parameters Set Via Input File
    // ***********************************************************************

    // number of cells on each direction
    int nx;
    int ny;

    // cell size in each direction
    amrex::Real dx;
    amrex::Real dy;

    // mesh will be broken into chunks of up to max_chunk_size
    int max_chunk_size;

    // number of time steps to take
    int n_time_steps;

    // size of time step
    amrex::Real dt;

    // how often to write a plotfile
    //     Optional argument. If left out no plot files will be written.
    int plot_interval;

    // Set parameter values from inputs file
    parse_input(nx, ny, dx, dy, max_chunk_size,
                n_time_steps, dt, plot_interval);

    // ***********************************************************************
    // Define simulation setup and geometry
    // ***********************************************************************

    // define lower and upper indices
    amrex::IntVect dom_lo(0,0);
    amrex::IntVect dom_hi(nx-1, ny-1);
    
    amrex::Box cell_centered_box(dom_lo, dom_hi);

    // Initialize the boxarray "cell_box_array" from the single box "domain"
    amrex::BoxArray cell_box_array;
    cell_box_array.define(cell_centered_box);

    // Break up boxarray "cell_box_array" into chunks no larger than "max_chunk_size" along a direction
    cell_box_array.maxSize(max_chunk_size);

    // Ncomp = number of components for each array
    int Ncomp = 1;

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    amrex::DistributionMapping distribution_mapping(cell_box_array);

    amrex::MultiFab psi(cell_box_array, distribution_mapping, Ncomp, Nghost);

    amrex::BoxArray x_face_box_array = amrex::convert(cell_box_array, {1,0});
    amrex::MultiFab v(x_face_box_array, distribution_mapping, Ncomp, Nghost); 

    amrex::BoxArray y_face_box_array = amrex::convert(cell_box_array, {0,1});
    amrex::MultiFab u(y_face_box_array, distribution_mapping, Ncomp, Nghost);  

    amrex::BoxArray surrounding_nodes_box_array = cell_box_array;
    surrounding_nodes_box_array.surroundingNodes();
    amrex::MultiFab p(surrounding_nodes_box_array, distribution_mapping, Ncomp, Nghost);

    //amrex::Print() << "distribution_mapping " << distribution_mapping << std::endl;
    //amrex::Print() << "x_face_box_array " << x_face_box_array << std::endl;
    //amrex::Print() << "y_face_box_array " << y_face_box_array << std::endl;
    //amrex::Print() << "surrounding_nodes_box_array " << surrounding_nodes_box_array << std::endl;

    amrex::Geometry geom;
    initialize_geometry(nx, ny, dx, dy, geom);
    
    // **********************************
    // Initialize Data
    // **********************************

    initialize_variables(geom, psi, p, u, v);

    amrex::Print() << "Initial: " << std::endl;
    amrex::Print() << "psi max: " << psi.max(0) << std::endl;
    amrex::Print() << "psi min: " << psi.min(0) << std::endl;
    amrex::Print() << "p max: " << p.max(0) << std::endl;
    amrex::Print() << "p min: " << p.min(0) << std::endl;
    amrex::Print() << "u max: " << u.max(0) << std::endl;
    amrex::Print() << "u min: " << u.min(0) << std::endl;
    amrex::Print() << "v max: " << v.max(0) << std::endl;
    amrex::Print() << "v min: " << v.min(0) << std::endl;

    // **********************************
    // Write initial plot file
    // **********************************

    amrex::Real time = 0.0;

    // Interpolate the values to the cell center for writing output
    amrex::MultiFab output_values(cell_box_array, distribution_mapping, 4, 0);

    if (plot_interval > 0)
    {
        int time_step = 0;
        write_output(psi, p, u, v, geom, time, time_step, output_values);
    }

    // **********************************************
    // Intermediate Values used in time stepping loop
    // **********************************************

    // cu on the y faces (same locaions as u)
    amrex::MultiFab cu(u.boxArray(), u.DistributionMap(), 1, u.nGrow());

    // cv on the x faces (same locations as v)
    amrex::MultiFab cv(v.boxArray(), v.DistributionMap(), 1, v.nGrow());

    // z on the cell centers (same locations as psi)
    amrex::MultiFab z(psi.boxArray(), psi.DistributionMap(), 1, psi.nGrow());

    // h on the nodal points (same locations as p)
    amrex::MultiFab h(p.boxArray(), p.DistributionMap(), 1, p.nGrow());

    amrex::MultiFab u_old(u.boxArray(), u.DistributionMap(), u.nComp(), u.nGrow());
    amrex::MultiFab v_old(v.boxArray(), v.DistributionMap(), v.nComp(), v.nGrow());
    amrex::MultiFab p_old(p.boxArray(), p.DistributionMap(), p.nComp(), p.nGrow());

    amrex::MultiFab u_new(u.boxArray(), u.DistributionMap(), u.nComp(), u.nGrow());
    amrex::MultiFab v_new(v.boxArray(), v.DistributionMap(), v.nComp(), v.nGrow());
    amrex::MultiFab p_new(p.boxArray(), p.DistributionMap(), p.nComp(), p.nGrow());

    amrex::MultiFab::Copy(u_old, u, 0, 0, u.nComp(), u.nGrow());
    amrex::MultiFab::Copy(v_old, v, 0, 0, v.nComp(), v.nGrow());
    amrex::MultiFab::Copy(p_old, p, 0, 0, p.nComp(), p.nGrow());

    // Constants used in time stepping loop
    double fsdx = 4.0/dx;
    double fsdy = 4.0/dy;
    double tdt = dt;

    for (int time_step = 0; time_step < n_time_steps; ++time_step)
    {
        // fill ghost cells and periodic ghost cells 
        u.FillBoundary(geom.periodicity());
        v.FillBoundary(geom.periodicity());
        p.FillBoundary(geom.periodicity());

        for (amrex::MFIter mfi(p); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();

            // Read only arrays
            const amrex::Array4<amrex::Real const>& p_array = p.const_array(mfi);
            const amrex::Array4<amrex::Real const>& u_array = u.const_array(mfi);
            const amrex::Array4<amrex::Real const>& v_array = v.const_array(mfi);

            // Write arrays
            const amrex::Array4<amrex::Real>& cu_array = cu.array(mfi);
            const amrex::Array4<amrex::Real>& cv_array = cv.array(mfi);
            const amrex::Array4<amrex::Real>& h_array =   h.array(mfi);
            const amrex::Array4<amrex::Real>& z_array =   z.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                cu_array(i,j,k) = 0.5*(p_array(i,j,k) + p_array(i+1,j,k))*u_array(i,j,k);
                cv_array(i,j,k) = 0.5*(p_array(i,j,k) + p_array(i,j+1,k))*v_array(i,j,k);
                z_array(i,j,k) = (fsdx*(v_array(i+1,j,k)-v_array(i,j,k)) + fsdy*(u_array(i,j+1,k)-u_array(i,j,k)))/(p_array(i,j,k)+p_array(i+1,j,k)+p_array(i,j+1,k)+p_array(i+1,j+1,k));
                h_array(i,j,k) = p_array(i,j,k) + 0.25*(u_array(i-1,j,k)*u_array(i-1,j,k) + u_array(i,j,k)*u_array(i,j,k) + v_array(i,j-1,k)*v_array(i,j-1,k) + v_array(i,j,k)*v_array(i,j,k));
            });
        }

        cu.FillBoundary(geom.periodicity());
        cv.FillBoundary(geom.periodicity());
        h.FillBoundary(geom.periodicity());
        z.FillBoundary(geom.periodicity());

        // defined here because tdt changes after first time step
        double tdts8 = tdt/8.0;
        double tdtsdx = tdt/dx;
        double tdtsdy = tdt/dy;

        for (amrex::MFIter mfi(p); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();

            // Read only arrays
            const amrex::Array4<amrex::Real const>& p_old_array = p_old.const_array(mfi);
            const amrex::Array4<amrex::Real const>& u_old_array = u_old.const_array(mfi);
            const amrex::Array4<amrex::Real const>& v_old_array = v_old.const_array(mfi);
            const amrex::Array4<amrex::Real const>& cu_array = cu.const_array(mfi);
            const amrex::Array4<amrex::Real const>& cv_array = cv.const_array(mfi);
            const amrex::Array4<amrex::Real const>& h_array =   h.const_array(mfi);
            const amrex::Array4<amrex::Real const>& z_array =   z.const_array(mfi);

            // Write arrays
            const amrex::Array4<amrex::Real>& p_new_array = p_new.array(mfi);
            const amrex::Array4<amrex::Real>& u_new_array = u_new.array(mfi);
            const amrex::Array4<amrex::Real>& v_new_array = v_new.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                u_new_array(i,j,k) = u_old_array(i,j,k) + tdts8 * (z_array(i,j-1,k)+z_array(i,j,k)) * (cv_array(i,j-1,k) + cv_array(i,j,k) + cv_array(i+1,j-1,k) + cv_array(i+1,j,k)) - tdtsdx * (h_array(i+1,j,k) - h_array(i,j,k));
                v_new_array(i,j,k) = v_old_array(i,j,k) - tdts8 * (z_array(i-1,j,k)+z_array(i,j,k)) * (cu_array(i-1,j,k) + cu_array(i-1,j+1,k) + cu_array(i,j,k) + cu_array(i,j+1,k)) - tdtsdy * (h_array(i,j+1,k) - h_array(i,j,k));
                p_new_array(i,j,k) = p_old_array(i,j,k) - tdtsdx * (cu_array(i,j,k) - cu_array(i-1,j,k)) - tdtsdy * (cv_array(i,j,k) - cv_array(i,j-1,k));
            });
        }

        u_new.FillBoundary(geom.periodicity());
        v_new.FillBoundary(geom.periodicity());
        p_new.FillBoundary(geom.periodicity());

        time = time + dt;

        if (time_step>0) {

            double alpha = 0.001;

            for (amrex::MFIter mfi(p); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();

                // Read only arrays
                const amrex::Array4<amrex::Real const>& p_array = p.const_array(mfi);
                const amrex::Array4<amrex::Real const>& u_array = u.const_array(mfi);
                const amrex::Array4<amrex::Real const>& v_array = v.const_array(mfi);
                const amrex::Array4<amrex::Real const>& p_new_array = p_new.const_array(mfi);
                const amrex::Array4<amrex::Real const>& u_new_array = u_new.const_array(mfi);
                const amrex::Array4<amrex::Real const>& v_new_array = v_new.const_array(mfi);

                // Write arrays
                const amrex::Array4<amrex::Real>& p_old_array = p_old.array(mfi);
                const amrex::Array4<amrex::Real>& u_old_array = u_old.array(mfi);
                const amrex::Array4<amrex::Real>& v_old_array = v_old.array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    amrex::Real u_old_temp = u_old_array(i,j,k);
                    amrex::Real v_old_temp = v_old_array(i,j,k);
                    amrex::Real p_old_temp = p_old_array(i,j,k);

                    u_old_array(i,j,k) = u_array(i,j,k) + alpha * (u_new_array(i,j,k) - 2.0*u_array(i,j,k) + u_old_temp);
                    v_old_array(i,j,k) = v_array(i,j,k) + alpha * (v_new_array(i,j,k) - 2.0*v_array(i,j,k) + v_old_temp);
                    p_old_array(i,j,k) = p_array(i,j,k) + alpha * (p_new_array(i,j,k) - 2.0*p_array(i,j,k) + p_old_temp);
                });
            }

            amrex::MultiFab::Copy(u, u_new, 0, 0, u.nComp(), u.nGrow());
            amrex::MultiFab::Copy(v, v_new, 0, 0, v.nComp(), v.nGrow());
            amrex::MultiFab::Copy(p, p_new, 0, 0, p.nComp(), p.nGrow());

        } else {
            tdt = tdt+tdt;

            amrex::MultiFab::Copy(u_old, u, 0, 0, u.nComp(), u.nGrow());
            amrex::MultiFab::Copy(v_old, v, 0, 0, v.nComp(), v.nGrow());
            amrex::MultiFab::Copy(p_old, p, 0, 0, p.nComp(), p.nGrow());

            amrex::MultiFab::Copy(u, u_new, 0, 0, u.nComp(), u.nGrow());
            amrex::MultiFab::Copy(v, v_new, 0, 0, v.nComp(), v.nGrow());
            amrex::MultiFab::Copy(p, p_new, 0, 0, p.nComp(), p.nGrow());

        }

        // Write a plotfile of the current data (plot_interval was defined in the inputs file)
        if (plot_interval > 0 && time_step%plot_interval == 0)
        {
            write_output(psi, p, u, v, geom, time, time_step, output_values);
        }

        //amrex::Print() << "Done with time step " << time_step << std::endl;
    }

    amrex::Print() << "Final: " << std::endl;
    amrex::Print() << "p max: " << p.max(0) << std::endl;
    amrex::Print() << "p min: " << p.min(0) << std::endl;
    amrex::Print() << "u max: " << u.max(0) << std::endl;
    amrex::Print() << "u min: " << u.min(0) << std::endl;
    amrex::Print() << "v max: " << v.max(0) << std::endl;
    amrex::Print() << "v min: " << v.min(0) << std::endl;

    }
    amrex::Finalize();
    return 0;
}