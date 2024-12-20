/*
 * A simplified single file version of the HeatEquation_EX0_C exmaple.
 * This code is designed to be used with Demo_Tutorial.rst.
 *
 */

#include <numbers>
#include <cmath>

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>


int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {

    // **********************************
    // DECLARE SIMULATION PARAMETERS
    // **********************************

    // number of cells on each side of the domain
    int n_cell;

    // size of each box (or grid)
    int max_grid_size;

    // total steps in simulation
    int nsteps;

    // how often to write a plotfile
    int plot_int;

    // time step
    amrex::Real dt;

    // **********************************
    // READ PARAMETER VALUES FROM INPUT DATA
    // **********************************
    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        amrex::ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // time step
        pp.get("dt",dt);
    }

    // **********************************
    // DEFINE SIMULATION SETUP AND GEOMETRY
    // **********************************

    // make BoxArray and Geometry
    // ba will contain a list of boxes that cover the domain
    // geom contains information such as the physical domain size,
    // number of points in the domain, and periodicity

    // define lower and upper indices
    amrex::IntVect dom_lo(0,0);
    amrex::IntVect dom_hi(n_cell-1, n_cell-1);

    // Make a single box that is the entire domain
    //amrex::Box domain(dom_lo, dom_hi);
    
    amrex::Box cell_centered_box(dom_lo, dom_hi);
    amrex::Box x_face_centered_box(dom_lo, dom_hi+1, amrex::IndexType({1,0}));
    amrex::Box y_face_centered_box(dom_lo, dom_hi+1, amrex::IndexType({0,1}));
    amrex::Box node_centered_box = amrex::surroundingNodes(cell_centered_box);

    amrex::Print() << "Cell centered box " << cell_centered_box << std::endl;
    amrex::Print() << "X face centered box " << x_face_centered_box << std::endl;
    amrex::Print() << "y face centered box " << y_face_centered_box << std::endl;
    amrex::Print() << "node centered box " << node_centered_box << std::endl;

    // This defines the physical box, [0,1] in each direction.
//    amrex::RealBox real_box({ 0., 0.},
//                     { 1., 1.});
    amrex::RealBox real_box({ 0., 0.},
                     { 2*std::numbers::pi, 2*std::numbers::pi});

    // This, a value of 0, says we are using Cartesian coordinates
    int coordinate_system = 0;

    // This sets the boundary conditions in each direction to periodic
    amrex::Array<int,AMREX_SPACEDIM> is_periodic {1,1};

    // This defines a Geometry object
    //amrex::Geometry geom(cell_centered_box, real_box, coordinate_system, is_periodic);
    amrex::Geometry geom(cell_centered_box, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::Print() << "geom " << geom << std::endl;

    
    amrex::Geometry geom2;
    {
      amrex::RealBox real_box({ 0., 0.},
                       { 2*std::numbers::pi, 2*std::numbers::pi});

      // This, a value of 0, says we are using Cartesian coordinates
      int coordinate_system = 0;

      // This sets the boundary conditions in each direction to periodic
      amrex::Array<int,AMREX_SPACEDIM> is_periodic {1,1};

      // This defines a Geometry object
      //amrex::Geometry geom(cell_centered_box, real_box, coordinate_system, is_periodic);
      //amrex::Geometry geom2(cell_centered_box, real_box, amrex::coordSys::Cartesian, is_periodic);
      geom2.define(cell_centered_box, real_box, amrex::CoordSys::cartesian, is_periodic);

      amrex::Print() << "geom2 " << geom2 << std::endl;

    }

    // Initialize the boxarray "ba" from the single box "domain"
    amrex::BoxArray ba;
    //ba.define(domain);
    ba.define(cell_centered_box);

    amrex::Print() << "ba before max size " << ba << std::endl;
    

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    amrex::Print() << "ba after max size " << ba << std::endl;


    // extract dx from the geometry object
    amrex::GpuArray<amrex::Real,2> dx = geom.CellSizeArray();

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    // Ncomp = number of components for each array
    int Ncomp = 1;

    // How Boxes are distrubuted among MPI processes
    amrex::DistributionMapping dm(ba);
    amrex::Print() << "dm " << dm << std::endl;

    // we allocate two phi multifabs; one will store the old state, the other the new.
    amrex::MultiFab phi_old(ba, dm, Ncomp, Nghost);
    amrex::MultiFab phi_new(ba, dm, Ncomp, Nghost);

    // time = starting time in the simulation
    amrex::Real time = 0.0;

    // **********************************
    // INITIALIZE DATA LOOP
    // **********************************

    // coeffiecent for initialization of stream function
    amrex::Real a = 1000000;

    // loop over boxes
    for (amrex::MFIter mfi(phi_old); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& phiOld = phi_old.array(mfi);

        // set phi = 1 + e^(-(r-0.5)^2)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            amrex::Real x_cell_center = (i+0.5) * dx[0];
            amrex::Real y_cell_center = (j+0.5) * dx[1];
            amrex::Real x_node = (i) * dx[0];
            amrex::Real y_node = (j) * dx[1];

            phiOld(i,j,k) = a*std::sin(x_cell_center)*std::sin(y_cell_center);
        });
    }

    // **********************************
    // WRITE INITIAL PLOT FILE
    // **********************************

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,5);
        WriteSingleLevelPlotfile(pltfile, phi_old, {"phi"}, geom, time, 0);
    }


//    // **********************************
//    // MAIN TIME EVOLUTION LOOP
//    // **********************************
//
//    for (int step = 1; step <= nsteps; ++step)
//    {
//        // fill periodic ghost cells
//        phi_old.FillBoundary(geom.periodicity());
//
//        // new_phi = old_phi + dt * Laplacian(old_phi)
//        // loop over boxes
//        for ( amrex::MFIter mfi(phi_old); mfi.isValid(); ++mfi )
//        {
//            const amrex::Box& bx = mfi.validbox();
//
//            const amrex::Array4<amrex::Real>& phiOld = phi_old.array(mfi);
//            const amrex::Array4<amrex::Real>& phiNew = phi_new.array(mfi);
//
//            // advance the data by dt
//            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
//            {
//
//                // **********************************
//                // EVOLVE VALUES FOR EACH CELL
//                // **********************************
//
//                //phiNew(i,j,k) = phiOld(i,j,k) + dt *
//                //    ( (phiOld(i+1,j,k) - 2.*phiOld(i,j,k) + phiOld(i-1,j,k)) / (dx[0]*dx[0])
//                //     +(phiOld(i,j+1,k) - 2.*phiOld(i,j,k) + phiOld(i,j-1,k)) / (dx[1]*dx[1])
//                //     +(phiOld(i,j,k+1) - 2.*phiOld(i,j,k) + phiOld(i,j,k-1)) / (dx[2]*dx[2])
//                //        );
//
//                phiNew(i,j,k) = phiOld(i,j,k) + dt *
//                    ( (phiOld(i+1,j,k) - 2.*phiOld(i,j,k) + phiOld(i-1,j,k)) / (dx[0]*dx[0])
//                     +(phiOld(i,j+1,k) - 2.*phiOld(i,j,k) + phiOld(i,j-1,k)) / (dx[1]*dx[1])
//                        );
//            });
//        }
//
//        // **********************************
//        // INCREMENT
//        // **********************************
//
//        // update time
//        time = time + dt;
//
//        // copy new solution into old solution
//        amrex::MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);
//
//        // Tell the I/O Processor to write out which step we're doing
//        amrex::Print() << "Advanced step " << step << "\n";
//
//
//        // **********************************
//        // WRITE PLOTFILE AT GIVEN INTERVAL
//        // **********************************
//
//        // Write a plotfile of the current data (plot_int was defined in the inputs file)
//        if (plot_int > 0 && step%plot_int == 0)
//        {
//            const std::string& pltfile = amrex::Concatenate("plt",step,5);
//            WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, step);
//        }
//    }


    }
    amrex::Finalize();
    return 0;
}


