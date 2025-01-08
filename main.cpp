/*
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

    // TODO: a hack to make sure that we use one process
    n_cell = 64;
    max_grid_size = 1000000000;

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
    //amrex::Box x_face_centered_box(dom_lo, dom_hi+1, amrex::IndexType({1,0}));
//    amrex::Box y_face_centered_box(dom_lo, dom_hi+1, amrex::IndexType({0,1}));
//    amrex::Box node_centered_box = amrex::surroundingNodes(cell_centered_box);

    //amrex::Box x_face_centered_box = amrex::convert(cell_centered_box, amrex::IndexType({1,0}));
    //amrex::Box y_face_centered_box = amrex::convert(cell_centered_box, amrex::IndexType({0,1}));

    //amrex::Print() << "Cell centered box " << cell_centered_box << std::endl;
    //amrex::Print() << "X face centered box " << x_face_centered_box << std::endl;
    //amrex::Print() << "y face centered box " << y_face_centered_box << std::endl;
//    amrex::Print() << "node centered box " << node_centered_box << std::endl;
    
    amrex::Geometry geom;
    {
      amrex::RealBox real_box({ 0., 0.},
                       { 2*std::numbers::pi, 2*std::numbers::pi});

      // This, a value of 0, says we are using Cartesian coordinates
      //int coordinate_system = 0;

      // This sets the boundary conditions in each direction to periodic
      amrex::Array<int,AMREX_SPACEDIM> is_periodic {1,1};

      // This defines a Geometry object
      //amrex::Geometry geom(cell_centered_box, real_box, coordinate_system, is_periodic);
      //amrex::Geometry geom2(cell_centered_box, real_box, amrex::coordSys::Cartesian, is_periodic);
      geom.define(cell_centered_box, real_box, amrex::CoordSys::cartesian, is_periodic);

      amrex::Print() << "geom " << geom << std::endl;
    }

    // Initialize the boxarray "ba" from the single box "domain"
    amrex::BoxArray ba;
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

    amrex::MultiFab psi_old(ba, dm, Ncomp, Nghost);
    amrex::MultiFab psi_new(ba, dm, Ncomp, Nghost);


    //// Set up face fields
    //    // build the flux multifabs
    //Array<MultiFab, AMREX_SPACEDIM> flux;
    //for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    //{
    //    // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
    //    BoxArray edge_ba = ba;
    //    edge_ba.surroundingNodes(dir);
    //    flux[dir].define(edge_ba, dm, 1, 0);
    //}

    amrex::BoxArray x_face_box_array = amrex::convert(ba, {1,0});
    amrex::Print() << "x_face_box_array " << x_face_box_array << std::endl;
    amrex::MultiFab v(x_face_box_array, dm, Ncomp, 0);

    amrex::BoxArray y_face_box_array = amrex::convert(ba, {0,1});
    amrex::Print() << "y_face_box_array " << y_face_box_array << std::endl;
    amrex::MultiFab u(y_face_box_array, dm, Ncomp, 0);

    amrex::BoxArray surrounding_nodes_box_array= ba;
    surrounding_nodes_box_array.surroundingNodes();
    amrex::Print() << "surrounding_nodes_box_array " << surrounding_nodes_box_array << std::endl;
    amrex::MultiFab p(surrounding_nodes_box_array, dm, Ncomp, 0);

    //BoxArray x_face_ba = ba;
    //x_face_ba.surroundingNodes(0);
    //amrex::MultiFab phi_new(ba, dm, Ncomp, Nghost);

    // time = starting time in the simulation
    amrex::Real time = 0.0;

    // **********************************
    // INITIALIZE DATA LOOP
    // **********************************

    // coeffiecent for initialization of stream function
    amrex::Real a = 1000000;

    // loop over cell centers
    for (amrex::MFIter mfi(psi_old); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& phiOld = psi_old.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            amrex::Real x_cell_center = (i+0.5) * dx[0];
            amrex::Real y_cell_center = (j+0.5) * dx[1];
            //amrex::Real x_node = (i) * dx[0];
            //amrex::Real y_node = (j) * dx[1];

            phiOld(i,j,k) = a*std::sin(x_cell_center)*std::sin(y_cell_center);
        });
    }

    psi_old.FillBoundary(geom.periodicity());

    amrex::Print() << "psi_old min: " << psi_old.min(0) << std::endl;
    amrex::Print() << "psi_old max: " << psi_old.max(0) << std::endl;

    // coeffiecent for initialization of p
    int N = n_cell; // Change to read into input file later... choose this name to correspond with the name from the python script
    double mesh_dx = 100000;
    double el = N*mesh_dx;
    amrex::Real pcf = (std::numbers::pi * std::numbers::pi * a * a)/(el * el);

    for (amrex::MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& p_array = p.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            amrex::Real x_node = i * dx[0];
            amrex::Real y_node = j * dx[1];

            p_array(i,j,k) = pcf * (std::cos(2*x_node) + std::cos(2*y_node)) + 5000;
        });
    }

    //p.FillBoundary();
    //p.FillBoundary(geom.periodicity());

    amrex::Print() << "p min: " << p.min(0) << std::endl;
    amrex::Print() << "p max: " << p.max(0) << std::endl;

    double mesh_dy = 100000;
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& u_array = u.array(mfi);
        const amrex::Array4<amrex::Real>& phi_old_array = psi_old.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            u_array(i,j,k) = -(phi_old_array(i,j,k)-phi_old_array(i,j-1,k))/mesh_dy;
        });
    }

    amrex::Print() << "u min: " << u.min(0) << std::endl;
    amrex::Print() << "u max: " << u.max(0) << std::endl;

    for (amrex::MFIter mfi(v); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& v_array = v.array(mfi);
        const amrex::Array4<amrex::Real>& phi_old_array = psi_old.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            v_array(i,j,k) = (phi_old_array(i,j,k)-phi_old_array(i-1,j,k))/mesh_dx;
        });
    }

    // for writing values to output file so everything is at the cell centers
    //amrex::MultiFab output_values(ba, dm, 1, Nghost);
    //amrex::MultiFab output_values(ba, dm, 1, 0);
    //for (amrex::MFIter mfi(output_values); mfi.isValid(); ++mfi)
    //{
    //    const amrex::Box& bx = mfi.validbox();

    //    const amrex::Array4<amrex::Real const>& phi_old_array = psi_old.const_array(mfi);
    //    const amrex::Array4<amrex::Real const>& p_array = p.const_array(mfi);

    //    const amrex::Array4<amrex::Real>& output_values_array = output_values.array(mfi);

    //    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    //    {
    //        //output_values_array(i,j,k) = phi_old_array(i,j,k);
    //        output_values_array(i,j,k) = (p_array(i,j,k) + p_array(i+1,j,k) + p_array(i,j+1,k) + p_array(i+1,j+1,k))/4.0;
    //        //output_values_array(i,j,k) = 100;
    //    });
    //}

    amrex::MultiFab output_values(ba, dm, 4, 0);
    for (amrex::MFIter mfi(output_values); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real const>& phi_old_array = psi_old.const_array(mfi);
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

    // **********************************
    // WRITE INITIAL PLOT FILE
    // **********************************

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,5);
        //WriteSingleLevelPlotfile(pltfile, psi_old, {"phi"}, geom, time, 0);
        //WriteSingleLevelPlotfile(pltfile, p, {"p"}, geom, time, 0);
        WriteSingleLevelPlotfile(pltfile, output_values, {"psi_old", "p", "u", "v"}, geom, time, 0);
    }


//    // **********************************
//    // MAIN TIME EVOLUTION LOOP
//    // **********************************
//
//    for (int step = 1; step <= nsteps; ++step)
//    {
//        // fill periodic ghost cells
//        psi_old.FillBoundary(geom.periodicity());
//
//        // new_phi = old_phi + dt * Laplacian(old_phi)
//        // loop over boxes
//        for ( amrex::MFIter mfi(psi_old); mfi.isValid(); ++mfi )
//        {
//            const amrex::Box& bx = mfi.validbox();
//
//            const amrex::Array4<amrex::Real>& phiOld = psi_old.array(mfi);
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
//        amrex::MultiFab::Copy(psi_old, phi_new, 0, 0, 1, 0);
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


