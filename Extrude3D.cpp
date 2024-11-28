#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BLFort.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "This script takes a 2-D AMReX plotfile and writes a 3-D plotfile with specified thickness along the third direction. Usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> n_cells=<number_of_cells> \n\tOptions:\n\tis_per=<L M N> ";
  exit(1);
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    if (argc < 2) {
      print_usage(argc,argv);
    }

    // ---------------------------------------------------------------------
    // Set default input values
    // ---------------------------------------------------------------------
    std::string infile         = "";
    std::string outfile        = "";
    Vector<int> is_per(AMREX_SPACEDIM,1);	   
    int finestLevel            = 0;
    int n_cells                = 1;  // Default to 1 cell thick if n_cells not specified
    
    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;
    
    pp.get("infile",infile);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    pp.query("n_cells", n_cells);
    outfile = infile + "_3D";

    //----------------------------------------------------------------------
    // Read in the plotfile using PlotFileData
    //----------------------------------------------------------------------
    PlotFileData plotfile_2D(infile);
    
    //----------------------------------------------------------------------
    // Get plotfile features to write to 3D
    //----------------------------------------------------------------------
    finestLevel = plotfile_2D.finestLevel();
    std::cout<<"Finest Level = "<<finestLevel<<"\n";
    finestLevel = finestLevel+1;

    int nCompIn = plotfile_2D.nComp();
    std::cout<<"\nNumber of components = "<<nCompIn<<"\n";

    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i) {
      destFillComps[i] = i;
    }

    Vector<std::string> plotVarNames = plotfile_2D.varNames();
    std::cout<<"\n Variables available in the plotfile:\n";
    
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      std::cout<<i+1<<".)"<<plotVarNames[i]<<"\n";
    }

    Vector<MultiFab> state(finestLevel);
    Vector<Geometry> geoms(finestLevel);
    Vector<BoxArray> grids(finestLevel);
    Vector<DistributionMapping> dmap(finestLevel);
    RealBox rb(&(plotfile_2D.probLo()[0]), 
               &(plotfile_2D.probHi()[0]));
    int coord = plotfile_2D.coordSys();
    IntVect nGrow;

    // ---------------------------------------------------------------------
    // Modify the domain in the z-direction to be n_cells thick
    // ---------------------------------------------------------------------
    const BoxArray ba2D = plotfile_2D.boxArray(0);
    BoxArray ba3D = ba2D;  // Start with 2D BoxArray, extend to 3D
    
    ba3D.enclosedCells();  // Ensure cells are defined at cell centers
    ba3D.maxSize({ba2D[0].length(0), ba2D[0].length(1), n_cells});  // Set thickness in z-direction

    // Define 3D geometry with n_cells in z-direction
    Box domain_3D(ba3D.minimalBox());
    domain_3D.setBig(2, n_cells - 1); // Set extent in z-direction
    
    Geometry geom3D(domain_3D, &rb, coord, is_per.data());
    
    // Set up data structures
    grids[0] = ba3D;
    dmap[0] = plotfile_2D.DistributionMap(0);
    geoms[0] = geom3D;
    nGrow = plotfile_2D.nGrowVect(0);
    
    state[0].define(grids[0], dmap[0], nCompIn, nGrow);

    // Copy 2D data into each z-slice of the 3D MultiFab
    Print() << "Writing data for level: 0" << std::endl;
    const MultiFab& mf_2D = plotfile_2D.get(0);
    for (MFIter mfi(state[0]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        FArrayBox& fab = state[0][mfi];
        const FArrayBox& fab_2D = mf_2D[mfi];
        for (int z = 0; z < n_cells; ++z) {
            Box slice = bx;
            slice.setSmall(2, z);
            slice.setBig(2, z);
            fab.copy(fab_2D, bx, 0, slice, 0, nCompIn);
        }
    }

    Real time = plotfile_2D.time();
    Vector<int> isteps(finestLevel, 0);
    Vector<IntVect> refRatios(finestLevel-1, {AMREX_D_DECL(2, 2, 2)});
    
    //----------------------------------------------------------------------
    // Assemble and write plotfile out
    //----------------------------------------------------------------------
    WriteMultiLevelPlotfile(outfile,
                            finestLevel,
                            GetVecOfConstPtrs(state),
                            plotVarNames,
                            geoms,
                            time,
                            isteps,
                            refRatios);
    std::cout<<"\n Wrote 3D plotfile "<<outfile<<"\n";				
  }
  amrex::Finalize();
}
