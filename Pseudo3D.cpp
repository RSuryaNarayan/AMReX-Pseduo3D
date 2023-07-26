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
  std::cerr << "This script takes a 2-D AMReX plotfile and writes a 3-D plotfile which is one cell thick along the third direction. Usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> \n\tOptions:\n\tis_per=<L M N> ";
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
	// Set defaults input values
	// ---------------------------------------------------------------------
	std::string infile         = "";
	std::string outfile        = "";
	Vector<int> is_per(AMREX_SPACEDIM,1);	   
	int finestLevel            = 0;
	
	// ---------------------------------------------------------------------
	// ParmParse
	// ---------------------------------------------------------------------
	ParmParse pp;

 	pp.get("infile",infile);
    	pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
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
	Vector<IntVect> refRatios(finestLevel-1);
	Vector<int> levelSteps(finestLevel);
	IntVect nGrow;

    	// Read data on all the levels
	for (int lev=0; lev<finestLevel; ++lev) {
	const BoxArray ba = plotfile_2D.boxArray(lev);
	grids[lev] = ba;
	dmap[lev] = plotfile_2D.DistributionMap(lev);
	geoms[lev] = Geometry(plotfile_2D.probDomain(lev),&rb,coord,&(is_per[0]));
	if (lev < finestLevel-1) {
		refRatios[lev] = IntVect(plotfile_2D.refRatio(lev));
	}
	levelSteps[lev] = plotfile_2D.levelStep(lev);
	nGrow = plotfile_2D.nGrowVect(lev);
	state[lev].define(grids[lev], dmap[lev], nCompIn, nGrow);
	Print() << "Writing data for level: " << lev << std::endl;
	state[lev] = plotfile_2D.get(lev);
	}
	
	Real time = plotfile_2D.time();
	
	//----------------------------------------------------------------------
	// Assemble and write plotfile out
	//----------------------------------------------------------------------
	WriteMultiLevelPlotfile (outfile,
				finestLevel,
				GetVecOfConstPtrs(state),
				plotVarNames,
				geoms,
				time,
				levelSteps,
				refRatios);
	std::cout<<"\n Wrote 3D plotfile "<<outfile<<"\n";				
    }
  amrex::Finalize();
}
