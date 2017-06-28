/*
 * Vox-FE (Version 2): Voxel-based Finite Element Modelling
 *
 * Interface developed as a ParaView plugin by
 * Richard Holbrey and Michael J. Fagan
 *
 * Dept of Mechanical Engineering
 * University of Hull
 * Cottingham Road
 * HU6 7RX
 *
 * Vox-FE (Version 1) was created by Andreas Bitternas, George Sisias, Jia Liu and others.
 */

#ifndef __voxfeStrainFilter_h
#define __voxfeStrainFilter_h


#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkDoubleArray.h>

#include "../itkReader/ReadVoxFEScript.h"
#include "voxfeStrainData.h"

/** \class voxfeStrainFilter
 *  \brief voxfeStrainFilter reads in tabled data from the solver. The data
 *  are used to overlay the points in each multi-block with vtk named arrays,
 *  so these can be used to colour the surface or block geometry.
 */
class voxfeStrainFilter : public vtkMultiBlockDataSetAlgorithm { 

  public:
    vtkTypeMacro(voxfeStrainFilter,vtkMultiBlockDataSetAlgorithm);
    static voxfeStrainFilter* New();
    
    void PrintSelf(ostream &os, vtkIndent indent);

    vtkSetStringMacro(StrainFile);
    vtkGetStringMacro(StrainFile);

    //vtkGetMacro(DisplacementsOnly, bool);
    //vtkSetMacro(DisplacementsOnly, bool);

    vtkGetMacro(BMUSolverOutput, bool);
    vtkSetMacro(BMUSolverOutput, bool);

    // Make sure the pipeline knows what type we expect as input
    int FillInputPortInformation( int port, vtkInformation* info );
    
    void SetInput0Connection(vtkAlgorithmOutput* algOutput);


  protected:
    voxfeStrainFilter();
    ~voxfeStrainFilter();
    void operator=(const voxfeStrainFilter&);  // Not implemented.
 
    /**  \brief Performs vtk array allocation. The input geometry is shallow-copied.
     */
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    /** \brief Read in the strain/displacement data.
     *  \param if BMUSolverOutput is true, actualNodeCount is used instead of nNodes,
     *    due to output from BMU solver (which is not aggregated).
     *  \return NULL if an error is encountered or a pointer to the table data.
     *  The variable strainTable is (re-)allocated and returned here to ensure deletion; for
     *  actual usage the mapping in mapGlobalToStrainData is preferred.
     */
    double* ReadStrainData( ifstream& fin, unsigned long dim[3], const unsigned long& actualNodeCount );

    void OutputIfCountError( const vtkIdType& expected, const vtkIdType& count, const int& blk );

    //try to force inline
    inline void GetGlobalNodeNumber( const double* pt, const double& voxelsize,
                                     const double& dimXY, const double& dimX,
                                     unsigned long& globalnode  );
    char* StrainFile;
    voxfeStrainData* strainData;
    
    double* strainTable;                                  ///< Stores strain/displacement data.
    map< unsigned long, double* > mapGlobalToStrainData;  ///< mapping from global node number to strain data.

    //bool DisplacementsOnly;  ///< Set to read just displacements (strains will be computed on import)
	bool BMUSolverOutput; // For use with old solver -- maybe removed in future...?

	rvsIntegerType* nodeMap; //To re-read node mapping from BMU solver

};


#endif

