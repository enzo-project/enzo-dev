/***********************************************************************
/
/  PROBLEM TYPE CLASS
/
/  written by: Matthew Turk, Oliver Hahn
/  date:       July, 2010
/
/  PURPOSE:
/
************************************************************************/
#ifdef NEW_PROBLEM_TYPES
#ifndef __PROBLEM_TYPES__
#define __PROBLEM_TYPES__

#include <string>

class EnzoProblemType
{
/* These will be overridden in the public namespace by implementations */
public:
    virtual int InitializeSimulation(FILE *fptr, FILE *Outfptr,
            HierarchyEntry &TopGrid, TopGridData &MetaData)
    {}
    virtual int InitializeGrid(grid *thisgrid,
            HierarchyEntry &TopGrid, TopGridData &MetaData)
    {}
protected:
    //.. constructor
    EnzoProblemType( )
    {}
    virtual ~EnzoProblemType()
    {}

    /*
       float **GetBaryonFields(grid *thisgrid) { return thisgrid->BaryonField; }
       int GetGridRank(grid *thisgrid) { return thisgrid->GridRank; }
       int OwnGrid(grid *thisgrid) {
       return (thisgrid->ProcessorNumber == MyProcessorNumber); }
       int *GetGridDimension(grid *thisgrid) { return thisgrid->GridDimension; }
       int GridIdentifyPhysicalQuantities(grid *thisgrid, int &DensNum, int &GENum,
       int &Vel1Num, int &Vel2Num, int &Vel3Num, int &TENum);

       int GridIdentifyPhysicalQuantities(grid *thisgrid, int &DensNum, int &GENum,
       int &Vel1Num, int &Vel2Num, int &Vel3Num, int &TENum,
       int &B1Num, int &B2Num, int &B3Num);

       int GridIdentifyPhysicalQuantities(grid *thisgrid, int &DensNum, int &GENum,
       int &Vel1Num, int &Vel2Num, int &Vel3Num, int &TENum,
       int &B1Num, int &B2Num, int &B3Num, int &PhiNum)
       { thisgrid->
     */

};

/*!
 * @brief implements abstract factory design pattern for plug-ins
 */
struct EnzoProblemType_creator
{
    //! create an instance of a plug-in
    virtual EnzoProblemType * create( ) const = 0;
    
    //! destroy an instance of a plug-in
    virtual ~EnzoProblemType_creator() { }
};

//! maps the name of a plug-in to a pointer of the factory pattern 
std::map<std::string, EnzoProblemType_creator *>& get_problem_types();

/*!
 * @brief concrete factory pattern for plug-ins
 */
template< class DerivedProblemType >
struct EnzoProblemType_creator_concrete : public EnzoProblemType_creator
{
    //! register the plug-in by its name
    EnzoProblemType_creator_concrete( const std::string& problem_type_name )
    {
        get_problem_types()[ problem_type_name ] = this;
    }
    
    //! create an instance of the plug-in
    EnzoProblemType * create( ) const
    {
        return new DerivedProblemType( );
    }
};

//! failsafe version to select the plug-in
EnzoProblemType *select_problem_type( std::string problem_type_name );

#endif
#endif
