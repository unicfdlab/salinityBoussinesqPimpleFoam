#ifndef CellProbe_H
#define CellProbe_H

#include "Time.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "autoPtr.H"

namespace Foam
{

class CellProbe
{

private:

    //-
    const Time& runTime_;
    
    //-
    const fvMesh& mesh_;

    //-
    const word probeName_;
    
    //-
    bool probe_;
    
    //-
    bool grabMeshMotion_;
    
    //-
    word fileNamePrefix_;
    
    //-
    word fileNameExt_;
    
    //-
    word processorName_;
    
    //-
    List<vector> cellsCoords_;
    
    //-
    List<label>  cellsLabels_;
    
    //-
    List<word>   fieldsNames_;
    
    //-
    label outFrequency_;
    
    //-
    label stepCounter_;
    
    //-
    fileName pathToProbes_;
    
    //-
    fileName probesOutputFileName_;
    
    //-
    fileName coordsOutputFileName_;
    
    //-
    autoPtr<OFstream> probesOutput_;
    
    //-
    autoPtr<OFstream> coordsOutput_;
    
    //-
    CellProbe();
    
    //-
    CellProbe (const CellProbe& );
    
public:
    
    //-
    CellProbe(const Time& runTime, const fvMesh& mesh);
    
    //-
    void probe();

};


}

#endif

//END_OF_FILE

