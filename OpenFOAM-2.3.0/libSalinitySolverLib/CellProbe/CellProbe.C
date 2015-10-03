 #include "CellProbe.H"
 #include "volFields.H"
 #include "OStringStream.H"
 #include "volFields.H"
 
Foam::CellProbe::CellProbe(const Time& runTime, const fvMesh& mesh)
:
    runTime_(runTime),
    mesh_(mesh),
    probeName_("probe"),
    probe_(false),
    grabMeshMotion_(false)
{

    IOdictionary cellProbeDict
    (
        IOobject
        (
            "cellProbeDict",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    probe_ = cellProbeDict.lookupOrDefault<Switch> ("probe", false);
    if (!probe_) return;

    if (Pstream::parRun())
    {
	OStringStream sstream;
	sstream << "processor" << Pstream::myProcNo();
	processorName_ = sstream.str();	
    }
    else
    {
	processorName_ = "serialRun";
    }
    
    cellProbeDict.lookup("grabMeshMotion") >> grabMeshMotion_;
    
    cellProbeDict.lookup("outFrequency") >> outFrequency_;
    stepCounter_ = 0;

    cellProbeDict.lookup("cellsCoords") >> cellsCoords_;
    //Info << cellsCoords << endl;

    cellProbeDict.lookup("fieldsNames") >> fieldsNames_;
    //Info << fieldsNames << endl;

    cellProbeDict.lookup("fileNamePrefix") >> fileNamePrefix_;
    cellProbeDict.lookup("fileNameExt") >> fileNameExt_;
    fileName probesFileName = fileNamePrefix_ + "probes." + fileNameExt_;
    fileName coordsFileName = fileNamePrefix_ + "coords." + fileNameExt_;

    pathToProbes_ = runTime_.rootPath() + "/" + runTime_.caseName() + "/cellProbes";
    probesOutputFileName_ = pathToProbes_ + "/" + probesFileName;
    coordsOutputFileName_ = pathToProbes_ + "/" + coordsFileName;

    bool done = mkDir (pathToProbes_, 0750);

    if (!done)
    {
        Info << "CellProbe: Directory doesn\'t exist! Can\'t create the new one!" << endl;
	Foam::FatalError.exit();
    }

    Info << "CellProbe path:   " << pathToProbes_ << endl;
    Info << "CellProbe probes: " << probesOutputFileName_ << endl;
    Info << "CellProbe coords: " << coordsOutputFileName_ << endl;

    probesOutput_.reset
    (
	new OFstream
	(
	    probesOutputFileName_
	)
    );

    coordsOutput_.reset
    (
	new OFstream
	(
	    coordsOutputFileName_
	)
    );

    cellsLabels_.resize(cellsCoords_.size());
    coordsOutput_() << "# x y z" << endl;
    
    
    
    forAll (cellsCoords_, jCell)
    {
	coordsOutput_() << "lookup point ";
        coordsOutput_() << jCell << " " << cellsCoords_[jCell].x()
                                 << " " << cellsCoords_[jCell].y()
                                 << " " << cellsCoords_[jCell].z()
                                 << endl;

        cellsLabels_[jCell] = mesh_.findCell(cellsCoords_[jCell]);
	
	if ( cellsLabels_[jCell] > -1 )
	{
//	    coordsOutput_() << "got cell : " << cellsLabels_[jCell] << "mesh len: " << mesh_.C().size() << endl;
//	    coordsOutput_() << "lookup cell  ";
//	    coordsOutput_() << jCell << " " << mesh_.C()[cellsLabels_[jCell]].x()
//                            << " " << mesh_.C()[cellsLabels_[jCell]].y()
//                            << " " << mesh_.C()[cellsLabels_[jCell]].z()
//                            << endl;
	}
    }

    probesOutput_() << "Time ";
    
    forAll(cellsLabels_, jCell)
    {
	if (cellsLabels_[jCell] > -1)
	{
	    forAll(fieldsNames_, jField)
	    {
		bool isScalar     = false;
		bool isVector     = false;
		bool isTensor     = false;
		bool isSymmTensor = false;
		
		{
		    isScalar     = mesh_.thisDb().foundObject<volScalarField>(fieldsNames_[jField]);
		    isVector     = mesh_.thisDb().foundObject<volVectorField>(fieldsNames_[jField]);
		    isTensor     = mesh_.thisDb().foundObject<volTensorField>(fieldsNames_[jField]);
		    isSymmTensor = mesh_.thisDb().foundObject<volSymmTensorField>(fieldsNames_[jField]);
		}
		if (isScalar)
		{
        	    probesOutput_() << probeName_ << "_" << jCell << "_" << fieldsNames_[jField] << " ";
		}
		else if (isVector)
		{
		    for (label j=0; j<3; j++)
		    {
			probesOutput_() << probeName_ << "_" << jCell << "_" << fieldsNames_[jField] << 
			vector::componentNames[j] << " ";
		    }
		}
		else if (isTensor)
		{
		    for (label j=0; j<9; j++)
		    {
			probesOutput_() << probeName_ << "_" << jCell << "_" << fieldsNames_[jField] << 
			tensor::componentNames[j] << " ";
		    }
		}
		else if (isSymmTensor)
		{
		    for (label j=0; j<6; j++)
		    {
			probesOutput_() << probeName_ << "_" << jCell << "_" << fieldsNames_[jField] << 
			symmTensor::componentNames[j] << " ";
		    }
		}
		else
		{
		    Info << "CellProbe: Can\'t find field " << fieldsNames_[jField] << endl;
		    Foam::FatalError.exit();
		}
	    
	    }
	}

    }
    probesOutput_() << endl;
}


void Foam::CellProbe::probe()
{
    if (!probe_) return;
    
    stepCounter_++;
    if (stepCounter_ != outFrequency_) return;
    stepCounter_ = 0;

    probesOutput_() << runTime_.timeName() << " ";
    
    Info << "Starting probe...." << endl;
    
    if (!grabMeshMotion_)
    {
	forAll (cellsCoords_, jCell)
	{
	    cellsLabels_[jCell] = mesh_.findCell(cellsCoords_[jCell]);
	}
    }
    
    forAll(cellsLabels_, jCell)
    {
	if (cellsLabels_[jCell] > -1)
	{
	    forAll(fieldsNames_, jField)
	    {
		bool isScalar     = false;
		bool isVector     = false;
		bool isTensor     = false;
		bool isSymmTensor = false;
		
		{
		    isScalar     = mesh_.thisDb().foundObject<volScalarField>(fieldsNames_[jField]);
		    isVector     = mesh_.thisDb().foundObject<volVectorField>(fieldsNames_[jField]);
		    isTensor     = mesh_.thisDb().foundObject<volTensorField>(fieldsNames_[jField]);
		    isSymmTensor = mesh_.thisDb().foundObject<volSymmTensorField>(fieldsNames_[jField]);
		}
		if (isScalar)
		{
		    volScalarField field = mesh_.thisDb().lookupObject<volScalarField>(fieldsNames_[jField]);
        	    probesOutput_() << field[cellsLabels_[jCell]] << " ";
		}
		else if (isVector)
		{
		    volVectorField field = mesh_.thisDb().lookupObject<volVectorField>(fieldsNames_[jField]);
		    vector val = field[cellsLabels_[jCell]];
		    for (label j=0; j<3; j++)
		    {
			probesOutput_() << val[j] << " ";
		    }
		}
		else if (isTensor)
		{
		    volTensorField field = mesh_.thisDb().lookupObject<volTensorField>(fieldsNames_[jField]);
		    tensor val = field[cellsLabels_[jCell]];
		    for (label j=0; j<9; j++)
		    {
			probesOutput_() << val[j] << " ";
		    }
		}
		else if (isSymmTensor)
		{
		    volSymmTensorField field = mesh_.thisDb().lookupObject<volSymmTensorField>(fieldsNames_[jField]);
		    symmTensor val = field[cellsLabels_[jCell]];
		    for (label j=0; j<6; j++)
		    {
			probesOutput_() << val[j] << " ";
		    }
		}
		else
		{
		    Info << "CellProbe: Can\'t find field " << fieldsNames_[jField] << endl;
		    Foam::FatalError.exit();
		}
    	    }
	}
    }
    
    probesOutput_() << endl;

}
//END_OF_FILE


