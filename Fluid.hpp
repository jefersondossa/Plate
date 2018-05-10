//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------FLUID-------------------------------------
//------------------------------------------------------------------------------

#ifndef FLUID_H
#define FLUID_H

#include "Element.hpp"
#include "Boundary.hpp"

// PETSc libraries
#include <metis.h>
#include <petscksp.h> 

/// Mounts the incompressible flow problem

template<int DIM>
class Fluid{
public:
    /// Defines the class Element locally
    typedef Element<DIM> Elements;

    /// Defines the class Node locally
    typedef typename Elements::Nodes  Node;

    /// Defines the class Boundary locally
    typedef Boundary<DIM> Boundaries;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       nodes_;

    /// Defines the vector of fluid elements
    std::vector<Elements *>   elements_;
 
    /// Defines the vector of fluid boundaries mesh nodes
    std::vector<Boundaries *> boundary_;

private:
    //FLUID VARIABLES
    std::string inputFile; //Fluid input file
    int numElem;           //Number of elements in fluid mesh 
    int numNodes;          //Number of nodes in velocity/quadratic mesh
    int numBoundaries;     //Number of fluid boundaries
    int numBoundElems;     //Number of elements in fluid boundaries
    double pressInf;       //Undisturbed pressure 
    double rhoInf;         //Density
    double tempInf;        //Temperature
    double viscInf;        //Viscosity
    double ktermInf;       //Thermal condutivity
    double velocityInf[3]; //Undisturbed velocity
    double fieldForces[3]; //Field forces (constant)
    idx_t* part_elem;      //Fluid Domain Decomposition - Elements
    idx_t* part_nodes;     //Fluid Domain Decomposition - Nodes
    int numTimeSteps;      //Number of Time Steps
    int printFreq;         //Printing frequence of output files
    double dTime;          //Time Step
    double integScheme;    //Time Integration Scheme
    int rank;
    int numFSIInterfaces;
    int iAux;
    double dragCoefficient;
    double liftCoefficient;

public:
    /// Reads the input file and perform the preprocessing operations
    /// @param std::string input file @param std::string mirror file
    void dataReading(std::string inputFile, std::string mirror);

    /// Performs the domain decomposition for parallel processing
    void domainDecompositionMETIS(); 

    /// Export the domain decomposition 
    /// @return pair with the elements and nodes domain decompositions
    std::pair<idx_t*,idx_t*> getDomainDecomposition(){
        return std::make_pair(part_elem,part_nodes);};

    /// Gets the number of time steps
    /// @return number of time steps
    int getNumberOfTimeSteps(){return numTimeSteps;};

    /// Gets the time step size
    /// @return time step size
    double getTimeStep(){return dTime;};

    /// Gets the number of fluid-structure interfaces
    /// @return number of fluid boundaries which composes the 
    /// fluid structure interface
    int getNumberofFSIInterfaces(){return numFSIInterfaces;};

    /// Mounts and solve the steady incompressible flow steady problem    
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int problem type: 1 - Stokes problem; 2 - Navier-Stokes problem.
    int solveSteadyProblem(int iterNumber, double tolerance, int problem_type);

    /// Mounts and solve the transient incompressible flow problem    
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int problem type: 1 - Stokes problem; 2 - Navier-Stokes problem.
    int solveTransientProblem(int iterNumber,double tolerance,int problem_type);

    /// Mounts and solve the steady Laplace problem
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    int solveSteadyLaplaceProblem(int iterNumber, double tolerance);

    /// Mounts and solve the transient incompressible flow problem for 
    /// for fluid structure interaction problems
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int problem type: 1 - Stokes problem; 2 - Navier-Stokes problem.
    int solveFSIFluid(int iterNumber,
                      double tolerance,
                      int problem_type);

    /// Print the results for Paraview post-processing
    /// @param int time step
    void printVelocity(int step);

    /// Gets the fluid model nodes and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model nodes information
    std::vector<Node *> getNodesVelocity(){return nodes_;}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<Elements *> getElements(){return elements_;}
};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//---------------------SUBDIVIDES THE FINITE ELEMENT DOMAIN---------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::domainDecompositionMETIS() {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numElem;
    idx_t numNd = numNodes;
    idx_t dd = 2;
    idx_t ssize = size;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[(4*dd-2)*numEl];
    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];


    for (idx_t i = 0; i < numEl+1; i++){
        elem_start[i]=(4*dd-2)*i;
    };
    for (idx_t jel = 0; jel < numEl; jel++){
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();        
        
        for (idx_t i=0; i<(4*dd-2); i++){
        elem_connec[(4*dd-2)*jel+i] = connec(i);
        };
    };

    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec, \
                              NULL, NULL, &one, &ssize, NULL, NULL,    \
                              &objval, part_elem, part_nodes);

    mirrorData << std::endl \
               << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
    for(int i = 0; i < numElem; i++){
        mirrorData << "process = " << part_elem[i] \
                   << ", element = " << i << std::endl;
    };

    mirrorData << std::endl \
               << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
    for(int i = 0; i < numNodes; i++){
        mirrorData << "process = " << part_nodes[i] \
                   << ", node = " << i << std::endl;
    };
    
    return;

};

//------------------------------------------------------------------------------
//----------------------------PRINT VELOCITY RESULTS----------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::printVelocity(int step) {

    //    std::cout << "Printing Velocity Results" << std::endl;
   

    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();
    std::string s = "saidaVel"+result+".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << numNodes
             << "\"  NumberOfCells=\"" << numElem
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodes; i++){
        typename Node::VecLocD x;
        x=nodes_[i]->getCoordinates();
        output_v << x(0) << " " << x(1) << " " << 0.0 << std::endl;        
    };
    output_v << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i=0; i<numElem; i++){
        typename Elements::Connectivity connec;
        connec=elements_[i]->getConnectivity();
        output_v << connec(0) << " " << connec(1) << " " << connec(2) << " " \
                 << connec(3) << " " << connec(4) << " " << connec(5) << \
            std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    int aux = 0;
    for (int i=0; i<numElem; i++){
        output_v << aux + 6 << std::endl;
        aux += 6;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;
    
    for (int i=0; i<numElem; i++){
        output_v << 22 << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << std::endl;
    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Velocity\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodes; i++){
        output_v << nodes_[i] -> getVelocity(0) << " "                \
                 << nodes_[i] -> getVelocity(1) << " " << 0. << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    // output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //          << "Name=\"Divergent\" format=\"ascii\">" << std::endl;

    // for (int i=0; i<numNodes; i++){
    //     output_v << nodes_[i] -> getVelocityDivergent() << std::endl;
    // };
    // output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Pressure\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodes; i++){
        output_v << nodes_[i] -> getPressure() << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Mesh mov\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodes; i++){       

        typename Node::VecLocD x, xp;
 
        x=nodes_[i]->getCoordinates();
        xp=nodes_[i]->getInitialCoordinates();

        output_v << x(0)-xp(0) << " " << x(1)-xp(1) << " " << 0. << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;


    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Gradient\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodes; i++){        
        output_v << nodes_[i] -> getGradientComponent(0) << " "
                 << nodes_[i] -> getGradientComponent(1) << " " 
                 << 0. << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Process\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElem; i++){
        output_v << part_elem[i] << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Jacobian\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElem; i++){
        output_v << elements_[i] -> getJacobian() << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;
    
    output_v << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << std::endl;
    
    // output_v << "  <FieldData>" << std::endl;

    // output_v << "      <DataArray type=\"Float64\" Name=\"Time\" NumberOfTuples=\"1\" "
    //          << " format=\"ascii\">" << std::endl;
    // output_v << step << std::endl;
    // output_v << "      </DataArray> " << std::endl;

    // output_v << "      <DataArray type=\"Float64\" Name=\"LiftCoefficient\" NumberOfTuples=\"1\" "
    //          << " format=\"ascii\">" << std::endl;
    // output_v << liftCoefficient << std::endl;
    // output_v << "      </DataArray> " << std::endl;

    // // output_v << "      <DataSet type=\"Float64\" Name=\"Drag Coefficient\" NumberOfTuples=\"1\" "
    // //          << " format=\"ascii\">" << std::endl;
    // // output_v << dragCoefficient << std::endl;
    // // output_v << "      </DataSet> " << std::endl;

    // output_v << "  </FieldData>" << std::endl
    output_v << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;

};

//------------------------------------------------------------------------------
//----------------------------READS FLUID INPUT FILE----------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::dataReading(std::string inputFile, std::string mirror) {

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      

    if (rank == 0)std::cout << "Reading fluid data.." << std::endl;

    std::string line;
    
    //Defines input and output files    
    std::ifstream inputData(inputFile.c_str());
    std::ofstream mirrorData(mirror.c_str());

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read number of nodes, elements, time steps and printing frequence
    inputData >> numElem >> numNodes >> numTimeSteps >> printFreq;
    mirrorData << "Number of Elements     = " << numElem << std::endl;
    mirrorData << "Number of Nodes        = " << numNodes << std::endl;
    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read undisturbed velocity and pressure components
    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;

    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read undisturbed density, temperature, viscosity and thermal condutivity
    inputData >> rhoInf >> tempInf >> viscInf >> ktermInf;

    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
    mirrorData << "Undisturbed Temperature= " << tempInf << std::endl;
    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;
    mirrorData << "Undisturbed Term. Cond.= " << ktermInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read time step lenght
    inputData >> dTime >> integScheme;

    mirrorData << "Time Step              = " << dTime << std::endl;
    mirrorData << "Time Integration Scheme= " << integScheme << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //Read field forces
    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];

    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl \
               << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    int dimension=2;

    int index = 0;

    //Read nodal coordinates
    nodes_.reserve(numNodes);

    for (int inode=0; inode<numNodes; inode++){
        typename Node::VecLocD x;
        for (int i=0; i<dimension; i++){
                inputData >> x(i);
        };
        Node *node = new Node(x, index++);
        nodes_.push_back(node);
    };

    mirrorData << "Nodal Coordinates" << std::endl;
    for (int i = 0 ; i<numNodes; i++){
        typename Node::VecLocD x;
        x=nodes_[i]->getCoordinates();       
        for (int j=0; j<dimension; j++){
            mirrorData << x(j) << " ";
        };
        mirrorData << std::endl;
        nodes_[i] -> setVelocity(velocityInf);
        nodes_[i] -> setPreviousVelocity(velocityInf);
    };

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read element connectivity
    elements_.reserve(numElem);
    index = 0;
    for (int jel=0; jel<numElem; jel++){
        typename Elements::Connectivity connect;

        for (int i=0; i<4*dimension-2; i++){
            inputData >> connect(i);
            connect(i)=connect(i)-1;
        };
        Elements *el = new Elements(index++,connect,nodes_);
        elements_.push_back(el);
    };    
        
    mirrorData << std::endl << "Element Connectivity" << std::endl;        
    for (int jel=0; jel<numElem; jel++){
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();       
        for (int i=0; i<4*dimension-2; i++){
            mirrorData << connec(i) << " ";
        };
        mirrorData << std::endl;
    };

    //Sets Viscosity, density, time step, time integration 
    //scheme and field forces values
    for (int i=0; i<numElem; i++){
        elements_[i] -> setViscosity(viscInf);
        elements_[i] -> setDensity(rhoInf);
        elements_[i] -> setTimeStep(dTime);
        elements_[i] -> setTimeIntegrationScheme(integScheme);
        elements_[i] -> setFieldForce(fieldForces);
    };
        
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);
    
    //Read number of boundary groups
    inputData >> numBoundaries >> numBoundElems;
    mirrorData << "Number of Boundary groups = " << numBoundaries << std::endl;
    mirrorData << "Number of Boundary Elements = " << numBoundElems <<std::endl;

    boundary_.reserve(numBoundElems);

    numFSIInterfaces=0;


    //Read boundary information: connectivity and prescribed value
    for (int ibound=0; ibound<numBoundaries; ibound++){

        getline(inputData,line);getline(inputData,line);getline(inputData,line);
        getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
        int constrain[3];
        double value[3];
        int numGroup;
        Boundaries::BoundConnect connectB;
        index = 0;
        
        inputData >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] \
                  >> value[0] >> value[1] >> value[2];

        if ((constrain[0] == 3) || 
            (constrain[1] == 3) || 
            (constrain[2] == 3)) numFSIInterfaces++;
        
        mirrorData << "Constrains = " << constrain[0] << " " << constrain[1] \
                   << " " << constrain[2] << std::endl;
        mirrorData << "Values = " << value[0] << " " << value[1]  \
                   << " " << value[2] << std::endl;
        
        getline(inputData,line);getline(inputData,line);getline(inputData,line);
        getline(inputData,line);getline(inputData,line);

        mirrorData << "Connectivity = " << std::endl;
        
        for (int i = 0; i < numGroup; i++){
            inputData >> connectB(0) >> connectB(1) >> connectB(2);
            
            connectB(0) -= 1;
            connectB(1) -= 1;
            connectB(2) -= 1;

            Boundaries *bound = new Boundaries(connectB, index++,  \
                                               constrain, value, ibound);
            boundary_.push_back(bound);

            mirrorData << connectB(0) << " " << connectB(1) \
                       << " " << connectB(2) << std::endl;
        };
    };

    // for (int ino = 0; ino < numNodes; ino++){
    //     nodes_[ino] -> setConstrains(0, 0, 0);
    //     nodes_[ino] -> setConstrains(1, 0, 0);
    // };
        

    //Sets boundary constrains
    for (int ibound = 0; ibound < numBoundElems; ibound++){
        
        Boundaries::BoundConnect connectB;
        connectB = boundary_[ibound] -> getBoundaryConnectivity();
        int no1 = connectB(0);
        int no2 = connectB(1);
        int no3 = connectB(2);

        nodes_[no1] -> setConstrainsLaplace(0,1,0);
        nodes_[no2] -> setConstrainsLaplace(0,1,0);
        nodes_[no3] -> setConstrainsLaplace(0,1,0);

        nodes_[no1] -> setConstrainsLaplace(1,1,0);
        nodes_[no2] -> setConstrainsLaplace(1,1,0);
        nodes_[no3] -> setConstrainsLaplace(1,1,0);
        
        if (boundary_[ibound] -> getConstrain(0) == 1){

            //Desfazer primeira parte do if para voltar a cond. cont. constante
            // if (boundary_[ibound] -> getConstrainValue(0) <= 1.){
                
            //     typename Node::VecLocD x;
            //     x = nodes_[no1]->getCoordinates();                 
                
            //     nodes_[no1] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
                
            //     x = nodes_[no2]->getCoordinates();                 
                
            //     nodes_[no2] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
                
            //     x = nodes_[no3]->getCoordinates();                 
                
            //     nodes_[no3] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
            //     //ate aqui
            // } else {
            nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
            nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
            nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
             // };
        };
        if (boundary_[ibound] -> getConstrain(1) == 1){
            nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
            nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
            nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
        };               
    };
    
    //Print nodal constrains
    for (int i=0; i<numNodes; i++){
        mirrorData<< "Constrains " << i \
                  << " " << nodes_[i] -> getConstrains(0) \
                  << " " << nodes_[i] -> getConstrainValue(0) \
                  << " " << nodes_[i] -> getConstrains(1) \
                  << " " << nodes_[i] -> getConstrainValue(1) << std::endl;
    }; 

    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<numBoundElems; i++){
        if ((boundary_[i] -> getConstrain(0) > 0) ||
            (boundary_[i] -> getConstrain(1) > 0)) {

            Boundaries::BoundConnect connectB;
            connectB = boundary_[i] -> getBoundaryConnectivity();

            for (int j=0; j<numElem; j++){
                typename Elements::Connectivity connect;
                connect = elements_[j] -> getConnectivity();
                
                int flag = 0;
                int side[3];
                for (int k=0; k<6; k++){
                    if ((connectB(0) == connect(k)) || 
                        (connectB(1) == connect(k)) ||
                        (connectB(2) == connect(k))){
                        side[flag] = k;
                        flag++;
                    };
                };
                
                if (flag == 3){
                    boundary_[i] -> setElement(j);
                    //Sets element index and side
                    if ((side[0]==4) || (side[1]==4) || (side[2]==4)){
                        boundary_[i] -> setElementSide(0);
                        elements_[boundary_[i]->getElement()] -> 
                            setElemSideInBoundary(0);
                    };
                    if ((side[0]==5) || (side[1]==5) || (side[2]==5)){
                        boundary_[i] -> setElementSide(1);
                        elements_[boundary_[i]->getElement()] -> 
                            setElemSideInBoundary(1);
                    };
                    if ((side[0]==3) || (side[1]==3) || (side[2]==3)){
                        boundary_[i] -> setElementSide(2);
                        elements_[boundary_[i]->getElement()] -> 
                            setElemSideInBoundary(2);
                    };
                };
            };
        };
    };


    domainDecompositionMETIS();

    iAux = 0;


};


//------------------------------------------------------------------------------
//--------------------------SOLVE STEADY FLUID PROBLEM--------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveSteadyProblem(int iterNumber, double tolerance,\
                                 int problem_type) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //    MatNullSpace      nullsp;

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
        
    //Updates SUPG Parameter
    for (int i = 0; i < numElem; i++){
        elements_[i] -> getParameterSUPG();
    };
        

    for (int inewton = 0; inewton < iterNumber; inewton++){

        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, \
                            2*numNodes+numNodes, 2*numNodes+numNodes, \
                            50,NULL,150,NULL,&A); CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodes);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);

        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;
        
        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                
                //Compute Element matrix
                if (problem_type == 1) elements_[jel] -> getSteadyStokes();
                
                if (problem_type == 2) {
                    if (inewton == 0) {
                        elements_[jel] -> getSteadyStokes();
                    } else {
                        elements_[jel] -> getSteadyNavierStokes();
                    };
                };
                
                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec;
                
                //Gets element connectivity, jacobian and rhs 
                connec = elements_[jel] -> getConnectivity();
                Ajac = elements_[jel] -> getJacNRMatrix();
                Rhs = elements_[jel] -> getRhsVector();
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                &Ajac(2*i  ,2*j  ),
                                                ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,2*j  ),
                                                ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,2*j+1),
                                                ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,2*j+1),
                                                ADD_VALUES);
                        };
                        
                        //Matrix Q and Qt
                        if (fabs(Ajac(2*i  ,12+j)) >= 1.e-8){
                                int dof_i = 2*connec(i);
                                int dof_j = 2*numNodes + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j, \
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i, \
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };
                        if (fabs(Ajac(2*i+1,12+j)) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*numNodes + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,12+j),
                                                ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(12+j,2*i+1),
                                                ADD_VALUES);
                        };
                        if (fabs(Ajac(12+i,12+j)) >= 1.e-8){
                            int dof_i = 2*numNodes + connec(i);
                            int dof_j = 2*numNodes + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(12+i,12+j),
                                                ADD_VALUES);
                        };
                    };
                    
                    //Rhs vector
                    if (fabs(Rhs(2*i  )) >= 1.e-8){
                        int dof_i = 2*connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
                                            ADD_VALUES);
                    };
                    
                    if (fabs(Rhs(2*i+1)) >= 1.e-8){
                        int dof_i = 2*connec(i)+1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
                                            ADD_VALUES);
                    };
                };
                for (int i=0; i<6; i++){
                    if (fabs(Rhs(12+i)) >= 1.e-8){
                        int dof_i = 2*numNodes + connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),
                                                ADD_VALUES);
                    };
                };
            };
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        // ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
        ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
                                100);CHKERRQ(ierr);
        
        ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
      
        ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
        
        ierr = PCSetType(pc, PCASM);CHKERRQ(ierr);
        
        // ierr = KSPSetPCSide(ksp, PC_RIGHT);
        //ierr = KSPSetType(ksp,KSPBCGS); CHKERRQ(ierr);

        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        

        // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE,0, NULL, &nullsp);
        // ierr = MatSetNullSpace(A, nullsp);
        // ierr = MatNullSpaceDestroy(&nullsp);

// #if defined(PETSC_HAVE_MUMPS)
//             ierr = KSPSetType(ksp,KSPPREONLY);
//             ierr = KSPGetPC(ksp,&pc);
//             ierr = PCSetType(pc, PCLU);
// #endif
            
//             ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//             ierr = KSPSetUp(ksp);
            



        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);

        std::cout << "GMRES Iterations = " << iterations << std::endl;
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
                
        //Updates nodal values
        double u_ [2];
        double p_;
        Ione = 1;

        for (int i = 0; i < numNodes; ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[1] = val;
            nodes_[i] -> incrementVelocity(0,u_[0]);
            nodes_[i] -> incrementVelocity(1,u_[1]);
        };
        for (int i = 0; i<numNodes; i++){
            Ii = 2*numNodes+i;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            nodes_[i] -> incrementPressure(p_);
        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        if(rank == 0){
            std::cout << "Du Norm = " << val << std::endl;
            std::cout<<"Iteration = " << inewton << std::endl;
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        if(val <= tolerance){
            break;            
        };
        
        
    };

    if (rank == 0) {
        //Computing velocity divergent

        for (int ino = 0; ino < numNodes; ino++){
            nodes_[ino] -> setVelocityDivergent(0.);
        };                
       
        for (int jel = 0; jel < numElem; jel++){   
            elements_[jel] -> computeVelocityDivergent();
        };

        printVelocity(1);
    };

    return 0;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveTransientProblem(int iterNumber, double tolerance,\
                                 int problem_type) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //    MatNullSpace      nullsp;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::string dl = "dragLift.dat";
    std::ofstream dragLift(dl.c_str());
    
    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
        

    for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){

        for (int i = 0; i < numElem; i++){
            elements_[i] -> getParameterSUPG();
        };

        
        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}
        
        for (int i = 0; i < numNodes; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodes_[i] -> getVelocity(0);
            u[1] = nodes_[i] -> getVelocity(1);
            
            uprev[0] = nodes_[i] -> getPreviousVelocity(0);
            uprev[1] = nodes_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodes_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodes_[i] -> setPreviousVelocity(u);
            
            //Computes predictor
            // u[0] += dTime * accel[0];
            // u[1] += dTime * accel[1];
            
            // nodes_[i] -> setVelocity(u);
            
            // //Updates acceleration
            // nodes_[i] -> setPreviousAcceleration(accel);
            
        };

        double duNorm=100.;
        
        for (int inewton = 0; inewton < iterNumber; inewton++){
            boost::posix_time::ptime t1 =                             
                               boost::posix_time::microsec_clock::local_time();
            
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2*numNodes+numNodes, 2*numNodes+numNodes,
                                100,NULL,300,NULL,&A); 
            CHKERRQ(ierr);
            
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodes);
            CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
            
            //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;

            for (int jel = 0; jel < numElem; jel++){   
                
                if (part_elem[jel] == rank) {
                    //Compute Element matrix
                    if (problem_type == 1)
                        elements_[jel] -> getTransientStokes();
                    
                    if (problem_type == 2){
                        if (iTimeStep < 2){
                            elements_[jel] -> getTransientNavierStokes();
                        } else {
                            elements_[jel] -> getTransientNavierStokes();
                        };
                    };

                    typename Elements::LocalMatrix Ajac;
                    typename Elements::LocalVector Rhs;
                    typename Elements::Connectivity connec;
                    
                    //Gets element connectivity, jacobian and rhs 
                    connec = elements_[jel] -> getConnectivity();
                    Ajac = elements_[jel] -> getJacNRMatrix();
                    Rhs = elements_[jel] -> getRhsVector();
                    
                    //Disperse local contributions into the global matrix
                    //Matrix K and C
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A, 1, &dof_i,1, &dof_j,
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };
                        
                            //Matrix Q and Qt
                            if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,12+j),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &Ajac(12+j,2*i+1),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                int dof_i = 2 * numNodes + connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(12+i,12+j),
                                                    ADD_VALUES);
                            };
                        };
                        
                        //Rhs vector
                        if (fabs(Rhs(2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(2*i  ),
                                                ADD_VALUES);
                        };
                        
                        if (fabs(Rhs(2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i)+1;
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(2*i+1),
                                                ADD_VALUES);
                        };
                        if (fabs(Rhs(12+i)) >= 1.e-15){
                            int dof_i = 2 * numNodes + connec(i);
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(12+i),
                                                ADD_VALUES);
                        };
                    };
                };
            }; //Elements
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            
            // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            



        //     ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
        //                             500);CHKERRQ(ierr);
            
        //     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            
        //     ierr = KSPGetPC(ksp,&pc);
            
        //     ierr = PCSetType(pc,PCASM);
            
        //     ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);

        //     ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);
            
        //        //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
            

        // //   //   ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nullsp);
        // // // ierr = MatSetNullSpace(A, nullsp);
        // // // ierr = MatNullSpaceDestroy(&nullsp);

 

#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp,KSPPREONLY);
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCLU);
#endif          
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);



            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);

            ierr = KSPGetTotalIterations(ksp, &iterations);            

            //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
            
            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            double p_;
            duNorm = 0.;
            double dpNorm = 0.;
            Ione = 1;
            
            for (int i = 0; i < numNodes; ++i){
                //if (nodes_[i] -> getConstrains(0) == 0){
                    Ii = 2*i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    duNorm += val*val;
                    nodes_[i] -> incrementVelocity(0,val);
                    //}; 
                
                    //if (nodes_[i] -> getConstrains(1) == 0){
                    Ii = 2*i+1;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    duNorm += val*val;
                    nodes_[i] -> incrementVelocity(1,val);
                    //};
            };
            
            for (int i = 0; i<numNodes; i++){
                Ii = 2*numNodes+i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                dpNorm += val*val;
                nodes_[i] -> incrementPressure(p_);
            };
            
            //Computes the solution vector norm
            //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
   
            boost::posix_time::ptime t2 =                               \
                               boost::posix_time::microsec_clock::local_time();
         
            if(rank == 0){
                boost::posix_time::time_duration diff = t2 - t1;

                std::cout << "Iteration = " << inewton 
                          << " (" << iterations << ")"  
                          << "   Du Norm = " << std::scientific << sqrt(duNorm) 
                          << " " << sqrt(dpNorm)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds()/1000. << std::endl;
            };
                      
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);

            if (sqrt(duNorm) <= tolerance) {
                break;
            };

            //Updates SUPG Parameter
            // for (int i = 0; i < numElem; i++){
            //     elements_[i] -> getParameterSUPG();
            // };
            
        };//Newton-Raphson

        dragCoefficient = 0.;
        liftCoefficient = 0.;

        for (int jel = 0; jel < numBoundElems; jel++){   

            ublas::bounded_vector<double,2> load;
            load.clear();

            if (boundary_[jel] -> getBoundaryGroup() == 0){               
                int iel = boundary_[jel] -> getElement();
                load = elements_[iel] -> getDragAndLiftForces();               
            };
            
            dragCoefficient += load(0) / (0.5 * rhoInf * velocityInf[0]);
            liftCoefficient += load(1) / (0.5 * rhoInf * velocityInf[0]);

        };

        if (rank == 0) {
        //Computing velocity divergent
            
            dragLift << iTimeStep * dTime << " " << dragCoefficient 
                     << " " << liftCoefficient << std::endl;
            
            // for (int ino = 0; ino < numNodes; ino++){
            //     nodes_[ino] -> setVelocityDivergent(0.);
            // };                
            
            // for (int jel = 0; jel < numElem; jel++){   
            //     elements_[jel] -> computeVelocityDivergent();
            // };

            printVelocity(iTimeStep);
        };
        
    };
    
    return 0;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE STEADY LAPLACE PROBLEM-------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveSteadyLaplaceProblem(int iterNumber, double tolerance) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);        

    for (int inewton = 0; inewton < iterNumber; inewton++){

        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2*numNodes, 2*numNodes,
                            30,NULL,30,NULL,&A); CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;
        
        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                
                //Compute Element matrix
                elements_[jel] -> getSteadyLaplace();
                                
                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec;
                
                //Gets element connectivity, jacobian and rhs 
                connec = elements_[jel] -> getConnectivity();
                Ajac = elements_[jel] -> getJacNRMatrix();
                Rhs = elements_[jel] -> getRhsVector();
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,2*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,2*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,2*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,2*j+1),ADD_VALUES);
                        };
                    };
                                        
                    //Rhs vector
                    if (fabs(Rhs(2*i  )) >= 1.e-8){
                        int dof_i = 2*connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);
                    };
                    
                    if (fabs(Rhs(2*i+1)) >= 1.e-8){
                        int dof_i = 2*connec(i)+1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),ADD_VALUES);
                    };
                };
            };
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
#endif
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);
        
        
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);

        //std::cout << "GMRES Iterations = " << iterations << std::endl;
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
                
        //Updates nodal values
        double u_ [2];
        Ione = 1;

        for (int i = 0; i < numNodes; ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[1] = val;
            nodes_[i] -> incrementCoordinate(0,u_[0]);
            nodes_[i] -> incrementCoordinate(1,u_[1]);
        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        if(rank == 0){
            std::cout << "MESH MOVING - ERROR = " << val 
                      << std::scientific <<  std::endl;
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        if(val <= tolerance){
            break;            
        };             
    };
    
    // for (int i=0; i<numElem; i++){
    //     elements_[i] -> computeNodalGradient();            
    // };

    if (rank == 0) {
        //Computing velocity divergent
        //      printVelocity(1);
    };

    return 0;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveFSIFluid(int iterNumber, double tolerance, int problem_type){

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //    MatNullSpace      nullsp;

    int rank;

    iAux++;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
        
    for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
        double duNorm=100.;
        
        //Updates SUPG Parameter
        for (int i = 0; i < numElem; i++){
            elements_[i] -> getParameterSUPG();
        };
             
        for (int inewton = 0; inewton < iterNumber; inewton++){
            boost::posix_time::ptime t1 =                             
                               boost::posix_time::microsec_clock::local_time();
            
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2*numNodes+numNodes, 2*numNodes+numNodes,
                                100,NULL,300,NULL,&A); 
            CHKERRQ(ierr);
            
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodes);
            CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
            
            //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;

            for (int jel = 0; jel < numElem; jel++){   
                
                if (part_elem[jel] == rank) {
                    //Compute Element matrix
                    if (problem_type == 1)
                        elements_[jel] -> getTransientStokes();
                    
                    if (problem_type == 2){
                        if (iTimeStep < 2){
                            elements_[jel] -> getTransientNavierStokes();
                        } else {
                            elements_[jel] -> getTransientNavierStokes();
                        };
                    };

                    typename Elements::LocalMatrix Ajac;
                    typename Elements::LocalVector Rhs;
                    typename Elements::Connectivity connec;
                    
                    //Gets element connectivity, jacobian and rhs 
                    connec = elements_[jel] -> getConnectivity();
                    Ajac = elements_[jel] -> getJacNRMatrix();
                    Rhs = elements_[jel] -> getRhsVector();
                    
                    //Disperse local contributions into the global matrix
                    //Matrix K and C
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A, 1, &dof_i,1, &dof_j,
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };
                        
                            //Matrix Q and Qt
                            if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,12+j),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &Ajac(12+j,2*i+1),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                int dof_i = 2 * numNodes + connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(12+i,12+j),
                                                    ADD_VALUES);
                            };
                        };
                        
                        //Rhs vector
                        if (fabs(Rhs(2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(2*i  ),
                                                ADD_VALUES);
                        };
                        
                        if (fabs(Rhs(2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i)+1;
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(2*i+1),
                                                ADD_VALUES);
                        };
                        if (fabs(Rhs(12+i)) >= 1.e-15){
                            int dof_i = 2 * numNodes + connec(i);
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(12+i),
                                                ADD_VALUES);
                        };
                    };
                };
            }; //Elements
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            
            // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            



        //     ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
        //                             500);CHKERRQ(ierr);
            
        //     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            
        //     ierr = KSPGetPC(ksp,&pc);
            
        //     ierr = PCSetType(pc,PCJACOBI);
            
        //     //ierr = KSPSetType(ksp,KSPBCGS); CHKERRQ(ierr);

        //     // ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
            
        //        //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
            

        //   //   ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nullsp);
        // // ierr = MatSetNullSpace(A, nullsp);
        // // ierr = MatNullSpaceDestroy(&nullsp);

 

#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp,KSPPREONLY);
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCLU);
#endif          
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);



            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);

            ierr = KSPGetTotalIterations(ksp, &iterations);            

            //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
            
            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            double p_;
            duNorm = 0.;
            double dpNorm = 0.;
            Ione = 1;
            
            for (int i = 0; i < numNodes; ++i){
                //if (nodes_[i] -> getConstrains(0) == 0){
                    Ii = 2*i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    duNorm += val*val;
                    nodes_[i] -> incrementVelocity(0,val);
                    //}; 
                
                    //if (nodes_[i] -> getConstrains(1) == 0){
                    Ii = 2*i+1;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    duNorm += val*val;
                    nodes_[i] -> incrementVelocity(1,val);
                    //};
            };
            
            for (int i = 0; i<numNodes; i++){
                Ii = 2*numNodes+i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                dpNorm += val*val;
                nodes_[i] -> incrementPressure(p_);
            };
            
            //Computes the solution vector norm
            //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
   
            boost::posix_time::ptime t2 =                               \
                               boost::posix_time::microsec_clock::local_time();
         
            if(rank == 0){
                boost::posix_time::time_duration diff = t2 - t1;

                std::cout << "Iteration = " << inewton 
                          << " (" << iterations << ")"  
                          << "   Du Norm = " << std::scientific << sqrt(duNorm) 
                          << " " << sqrt(dpNorm)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds()/1000. << std::endl;
            };
                      
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);

            if (sqrt(duNorm) <= tolerance) {
                break;
            };


   
        };//Newton-Raphson
            
    };
    
    return 0;
};

#endif
