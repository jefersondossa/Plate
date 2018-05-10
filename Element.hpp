//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef ELEMENT_H
#define ELEMENT_H

#include "Node.hpp"
#include "BoundaryIntegrationQuadrature.hpp"
#include "IntegrationQuadrature.hpp"
#include "IntegrationQuadrature11.hpp"
#include "PartitionedQuadrature.hpp"

/// Defines the fluid element object and all the element information

template<int DIM>
class Element{
 
public:
    /// Defines the class Node locally
    typedef Node<DIM>                                           Nodes;

    /// Defines a type to store the element mesh connectivity
    typedef ublas::bounded_vector<int, 4*DIM-2>                 Connectivity;

    /// Defines a blas-type vector with dimension = DIM
    typedef ublas::bounded_vector<double, DIM>                  DimVector;

    /// Defines a blas-type matrix with dimension = DIM x DIM
    typedef ublas::bounded_matrix<double, DIM, DIM>             DimMatrix;

    /// Defines the vector which contains the element nodal coordinates
    typedef ublas::bounded_matrix<double, 4*DIM-2, DIM>         LocalNodes;

    /// Defines the local vector type with dimension 15 for DIM=2 
    /// and 34 for DIM=3
    typedef ublas::bounded_vector<double, 22*DIM-26>            LocalVector;

    /// Defines the local matrix type with dimension 15x15
    /// for DIM=2 and 34x34 for DIM=3
    typedef ublas::bounded_matrix<double, 22*DIM-26, 22*DIM-26> LocalMatrix;

    ///Defines the partitioned integration quadrature rule class locally
    //typedef PartQuadrature<DIM>                                 SpecialQuad;
    typedef IntegQuadratureSpecial<DIM>                         SpecialQuad;

    /// Defines the normal integration quadrature rule class locally
    typedef IntegQuadrature<DIM>                                NormalQuad;

    /// Defines the boundary integration quadrature rule class locally
    typedef BoundaryIntegQuadrature<DIM>                        BoundaryQuad;
    
    /// Define the type VecLocD from class Node locally
    typedef typename Nodes::VecLocD                             VecLoc;

private:
    QuadShapeFunction<DIM> shapeQuad; //Quadratic shape function
    BoundShapeFunction<DIM>shapeBound;//Boundary shape function
    SpecialQuad            sQuad;     //Integration quadrature
    NormalQuad             nQuad;     //Integration quadrature
    std::vector<Nodes *>   nodes_;    //Velocity nodes
    Connectivity  connect_;           //Velocity mesh connectivity 
    int           index_;             //Element index
    LocalNodes    localNodes_;        //Nodal coordinates - velocity
    LocalNodes    localNodesBoundary_;//Nodal coordinates - velocity
    double        djac_;              //Jacobian determinant
    double        weight_;            //Integration point weight
    DimMatrix     ainv_;              //Inverse of Jacobian transform
    double        visc_;              //Fluid dynamic viscosity
    double        dens_;              //Fluid Density
    double        dTime_;             //Time Step
    double        timeScheme_;        //Time Integration scheme
    LocalMatrix   diffMatrix;         //Diffusion / Viscosity matrix
    LocalMatrix   laplMatrix;         //Laplacian problem matrix
    LocalMatrix   lagrMultMatrix;     //Lagrange Multiplier matrix
    LocalMatrix   jacobianNRMatrix;   //Newton's method jacobian
    double        u_, v_, w_, p_;     //Interpolated velocity and pressure
    double        uPrev_, vPrev_, wPrev_, pPrev_;
    double        ax_, ay_, az_;      //Interpolated acceleration
    double        umesh_, vmesh_, wmesh_; //Interpolated mesh velocity
    double        du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, 
                  dw_dx, dw_dy, dw_dz;//Interpolated fluid spatial derivatives
    double        dp_dx, dp_dy, dp_dz;//Interpolated pressure spatial derivative
    double        dumesh_dx, dumesh_dy, dumesh_dz,//Interpolated mesh 
                  dvmesh_dx, dvmesh_dy, dvmesh_dz,//velocity derivatives
                  dwmesh_dx, dwmesh_dy, dwmesh_dz;
    double        lagMx_,lagMy_;
    double        tSUPG_;             //SUPG parameter
    double        tPSPG_;             //PSPG parameter
    double        tLSIC_;             //LSIC parameter
    double        tARLQ_;             //Arlequin Stabilization parameter
    double        tSUGN1_;
    double        tSUGN2_;
    double        tSUGN3_;
    double        hRGN_;
    DimVector     fieldForce_;        //Element field force
    LocalVector   externalForces;     //Field forces integrated over
                                      //the element domain
    LocalVector   rhsVector;          //RHS vector of Newton's method
    LocalVector   rhsVectorLM;        //Lag Mult RHS vector of Newton's method
    DimVector     xK, XK;
    double        dCk[DIM+1], dck[DIM+1];
    LocalNodes    gradient_;
    int           sideBoundary_;
    double        meshMovingParameter;
    static const double k1;
    static const double k2;

    typename SpecialQuad::PointWeight   intPointWeightFunction;
    typename SpecialQuad::PointWeight   intPointDistGlueZone;
    typename SpecialQuad::PointLogical  intPointGlueZone;
    bool          glueZone;
    bool          compressibility;   //True if the element is compressible
    bool          model; //true for local and false for global
    bool          FSIInterface;
    
    std::vector<ublas::bounded_vector<double,DIM> > di;

    //Second derivatives of velocity shape functions
    typename QuadShapeFunction<DIM>::ValueDDeriv ddphi_dx; 
    //First derivatives of velocity shape functions
    typename QuadShapeFunction<DIM>::ValueDeriv  dphi_dx, dphiL_dx;
    //Values of velocity shape functins
    typename QuadShapeFunction<DIM>::Values      phi_;     
    //Values of velocity shape functins
    typename BoundShapeFunction<DIM>::Values     phib_;     
    //Values of velocity shape functins
    typename BoundShapeFunction<DIM>::ValueDeriv dphib_;     

    std::vector<VecLoc>                 intPointCorrespXsi;
    std::vector<VecLoc>                 intPointCoordinates;
    std::vector<int>                    intPointCorrespElem;

public:
    /// fluid element constructor
    /// @param int element index @param Connectivity element connectivity
    /// @param vector<Nodes> 
    Element(int index, Connectivity& connect, std::vector<Nodes *> nodes){
        index_ = index;
        connect_ = connect;
        nodes_ = nodes;

        localNodes_.clear();  
        djac_ = 0.;             weight_ = 0.;             ainv_.clear();
        visc_ = 0.;             dens_ = 0.;               dTime_ = 0.; 
        timeScheme_ = 0.;       laplMatrix.clear();
        jacobianNRMatrix.clear();lagrMultMatrix.clear();
        u_ = 0.;          v_ = 0.;          w_ = 0;          p_ = 0.; 
        umesh_ = 0.;      vmesh_ = 0.;      wmesh_ = 0.; 
        du_dx = 0.;       du_dy = 0.;       du_dz = 0.; 
        dv_dx = 0.;       dv_dy = 0.;       dv_dz = 0.;
        dw_dx = 0.;       dw_dy = 0.;       dw_dz = 0.;      tSUPG_ = 0.;
        dumesh_dx = 0.;       dumesh_dy = 0.;       dumesh_dz = 0.; 
        dvmesh_dx = 0.;       dvmesh_dy = 0.;       dvmesh_dz = 0.;
        dwmesh_dx = 0.;       dwmesh_dy = 0.;       dwmesh_dz = 0.;
        fieldForce_.clear();     externalForces.clear();     rhsVector.clear();
        rhsVectorLM.clear();
        xK.clear();              XK.clear();                 di.clear();
        glueZone = false;        FSIInterface = false;
        intPointGlueZone.clear();
        sideBoundary_ = 0;

        for (int i=0; i < (sQuad.end() - sQuad.begin()); i++){
            intPointWeightFunction(i) = 1.;
        };

        setLocalNodes();
        getIntegPointCoordinates();

        DimVector xsi;        
        xsi.clear();        
        getJacobianMatrix(xsi);
    };

    //........................Element basic information.........................
    /// Clear all element variables
    void clearVariables();

    /// Sets the element connectivity
    /// @param Connectivity element connectivity
    void setConnectivity(Connectivity& connect){connect_ = connect;};

    /// Gets the element connectivity
    /// @return element connectivity
    Connectivity getConnectivity(){return connect_;};

    /// Sets the element density
    /// @param double element density
    void setDensity(double& dens){dens_ = dens;}

    /// Sets the element viscosity
    /// @param double element viscosity
    void setViscosity(double& visc){visc_ = visc;}

    /// Sets the time step size
    /// @param double time step size
    void setTimeStep(double& dt){dTime_ = dt;}

    /// Sets the time integration scheme
    /// @param double time integration scheme: 0.0 - Explicit forward Euler;
    /// 1.0 - Implicit backward Euler;
    /// 0.5 - Implicit Trapezoidal Rule.
    void setTimeIntegrationScheme(double& b){timeScheme_ = b;}

    /// Sets the body forces
    /// @param double* body forces
    void setFieldForce(double* ff);

    /// Sets the element vector of local nodes
    void setLocalNodes();

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void getJacobianMatrix(ublas::bounded_vector<double, DIM>& xsi);

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void getSpatialDerivatives(ublas::bounded_vector<double, DIM>& xsi);

    /// Compute and stores the interpolated velocities, mesh velocities, 
    /// previous mesh velocity, acceleration, among others
    void getVelAndDerivatives();

    /// Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getParameterSUPG();
    void getParameterSUPG2();

    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian(){return djac_;};

    /// Compute and store the velocity divergent
    void computeVelocityDivergent();

    /// Compute and store the drag and lift forces at the element boundary
    ublas::bounded_vector<double,DIM> getDragAndLiftForces();

    /// Compute and store the boundary forces
    ublas::bounded_vector<double,DIM> getBoundaryLoad(DimVector xsi);

    /// Gets the spatial jacobian matrix
    /// @param bounded_vector integration point coordinates
    /// @return spatial jacobian matrix
    DimMatrix getJacobianMatrixValues(ublas::bounded_vector<double, DIM>& xsi){
        getJacobianMatrix(xsi);
        return ainv_;};

    /// Gets the nodal gradient value for potential problem
    /// @param int direction @param int element node 
    /// @return nodal gradient value
    double getNodalGradientValue(int dir, int node){return gradient_(node,dir);}

    /// Sets the element side in boundary
    /// @param int side in boundary
    void setElemSideInBoundary(int side){sideBoundary_ = side;};

    /// Gets the element side in boundary
    /// @return side in boundary
    int getElemSideInBoundary(){return sideBoundary_;};

    /// Gets the results in a specific point for the potential problem
    /// @param VecLoc point adimensional coordinates
    /// @return potential results for the specific point
    std::pair<double, ublas::bounded_vector<double,DIM> >
    getPotentialResultsInPoint(VecLoc xsi);

    /// Sets the mesh moving weighting parameter for solving the Laplace problem
    /// @param double parameter value
    void setMeshMovingParameter(double value) {meshMovingParameter = value;};

    /// Gets the mesh moving weighting parameter
    /// @return mesh moving weighting parameter
    double getMeshMovingParameter(){return meshMovingParameter;};

    //.................Element intersection and correspondence..................
    /// Gets the element intersection parameters 
    /// @param minimum coordinates @param maximum coordinates 
    /// @param minimum inner product @param maximum inner product
    /// @param side lenght
    void setIntersectionParameters(DimVector x, DimVector X, 
                                   double *Dk, double *dk,
                           std::vector<ublas::bounded_vector<double,DIM> > dii);

    /// Gets the coordinates intersection parameters
    /// @return minimum and maximum coordinates
    std::pair<DimVector,DimVector> getXIntersectionParameter()
    {return std::make_pair(xK,XK);};

    /// Gets the inner product intersection parameters
    /// @return minimum and maximum inner products
    std::pair<double*,double*> getDIntersectionParameter()
    {return std::make_pair(dck,dCk);};

    /// Gets the side lenght intersection parameter
    /// @return side lenght intersection parameter
    std::vector<ublas::bounded_vector<double,DIM> >getDiIntersectionParameter()
    {return di;};

    //.............................Model functions..............................
    /// Sets if the element is in the gluing zone
    /// @param glueZone: if true is in the glue zone
    void setGlueZone(){glueZone = true;}

    /// Sets which model the fluid element belongs
    /// @param bool model: true = fine; false = coarse.
    void setModel(bool m){model = m;};

    /// Sets if the element belongs to the fluid structure interface
    void setFSIInterface(){FSIInterface = true;};

    /// Sets if the element must be incompressible
    /// If true, the element is compressible else it is incompressible
    void setCompressibility(){compressibility = true;};

    /// Clears the element compressibility condition setting false. 
    /// @see Element::setCompressibility()
    void clearCompressibility(){compressibility = false;};

    /// Gets the element compressibility constrain
    /// @return compressibility constrain @see Element::setCompressibility()
    bool getCompressibility(){return compressibility;};

    //......................Integration Points Information......................
    /// Gets the number of integration points of the special quadrature rule
    /// @retunr number of integration point of the special quadrature rule
    int getNumberOfIntegrationPoints(){return (sQuad.end() - sQuad.begin());};

    /// Sets the integration point correspondence to the overlapped mesh
    /// @param int element correspondent @param VecLoc Adimensional coordinates
    void setIntegrationPointCorrespondence(int elem, VecLoc x){
        intPointCorrespElem.push_back(elem);
        intPointCorrespXsi.push_back(x); };

    /// Gets the integration point correspondence - element
    /// @return overlapped element correspondence
    int getIntegPointCorrespondenceElement(int index)
    {return intPointCorrespElem[index];};

    /// Compute and store the integration points global coordinates
    void getIntegPointCoordinates();

    /// Gets the integration point global coordinates
    /// @param int integration point index @return integration point coordinates
    VecLoc getIntegPointCoordinatesValue(int index)
    {return intPointCoordinates[index];};

    /// Sets the integration point energy weight function
    /// @param int integration point index 
    /// @param double energy weight function value
    void setIntegPointWeightFunction(int index, double val)
    {intPointWeightFunction(index) = val;};

    /// Gets the integration point energy weight function
    /// @param int integration point index @return energy weight function value
    double getIntegPointWeightFunction(int index)
    {return intPointWeightFunction(index);};

    /// Sets if the integration point is the gluing zone
    /// @param int integration point index 
    void setIntegPointInGlueZone(int index){intPointGlueZone(index) = true;};

    /// Gets true if the integration point is in the gluing zone
    /// @param int integration point index @return If is in the gluing zone
    bool getIntegPointInGlueZone(int index){return intPointGlueZone(index);};

    /// Sets the integration point signaled distance function
    /// @param int integration point index
    /// @param signaled distance function valye
    void setIntegPointDistFunction(int index, double val)
    {intPointDistGlueZone(index) = val;};

    /// Gets the integration point signaled distance function value
    /// @return integration point signaled distance function value
    double getIntegPointDistFunction(int index)
    {return intPointDistGlueZone(index);};

    //.......................Element vectors and matrices.......................
    /// Compute and store the element matrix for the incompressible flow problem
    /// @param int integration point index
    void getElemMatrix(int index);

    /// Compute and store the element matrix for the Laplace/Poisson problem
    void getElemLaplMatrix();

    /// Sets the boundary conditions for the incompressible flow problem
    void setBoundaryConditions();

    /// Sets the boundary conditions for the Laplace/Poisson problem
    void setBoundaryConditionsLaplace();

    ///Compute and store the residual vector for the incompressible flow problem
    /// @param int integration point index
    void getResidualVector(int index);

    /// Compute and store the residual vector for the Laplace/Poisson problem
    void getResidualVectorLaplace();

    /// Gets the Newton-Raphson's jacobian matrix
    /// @return Newton-Raphson's jacobian matrix
    LocalMatrix getJacNRMatrix(){return jacobianNRMatrix;};

    /// Gets the Lagrange multiplier operator matrix
    /// @return Lagrange multiplier operator matrix
    LocalMatrix getLagrMultMatrix(){return lagrMultMatrix;};

    /// Gets the residual vector
    /// @return residual vecot
    LocalVector getRhsVector(){return rhsVector;};

    /// Gets the Lagrange multiplier residual vector
    /// @return Lagrange multiplier residual vector
    LocalVector getRhsVectorLagMult(){return rhsVectorLM;};

    /// Apply the boundary conditions and returns the matrix and residual vector
    /// for the Lagrange multiplier operator matrix, used when computing the 
    /// operator term from different meshes
    /// @param LocalMatrix Lagrange multiplier operator matrix
    /// @return residual vector and Lagrange multiplier operator matrix with
    /// boundary conditions applied
    std::pair<LocalVector,LocalMatrix> 
    getRhsVectorAndBoundaryConditions(LocalMatrix Ajac){
        jacobianNRMatrix = Ajac;
        setBoundaryConditionsLagrangeMultipliers();       
        return std::make_pair(rhsVector,jacobianNRMatrix);    };

    /// Sets the boundary conditions to the Lagrange multiplier operator
    void setBoundaryConditionsLagrangeMultipliers();

    /// Compute and store the nodal gradient value (for potential problems)
    void computeNodalGradient();

    /// Compute and store the Lagrange multiplier operator when integrating 
    /// the same mesh portion
    void getLagrangeMultipliersSameMesh();

    /// Compute and store the Lagrange multiplier operator when integrationg
    /// the different mesh portion
    /// @param int element of the coarse mesh (used to verify which integration
    /// point belongs to the coarse mesh element)
    void getLagrangeMultipliersDifferentMesh(int ielem);

    /// Compute and store the stabilization for the Lagrange multiplier operator
    /// for the same mesh portion
    void getLMStabilizationSameMesh();

    /// Compute and store the stabilization for the Lagrange multiplier operator
    /// for the different mesh portion
    /// @param int element of the coarse mesh 
    /// @see Element::getLagrangeMultipliersDifferentMesh(int ielem)
    void getStabCoarse(int ielem);

    //...............................Problem type...............................
    /// Compute the Steady Stokes problem matrices and vectors
    void getSteadyStokes();

    /// Compute the Steady Navier-Stokes problem matrices and vectors
    void getSteadyNavierStokes();

    /// Compute the Transient Stokes problem matrices and vectors
    void getTransientStokes();

    /// Compute the Transient Navier-Stokes problem matrices and vectors
    void getTransientNavierStokes();

    /// Compute the Steady Laplace problem matrices and vectors 
    /// (usually for the mesh moving step)
    void getSteadyLaplace();

    /// Compute the Transient Laplace problem matrices and vectors
    /// (usually for the mesh moving step)
    void getTransientLaplace();

};

template<>
double const Element<2>::k1 = 1.e0;
template<>
double const Element<2>::k2 = 0.0;


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setIntersectionParameters(DimVector x, DimVector X, 
                                           double *Dk, double *dk,
                           std::vector<ublas::bounded_vector<double,2> > dii) {
    xK = x; XK = X; 
    dCk[0] = Dk[0];
    dCk[1] = Dk[1]; 
    dCk[2] = Dk[2];  
    
    dck[0] = dk[0];
    dck[1] = dk[1];
    dck[2] = dk[2];   
    
    di = dii;
    
    return;
};

//------------------------------------------------------------------------------
//----------------------------SET ELEMENT FIELD FORCE---------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setFieldForce(double *ff) {
    fieldForce_(0) = ff[0];
    fieldForce_(1) = ff[1];
    
    return;
};

template<>
void Element<3>::setFieldForce(double *ff) {
    fieldForce_(0) = ff[0];
    fieldForce_(1) = ff[1];
    fieldForce_(2) = ff[2];
    
    return;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE INTEGRATION POINT COORDINATE-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getIntegPointCoordinates(){

    typename SpecialQuad::NodalValuesQuad nodalCoordx, nodalCoordy;
    VecLoc x; 
    
    for (int i = 0; i < 6; i++){
        nodalCoordx(i) = localNodes_(i,0);
        nodalCoordy(i) = localNodes_(i,1);
    };

    for (int i = 0; i < sQuad.end() - sQuad.begin(); i++){
        x(0) = 0.; x(1) = 0.;

        x(0) = sQuad.interpolateQuadraticVariable(nodalCoordx, i);
        x(1) = sQuad.interpolateQuadraticVariable(nodalCoordy, i);

        intPointCoordinates.push_back(x);
    };

    return;
};

template<>
void Element<3>::getIntegPointCoordinates(){

    typename SpecialQuad::NodalValuesQuad nodalCoordx, nodalCoordy, nodalCoordz;
    VecLoc x; 
    
    for (int i = 0; i < 10; i++){
        nodalCoordx(i) = localNodes_(i,0);
        nodalCoordy(i) = localNodes_(i,1);
    };

    for (int i = 0; i < sQuad.end() - sQuad.begin(); i++){
        x(0) = 0.; x(1) = 0.; x(2) = 0.;

        x(0) = sQuad.interpolateQuadraticVariable(nodalCoordx, i);
        x(1) = sQuad.interpolateQuadraticVariable(nodalCoordy, i);
        x(2) = sQuad.interpolateQuadraticVariable(nodalCoordz, i);

        intPointCoordinates.push_back(x);
    };

    return;
};

//------------------------------------------------------------------------------
//------------------------SEARCH THE ELEMENT LOCAL NODES------------------------
//------------------------------------------------------------------------------
//Creates a vector with the local nodes coordinates
template<>
void Element<2>::setLocalNodes() {
    
    typename Nodes::VecLocD x;

    for (int i=0; i<6; i++){
        x = nodes_[connect_(i)] -> getCoordinates();
        localNodes_(i,0) = x(0);
        localNodes_(i,1) = x(1);
    };

    getIntegPointCoordinates();

    return;
};

template<>
void Element<3>::setLocalNodes() {
    
    typename Nodes::VecLocD x;

    for (int i=0; i<10; i++){
        x = nodes_[connect_(i)] -> getCoordinates();
        localNodes_(i,0) = x(0);
        localNodes_(i,1) = x(1);
        localNodes_(i,2) = x(2);
    };

    return;
};

//------------------------------------------------------------------------------
//---------------------------CLEAR ELEMENT VARIABLES----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::clearVariables(){

    djac_ = 0.;               weight_ = 0.;              ainv_.clear();
    
    diffMatrix.clear();        gradient_.clear();
    laplMatrix.clear();
    jacobianNRMatrix.clear();  lagrMultMatrix.clear();
    
    u_ = 0.;          v_ = 0.;          w_ = 0;          p_ = 0.; 
    umesh_ = 0.;      vmesh_ = 0.;      wmesh_ = 0.; 
    du_dx = 0.;       du_dy = 0.;       du_dz = 0.; 
    dv_dx = 0.;       dv_dy = 0.;       dv_dz = 0.;
    dw_dx = 0.;       dw_dy = 0.;       dw_dz = 0.; 
    tSUPG_ = 0.;      tPSPG_ = .0;      tLSIC_ = 0.;
    
    externalForces.clear();   rhsVector.clear();
    localNodes_.clear();      rhsVectorLM.clear();

    xK.clear();              XK.clear();                 di.clear();

    ddphi_dx.clear();        dphi_dx.clear();            phi_.clear();     

    intPointCorrespXsi.clear();    intPointCoordinates.clear();
    intPointCorrespElem.clear();

    glueZone = false;        compressibility = false;
    intPointGlueZone.clear();

    for (int i=0; i < (sQuad.end() - sQuad.begin()); i++){
        intPointWeightFunction(i) = 1.;
    };
    
    setLocalNodes();
    getIntegPointCoordinates();

    return;
}; 

template<>
void Element<3>::clearVariables(){

    djac_ = 0.;               weight_ = 0.;              ainv_.clear();
    
    diffMatrix.clear();        gradient_.clear();
    laplMatrix.clear();
    jacobianNRMatrix.clear();  lagrMultMatrix.clear();
    
    u_ = 0.;          v_ = 0.;          w_ = 0;          p_ = 0.; 
    umesh_ = 0.;      vmesh_ = 0.;      wmesh_ = 0.; 
    du_dx = 0.;       du_dy = 0.;       du_dz = 0.; 
    dv_dx = 0.;       dv_dy = 0.;       dv_dz = 0.;
    dw_dx = 0.;       dw_dy = 0.;       dw_dz = 0.; 
    tSUPG_ = 0.;      tPSPG_ = 0.;      tLSIC_ = 0.;
    
    externalForces.clear();   rhsVector.clear();  rhsVectorLM.clear();
    
    glueZone = false;
    intPointGlueZone.clear();
    
    for (int i=0; i < (sQuad.end() - sQuad.begin()); i++){
        intPointWeightFunction(i) = 1.;
    };

    setLocalNodes();
    getIntegPointCoordinates();

    return;
}; 

//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getJacobianMatrix(ublas::bounded_vector<double,2>& xsi) {

    //Computes the spatial Jacobian matrix and its inverse

    typename QuadShapeFunction<2>::ValueDeriv dphi;
    
    shapeQuad.evaluateGradient(xsi,dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    for (int i=0; i<6; i++){
        dx_dxsi1 += localNodes_(i,0) * dphi(0,i);
        dx_dxsi2 += localNodes_(i,0) * dphi(1,i);
        dy_dxsi1 += localNodes_(i,1) * dphi(0,i);
        dy_dxsi2 += localNodes_(i,1) * dphi(1,i);        
    };

    //Computing the jacobian determinant
    djac_ = dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1;
    //    djac_ = fabs(djac_);

    //Computing Jacobian inverse
    ainv_(0,0) =  dy_dxsi2 / djac_;
    ainv_(0,1) = -dy_dxsi1 / djac_;
    ainv_(1,0) = -dx_dxsi2 / djac_;
    ainv_(1,1) =  dx_dxsi1 / djac_;

    return;
    
};

template<>
void Element<3>::getJacobianMatrix(ublas::bounded_vector<double,3>& xsi) {

    //Computes the spatial Jacobian matrix and its inverse

    typename QuadShapeFunction<3>::ValueDeriv dphi;
    
    shapeQuad.evaluateGradient(xsi,dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dx_dxsi3 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;
    double dy_dxsi3 = 0.;
    double dz_dxsi1 = 0.;
    double dz_dxsi2 = 0.;
    double dz_dxsi3 = 0.;

    for (int i=0; i<10; i++){
        dx_dxsi1 += localNodes_(i,0) * dphi(0,i);
        dx_dxsi2 += localNodes_(i,0) * dphi(1,i);
        dx_dxsi3 += localNodes_(i,0) * dphi(2,i);
        dy_dxsi1 += localNodes_(i,1) * dphi(0,i);
        dy_dxsi2 += localNodes_(i,1) * dphi(1,i);        
        dy_dxsi3 += localNodes_(i,1) * dphi(2,i);       
        dz_dxsi1 += localNodes_(i,2) * dphi(0,i);
        dz_dxsi2 += localNodes_(i,2) * dphi(1,i);        
        dz_dxsi3 += localNodes_(i,2) * dphi(2,i);        
    };
    
    //Computing the jacobian determinant
    djac_ = dx_dxsi1 * dy_dxsi2 * dz_dxsi3 - dx_dxsi1 * dy_dxsi3 * dz_dxsi2 
          + dx_dxsi2 * dy_dxsi3 * dz_dxsi1 - dx_dxsi2 * dy_dxsi1 * dz_dxsi3 
          + dx_dxsi3 * dy_dxsi1 * dz_dxsi2 - dx_dxsi3 * dy_dxsi2 * dz_dxsi1;

    //Computing Jacobian inverse
    ainv_(0,0) = (dy_dxsi2 * dz_dxsi3 - dy_dxsi3 * dz_dxsi2) / djac_;
    ainv_(0,1) = (dx_dxsi3 * dz_dxsi2 - dx_dxsi2 * dz_dxsi3) / djac_;
    ainv_(0,2) = (dx_dxsi2 * dy_dxsi3 - dx_dxsi3 * dy_dxsi2) / djac_;
    ainv_(1,0) = (dy_dxsi3 * dz_dxsi1 - dy_dxsi1 * dz_dxsi3) / djac_;
    ainv_(1,1) = (dx_dxsi1 * dz_dxsi3 - dx_dxsi3 * dz_dxsi1) / djac_;
    ainv_(1,2) = (dx_dxsi3 * dy_dxsi1 - dx_dxsi1 * dy_dxsi3) / djac_;
    ainv_(2,0) = (dy_dxsi1 * dz_dxsi2 - dy_dxsi2 * dz_dxsi1) / djac_;
    ainv_(2,1) = (dx_dxsi2 * dz_dxsi1 - dx_dxsi1 * dz_dxsi2) / djac_;
    ainv_(2,2) = (dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1) / djac_;

    return;
    
};

//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSpatialDerivatives(ublas::bounded_vector<double,2>& xsi) {
    
    typename QuadShapeFunction<2>::ValueDDeriv ddphi;
    typename QuadShapeFunction<2>::ValueDeriv  dphi;
    
    shapeQuad.evaluateGradient(xsi,dphi);
    shapeQuad.evaluateHessian(xsi,ddphi);
    
    dphi_dx.clear();
    ddphi_dx.clear();

    //Quadratic shape functions spatial first derivatives
    noalias(dphi_dx) = prod(ainv_,dphi);

    //Quadratic shape functions spatial second derivatives
    for (int i=0; i<6; i++){
        for (int j=0; i<2; i++){
            for (int k=0; i<2; i++){
                for (int l=0; i<2; i++){
                    for (int m=0; i<2; i++){                 
                        ddphi_dx(j,l)(i) += ainv_(j,k) * ainv_(l,m) * 
                                            ddphi(k,m)(i);
                    };
                };
            };
        };
    };

    return;
};

template<>
void Element<3>::getSpatialDerivatives(ublas::bounded_vector<double,3>& xsi) {

    typename QuadShapeFunction<3>::ValueDDeriv ddphi;    
    typename QuadShapeFunction<3>::ValueDeriv  dphi;

    shapeQuad.evaluateGradient(xsi,dphi);
    shapeQuad.evaluateHessian(xsi,ddphi);

    dphi_dx.clear();
    ddphi_dx.clear();

    //Quadratic shape functions spatial first derivatives
    for (int i=0; i<10; i++){
        dphi_dx(0,i) += ainv_(0,0) * dphi(0,i) 
                      + ainv_(1,0) * dphi(1,i) 
                      + ainv_(2,0) * dphi(2,i);
 
        dphi_dx(1,i) += ainv_(0,1) * dphi(0,i) 
                      + ainv_(1,1) * dphi(1,i) 
                      + ainv_(2,1) * dphi(2,i);         
 
       dphi_dx(2,i) += ainv_(0,2) * dphi(0,i) 
                      + ainv_(1,2) * dphi(1,i) 
                      + ainv_(2,2) * dphi(2,i);         
    };

    //Quadratic shape functions spatial second derivatives
    for (int i=0; i<10; i++){
        for (int j=0; i<3; i++){
            for (int k=0; i<3; i++){
                for (int l=0; i<3; i++){
                    for (int m=0; i<3; i++){                 
                        ddphi_dx(j,l)(i) += ainv_(j,k) * ainv_(l,m) * 
                                            ddphi(k,m)(i);
                    };
                };
            };
        };
    };

    return;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getVelAndDerivatives() {

    u_ = 0.;
    v_ = 0.;

    ax_ = 0.;
    ay_ = 0.;

    uPrev_ = 0.;
    vPrev_ = 0.;

    pPrev_ = 0.;

    umesh_ = 0.;
    vmesh_ = 0.;

    p_ = 0.;

    du_dx = 0.;
    du_dy = 0.;
    dv_dx = 0.;
    dv_dy = 0.;

    dp_dx = 0.;
    dp_dy = 0.;

    lagMx_ = 0.;
    lagMy_ = 0.;

    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){
        u_ += nodes_[connect_(i)] -> getVelocity(0) * phi_(i);
        v_ += nodes_[connect_(i)] -> getVelocity(1) * phi_(i);

        uPrev_ += nodes_[connect_(i)] -> getPreviousVelocity(0) * phi_(i);
        vPrev_ += nodes_[connect_(i)] -> getPreviousVelocity(1) * phi_(i);

        // ax_ += nodes_[connect_(i)] -> getAcceleration(0) * phi_(i);
        // ay_ += nodes_[connect_(i)] -> getAcceleration(1) * phi_(i);
        
        du_dx += nodes_[connect_(i)] -> getVelocity(0) * dphi_dx(0,i);
        du_dy += nodes_[connect_(i)] -> getVelocity(0) * dphi_dx(1,i);
        dv_dx += nodes_[connect_(i)] -> getVelocity(1) * dphi_dx(0,i);
        dv_dy += nodes_[connect_(i)] -> getVelocity(1) * dphi_dx(1,i);

        p_ += nodes_[connect_(i)] -> getPressure() * phi_(i);
        
        pPrev_ += nodes_[connect_(i)] -> getPreviousPressure() * phi_(i);

        dp_dx += nodes_[connect_(i)] -> getPressure() * dphi_dx(0,i);
        dp_dy += nodes_[connect_(i)] -> getPressure() * dphi_dx(1,i);

        umesh_ += nodes_[connect_(i)] -> getMeshVelocity(0) * phi_(i);
        vmesh_ += nodes_[connect_(i)] -> getMeshVelocity(1) * phi_(i);
        
        dumesh_dx += nodes_[connect_(i)] -> getMeshVelocity(0) * dphi_dx(0,i);
        dumesh_dy += nodes_[connect_(i)] -> getMeshVelocity(0) * dphi_dx(1,i);
        dvmesh_dx += nodes_[connect_(i)] -> getMeshVelocity(1) * dphi_dx(0,i);
        dvmesh_dy += nodes_[connect_(i)] -> getMeshVelocity(1) * dphi_dx(1,i);
        // std::cout << "velocity " << phi_(i) << " " << nodes_[connect_(i)] -> getVelocity(0) << std::endl;

        lagMx_ += nodes_[connect_(i)] -> getLagrangeMultiplier(0) * phi_(i);
        lagMy_ += nodes_[connect_(i)] -> getLagrangeMultiplier(1) * phi_(i);

    };  

    ax_ = (u_ - uPrev_) / dTime_;
    ay_ = (v_ - vPrev_) / dTime_;

    return;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
ublas::bounded_vector<double,2>
Element<2>::getBoundaryLoad(DimVector xsi) {

    setLocalNodes();

    ublas::bounded_vector<double,2>            load;
    DimMatrix                                  shearStress;
    typename QuadShapeFunction<2>::ValueDeriv  dphi;
    DimVector                                  t_vector, n_vector;
    ublas::identity_matrix<double>             ident(2);

    //Computes the shape functions        
    shapeQuad.evaluate(xsi,phi_);

    //Computes the jacobian matrix
    getJacobianMatrix(xsi);
    
    //Computes spatial derivatives
    getSpatialDerivatives(xsi);
    
    shapeQuad.evaluateGradient(xsi,dphi); 

    //Interpolates velocity and its derivatives values
    getVelAndDerivatives();
   
   
    shearStress(0,0) = 2. * visc_ * du_dx;
    shearStress(0,1) = visc_ * (du_dy + dv_dx);
    shearStress(1,0) = visc_ * (du_dy + dv_dx);
    shearStress(1,1) = 2. * visc_ * dv_dy;
    
    t_vector.clear();
    n_vector.clear();

    if (sideBoundary_ == 1){
        for (int i = 0; i < 6; i++){
            t_vector(0) -= dphi(1,i) * localNodes_(i,0);
            t_vector(1) -= dphi(1,i) * localNodes_(i,1);
        };        
    };

    if (sideBoundary_ == 2){
        for (int i = 0; i < 6; i++){
            t_vector(0) += dphi(0,i) * localNodes_(i,0);
            t_vector(1) += dphi(0,i) * localNodes_(i,1);
        };        
    };

    if (sideBoundary_ == 0){
        std::cout << "VERIFICAR VETOR NORMAL - getBoundaryLoad" << std::endl;
    };

    n_vector(0) =  t_vector(1) / norm_2(t_vector);
    n_vector(1) = -t_vector(0) / norm_2(t_vector);

    load.clear();
    load = -p_ * prod(ident,n_vector) + prod(shearStress,n_vector);

    //std::cout << "N Vector " << sideBoundary_ << " " <<  n_vector(0) << " " << n_vector(1) << " " << load(0) << " " << load(1) << " " << p_ << std::endl;
    
    return load;
};


//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
ublas::bounded_vector<double,2> Element<2>::getDragAndLiftForces() {

    setLocalNodes();

    ublas::bounded_vector<double, 3> nodesb_; 
    if(sideBoundary_ == 0){
        nodesb_(0) = connect_(1); 
        nodesb_(1) = connect_(4); 
        nodesb_(2) = connect_(2); 
        for (int i=0; i<2; i++){
            localNodesBoundary_(0,i) = localNodes_(1,i);
            localNodesBoundary_(1,i) = localNodes_(4,i);
            localNodesBoundary_(2,i) = localNodes_(2,i);
        };
    }else{
        if(sideBoundary_ == 1){
            nodesb_(0) = connect_(2); 
            nodesb_(1) = connect_(5); 
            nodesb_(2) = connect_(0); 
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = localNodes_(2,i);
                localNodesBoundary_(1,i) = localNodes_(5,i);
                localNodesBoundary_(2,i) = localNodes_(0,i);
            };
        }else{
            nodesb_(0) = connect_(0);
            nodesb_(1) = connect_(3);
            nodesb_(2) = connect_(1);
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = localNodes_(0,i);
                localNodesBoundary_(1,i) = localNodes_(3,i);
                localNodesBoundary_(2,i) = localNodes_(1,i);
            };
        };        
    };


    BoundaryQuad           bQuad;     //Boundary Integration Quadrature
    std::pair<BoundaryQuad::PointCoord,BoundaryQuad::PointWeight> gaussQuad;
    DimVector                                  n_vector;
    DimMatrix                                  shearStress;
    ublas::identity_matrix<double>             ident(2);
    ublas::bounded_vector<double,2>            load;
    
    load.clear();    
    n_vector.clear();
    gaussQuad.first.clear();
    gaussQuad.second.clear();
    shearStress.clear();

    gaussQuad = bQuad.GaussQuadrature();
    
    int index = 0;
    for(typename BoundaryQuad::QuadratureListIt it = bQuad.begin(); 
        it != bQuad.end(); it++){
        
        phib_ = shapeBound.getShapeFunction(gaussQuad.first(index));
        dphib_ = shapeBound.getShapeFunctionDerivative(gaussQuad.first(index));

        double Tx=0.; double Ty = 0.;
        
        for (int i=0; i<3; i++){
            Tx += localNodesBoundary_(i,0) * dphib_(i);
            Ty += localNodesBoundary_(i,1) * dphib_(i);
        };

        double jacb_ = sqrt(Tx*Tx + Ty*Ty);
        
        n_vector(0) =  Ty / jacb_;
        n_vector(1) = -Tx / jacb_;
        
        p_ = 0.;
        du_dx = 0.; du_dy = 0.; dv_dx = 0.; dv_dy = 0.;
        for (int i=0; i<3; i++){
            p_ += nodes_[nodesb_(i)] -> getPressure() * phib_(i);

            du_dx += nodes_[nodesb_(i)] -> getVelocity(0) * Tx;
            du_dy += nodes_[nodesb_(i)] -> getVelocity(0) * Ty;
            dv_dx += nodes_[nodesb_(i)] -> getVelocity(1) * Tx;
            dv_dy += nodes_[nodesb_(i)] -> getVelocity(1) * Ty;
        };
        
        shearStress(0,0) = 2. * visc_ * du_dx;
        shearStress(0,1) = visc_ * (du_dy + dv_dx);
        shearStress(1,0) = visc_ * (du_dy + dv_dx);
        shearStress(1,1) = 2. * visc_ * dv_dy;
        
        load += (-p_ * prod(ident,n_vector) + prod(shearStress,n_vector)) 
            * jacb_ * gaussQuad.second(index);

        //  std::cout << "Normal Vector " << n_vector(0) <<" " << n_vector(1) << std::endl;
        
        index++;

    };
    
    return load;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getParameterSUPG() {

    DimVector   r, s;

    tSUPG_ = 0.;
    tSUGN1_ = 0.;
    tSUGN2_ = 0.;
    tSUGN3_ = 0.;
    hRGN_ = 0.;
    double hUGN_ = 0.;
    
    r.clear(); s.clear();

    double u__ = 0.;
    double v__ = 0.;

    for (int i = 0; i < 6; i++){
        double ua = nodes_[connect_(i)] -> getVelocity(0);
        double va = nodes_[connect_(i)] -> getVelocity(1);
        
        u__ += ua * phi_(i);
        v__ += va * phi_(i);
    };

    double uNorm = sqrt(u__ * u__ + v__ * v__);
    if(uNorm > 1.e-10){
        s(0) = u__ / uNorm;
        s(1) = v__ / uNorm;
    }else{
        s(0) = 1. / sqrt(2.);
        s(1) = 1. / sqrt(2.);
    };

    for (int i = 0; i < 6; i++){
        double ua = nodes_[connect_(i)] -> getVelocity(0);
        double va = nodes_[connect_(i)] -> getVelocity(1);
        
        r(0) += sqrt(ua * ua + va * va) * dphi_dx(0,i);
        r(1) += sqrt(ua * ua + va * va) * dphi_dx(1,i);
    };

    double rNorm = norm_2(r);

    if (rNorm >= 1.e-10){
        r /= rNorm;
    }else{
        r(0) = 1. / sqrt(2.);
        r(1) = 1. / sqrt(2.);
    };
    
    for (int i = 0; i < 6; i++){
        hRGN_ += fabs(r(0) * dphi_dx(0,i) + r(1) * dphi_dx(1,i));
        hUGN_ += fabs(s(0) * dphi_dx(0,i) + s(1) * dphi_dx(1,i));        
    };

    if (hRGN_ >= 1.e-10){
        hRGN_ = 2. / hRGN_;
    }else{
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10){
        hUGN_ = 2. / hUGN_;
    }else{
        hUGN_ = 2. / 1.e-10;
    };    

    if (uNorm >= 1.e-10){
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }else{
        tSUGN1_ = hUGN_ / 2.e-10;
    };
              
    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

   
    if (fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;

 
    //Computing tSUPG parameter
    tSUPG_ = 1. / sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3_ * tSUGN3_));
   
    tPSPG_ = tSUPG_;
    tLSIC_ = tSUPG_ * uNorm;
    tARLQ_ = tSUPG_;
 
    return;
};

//------------------------------------------------------------------------------
//----------------------ELEMENT DIFFUSION/VISCOSITY MATRIX----------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getElemMatrix(int index){
    
    tARLQ_ = 0.;
    //tSUPG_ = 0.;
    //tPSPG_ = 0.;
    //tLSIC_ = 0.;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            double acelx = (u_ - uPrev_) / dTime_;
            double acely = (v_ - vPrev_) / dTime_;
            double wSUPGi = (u_ - umesh_) * dphi_dx(0,i) +
                (v_ - vmesh_) * dphi_dx(1,i);
            double wSUPGj = (u_ - umesh_) * dphi_dx(0,j) +      
                (v_ - vmesh_) * dphi_dx(1,j);
            
            //Mass matrix (use for both directions)
            double mM =  phi_(i) * phi_(j) * dens_ + 
                wSUPGi * phi_(j) * tSUPG_ * dens_;
            double sMxx = dphi_dx(0,i) * phi_(j) * acelx * tSUPG_ * dens_;
            double sMxy = dphi_dx(1,i) * phi_(j) * acelx * tSUPG_ * dens_;
            double sMyx = dphi_dx(0,i) * phi_(j) * acely * tSUPG_ * dens_;
            double sMyy = dphi_dx(1,i) * phi_(j) * acely * tSUPG_ * dens_;
                       
            //Difusion matrix (viscosity)
            double Kxx = 2. * dphi_dx(0,i) * dphi_dx(0,j) * visc_ + 
                dphi_dx(1,i) * dphi_dx(1,j) * visc_;
            double Kxy = dphi_dx(1,i) * dphi_dx(0,j) * visc_;
            double Kyx = dphi_dx(0,i) * dphi_dx(1,j) * visc_;
            double Kyy = 2. * dphi_dx(1,i) * dphi_dx(1,j) * visc_ + 
                dphi_dx(0,i) * dphi_dx(0,j) * visc_;
            
            //Convection matrix
            double Cxx = (dphi_dx(0,j) * (u_ - umesh_) + 
                          dphi_dx(1,j) * (v_ - vmesh_)) * phi_(i) * dens_
                + wSUPGi * wSUPGj * tSUPG_ * dens_;
            
            double Cyy  = Cxx;
            
            double Cuu = phi_(i) * du_dx * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * du_dx * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(0,i) * phi_(j) * 
                (u_ * du_dx + v_ * du_dy) * dens_;
            
            double Cuv = phi_(i) * du_dy * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * du_dy * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(1,i) * phi_(j) * 
                (u_ * du_dx + v_ * du_dy) * dens_;
            
            double Cvu = phi_(i) * dv_dx * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * dv_dx * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(0,i) * phi_(j) * 
                (u_ * dv_dx + v_ * dv_dy) * dens_;
            
            double Cvv = phi_(i) * dv_dy * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * dv_dy * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(1,i) * phi_(j) * 
                (u_ * dv_dx + v_ * dv_dy) * dens_;
            
            double Quu = dphi_dx(0,i) * dp_dx * phi_(j) * tSUPG_;
            double Quv = dphi_dx(1,i) * dp_dx * phi_(j) * tSUPG_;
            double Qvu = dphi_dx(0,i) * dp_dy * phi_(j) * tSUPG_;
            double Qvv = dphi_dx(1,i) * dp_dy * phi_(j) * tSUPG_;
            
            double KLSxx = dphi_dx(0,i) * dphi_dx(0,j) * tLSIC_ * dens_;
            double KLSxy = dphi_dx(0,i) * dphi_dx(1,j) * tLSIC_ * dens_;
            double KLSyx = dphi_dx(1,i) * dphi_dx(0,j) * tLSIC_ * dens_;
            double KLSyy = dphi_dx(1,i) * dphi_dx(1,j) * tLSIC_ * dens_;

            
            jacobianNRMatrix(2*i  ,2*j  ) += (mM + timeScheme_ * dTime_ * 
                                              (sMxx + Kxx + KLSxx / 
                                               timeScheme_ + Cxx + Cuu + Quu))
                * weight_ * djac_ * intPointWeightFunction(index);

            jacobianNRMatrix(2*i+1,2*j+1) += (mM + timeScheme_ * dTime_ *
                                              (sMyy + Kyy + KLSyy / 
                                               timeScheme_ + Cyy + Cvv + Qvv))
                * weight_ * djac_ * intPointWeightFunction(index);

            jacobianNRMatrix(2*i  ,2*j+1) += (timeScheme_ * dTime_ * 
                (Kxy + sMxy + Cuv + Quv) + KLSxy * dTime_)
                * weight_ * djac_ * intPointWeightFunction(index);

            jacobianNRMatrix(2*i+1,2*j  ) += (timeScheme_ * dTime_ *
                (Kyx + sMyx + Cvu + Qvu) + KLSyx * dTime_)
                * weight_ * djac_ * intPointWeightFunction(index); 

            
            //SINAL DA PARCELA QUE MULTIPLICA O TSUPG ESTA
            //COM SINAL TROCADO NA FORMULAÇAO DO TEZDUYAR
            //multipy pressure direction x
            double QSUPGx = - (dphi_dx(0,i) * phi_(j) + 
                               (dphi_dx(0,i) * (u_ - umesh_) + 
                                dphi_dx(1,i) * (v_ - vmesh_)) * 
                               dphi_dx(0,j) * tSUPG_);
            //multiply pressure direction y
            double QSUPGy = - (dphi_dx(1,i) * phi_(j) +
                               (dphi_dx(0,i) * (u_ - umesh_) + 
                                dphi_dx(1,i) * (v_ - vmesh_)) *   
                               dphi_dx(1,j) * tSUPG_);
            //multiply velocity direction x
            double Qx = dphi_dx(0,i) * phi_(j);
            //multiply velocity direction y
            double Qy = dphi_dx(1,i) * phi_(j);
            
            
            jacobianNRMatrix(12+j,2*i  ) += Qx * dTime_ * weight_ * djac_ 
                * intPointWeightFunction(index);
            jacobianNRMatrix(12+j,2*i+1) += Qy * dTime_ * weight_ * djac_ 
                * intPointWeightFunction(index);
            jacobianNRMatrix(2*i  ,12+j) += QSUPGx * dTime_ * weight_ * djac_
                * intPointWeightFunction(index);
            jacobianNRMatrix(2*i+1,12+j) += QSUPGy * dTime_ * weight_ * djac_
                * intPointWeightFunction(index);


            double Hx = dphi_dx(0,i) * phi_(j) * tPSPG_;
            double Hy = dphi_dx(1,i) * phi_(j) * tPSPG_;
            
            double Gx = dphi_dx(0,i) * ((u_ - umesh_) * dphi_dx(0,j) + 
                                        (v_ - vmesh_) * dphi_dx(1,j)) * 
                tPSPG_;
            double Gy = dphi_dx(1,i) * ((u_ - umesh_) * dphi_dx(0,j) + 
                                        (v_ - vmesh_) * dphi_dx(1,j)) * 
                tPSPG_;
            
            double Guu = (dphi_dx(0,i) * du_dx * phi_(j) + 
                          dphi_dx(1,i) * dv_dx * phi_(j)) * tPSPG_;
            double Gvv = (dphi_dx(0,i) * du_dy * phi_(j) + 
                          dphi_dx(1,i) * dv_dy * phi_(j)) * tPSPG_;
            
            double Q = (dphi_dx(0,i) * dphi_dx(0,j) + 
                        dphi_dx(1,i) * dphi_dx(1,j)) * tPSPG_ / (dens_);
            

            jacobianNRMatrix(12+j,2*i  ) += (Hx + (Gx + Guu) * 
                                              timeScheme_ * dTime_)
                * weight_ * djac_ * intPointWeightFunction(index);
            jacobianNRMatrix(12+j,2*i+1) += (Hy + (Gy + Gvv) * 
                                              timeScheme_ * dTime_)
                * weight_ * djac_ * intPointWeightFunction(index);
            jacobianNRMatrix(12+j,12+i) += Q * dTime_ * weight_ * djac_
                * intPointWeightFunction(index);
        };
    };

   

    return;
};


//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setBoundaryConditions(){

    for (int i = 0; i < 6; i++){

        if ((nodes_[connect_(i)] -> getConstrains(0) == 1) ||
            (nodes_[connect_(i)] -> getConstrains(0) == 3))  {
            for (int j = 0; j < 15; j++){
                jacobianNRMatrix(2*i  ,j) = 0.;
                jacobianNRMatrix(j,2*i  ) = 0.;
            };
            jacobianNRMatrix(2*i  , 2*i  ) = 1.;
            rhsVector(2*i  ) = 0.;
        };

        if ((nodes_[connect_(i)] -> getConstrains(1) == 1) ||
            (nodes_[connect_(i)] -> getConstrains(1) == 3)) {
            for (int j = 0; j < 15; j++){
                jacobianNRMatrix(2*i+1,j) = 0.;
                jacobianNRMatrix(j,2*i+1) = 0.;
            };
            jacobianNRMatrix(2*i+1, 2*i+1) = 1.;
            rhsVector(2*i+1) =  0.;
        };
    };

    // if (compressibility){
    //     jacobianNRMatrix(12,12) = 1.;
    //     jacobianNRMatrix(13,13) = 1.;
    //     jacobianNRMatrix(14,14) = 1.;

    //     rhsVector(12) = 0.;
    //     rhsVector(13) = 0.;
    //     rhsVector(14) = 0.;

    // };


    // typename Nodes::VecLocD x;

    // for (int i = 0; i < 6; i++){
    //     //if(model){
    //     x = nodes_[connect_(i)] -> getCoordinates();
    //     if((x(0) > 0.999) && (x(1) > 0.999)){
    //         // std::cout << "AQUI  " << index_ << std::endl;
            
    //         for (int j = 0; j < 18; j++){
    //             jacobianNRMatrix(12+i,j) = 0.;
    //             jacobianNRMatrix(j,12+i) = 0.;
    //         };
    //         jacobianNRMatrix(12+i,12+i) = 1.;
    //         rhsVector(12+i) =  0.;
    //         //};
    //     };
    // };

    return;
};

//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setBoundaryConditionsLaplace(){

    for (int i = 0; i < 6; i++){

        if (nodes_[connect_(i)] -> getConstrainsLaplace(0) == 1) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i  ,j) = 0.;
                //jacobianNRMatrix(j,2*i  ) = 0.;
            };
            jacobianNRMatrix(2*i  , 2*i  ) = 1.;
            //rhsVector(2*i  ) = 0.;
        };

        if (nodes_[connect_(i)] -> getConstrainsLaplace(1) == 1) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i+1,j) = 0.;
                //jacobianNRMatrix(j,2*i+1) = 0.;
            };
            jacobianNRMatrix(2*i+1, 2*i+1) = 1.;
            //rhsVector(2*i+1) =  0.;
        };
    };
    
    return;
};


//------------------------------------------------------------------------------
//---------------APPLY BOUNDARY CONDITIONS TO LAGRANGE MULTIPLIERS--------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setBoundaryConditionsLagrangeMultipliers(){

    rhsVector.clear();
    rhsVectorLM.clear();
    
    LocalVector U_;

    U_.clear();

    for (int i=0; i<6; i++){
        U_(2*i  ) = nodes_[connect_(i)] -> getVelocity(0);
        U_(2*i+1) = nodes_[connect_(i)] -> getVelocity(1);
    };

    noalias(rhsVector) = - prod(jacobianNRMatrix,U_); 
    
    // for (int i = 0; i < 6; i++){
    //     if (nodes_[connect_(i)] -> getConstrains(0) == 1) {
    //         for (int j = 0; j < 15; j++){
    //             //jacobianNRMatrix(2*i  ,j) = 0.;
    //             jacobianNRMatrix(j,2*i  ) = 0.;
    //         };
    //         //rhsVector(2*i  ) = 0.0;
    //     };

    //     if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
    //         // std::cout << "aqui" << std::endl;
    //         for (int j = 0; j < 15; j++){
    //             //jacobianNRMatrix(2*i+1,j) = 0.;
    //             jacobianNRMatrix(j,2*i+1) = 0.;
    //         };
    //         //rhsVector(2*i+1) = 0.0;
    //     };
    // };

    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getResidualVector(int index){
    
    double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;

    double duna_dx = 0.;
    double duna_dy = 0.;
    double dvna_dx = 0.;
    double dvna_dy = 0.;

    for (int i = 0; i < 6; i++){
        duna_dx += una_ * dphi_dx(0,i);
        duna_dy += una_ * dphi_dx(1,i);        
        dvna_dx += vna_ * dphi_dx(0,i);
        dvna_dy += vna_ * dphi_dx(1,i);
    };

    for (int i = 0; i < 6; i++){
        //!!!! MULTIPLICAR NO FIM PELO PESO DE GAUSS E JACOBIANO
        double mx = phi_(i) * (u_ - uPrev_) * dens_ + 
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) * 
            (u_ - uPrev_) * tSUPG_ * dens_;
        double my = phi_(i) * (v_ - vPrev_) * dens_ + 
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) * 
            (v_ - vPrev_) * tSUPG_ * dens_;
        
        double Kx = 2. * dphi_dx(0,i) * du_dx * visc_ + 
            dphi_dx(1,i) * du_dy * visc_ + 
            dphi_dx(1,i) * dv_dx * visc_;
        double Ky= dphi_dx(0,i) * du_dy * visc_ + 
            2. * dphi_dx(1,i) * dv_dy * visc_ + 
            dphi_dx(0,i) * dv_dx * visc_;            

        double KLSx = dphi_dx(0,i) * (du_dx + dv_dy) * tLSIC_ * dens_;
        double KLSy = dphi_dx(1,i) * (du_dx + dv_dy) * tLSIC_ * dens_;

        double Cx = (du_dx * (una_ - umesh_) + 
                     du_dy * (vna_ - vmesh_)) * phi_(i) * dens_ +
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) *
            ((una_ - umesh_) * du_dx + (vna_ - vmesh_) * du_dy) * tSUPG_ *dens_;
        double Cy = (dv_dx * (una_ - umesh_) + 
                     dv_dy * (vna_ - vmesh_)) * phi_(i) * dens_ +
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) *
            ((una_ - umesh_) * dv_dx + (vna_ - vmesh_) * dv_dy) * tSUPG_ *dens_;

        mx *= intPointWeightFunction(index);
        my *= intPointWeightFunction(index);
        Kx *= intPointWeightFunction(index);
        Ky *= intPointWeightFunction(index);
        KLSx *= intPointWeightFunction(index);
        KLSy *= intPointWeightFunction(index);
        Cx *= intPointWeightFunction(index);
        Cy *= intPointWeightFunction(index);
  
       
        double Px = - (dphi_dx(0,i) * p_) - ((dphi_dx(0,i) * (una_ - umesh_) +
                                              dphi_dx(1,i) * (vna_ - vmesh_))
                                             * dp_dx * tSUPG_);
        double Py = - (dphi_dx(1,i) * p_) - ((dphi_dx(0,i) * (una_ - umesh_) +
                                              dphi_dx(1,i) * (vna_ - vmesh_))
                                             * dp_dy * tSUPG_);
           
        Px *= intPointWeightFunction(index);
        Py *= intPointWeightFunction(index);

        double Q = ((du_dx + dv_dy) * phi_(i)) +
            (dphi_dx(0,i) * dp_dx + dphi_dx(1,i) * dp_dy) * tPSPG_ / dens_ +
            dphi_dx(0,i) * (u_ - uPrev_) / dTime_ * tPSPG_ +
            dphi_dx(1,i) * (v_ - vPrev_) / dTime_ * tPSPG_ +
            dphi_dx(0,i) * ((una_ - umesh_) * du_dx +
                            (vna_ - vmesh_) * du_dy) * tPSPG_ +
            dphi_dx(1,i) * ((una_ - umesh_) * dv_dx +
                            (vna_ - vmesh_) * dv_dy) * tPSPG_;

        Q *= intPointWeightFunction(index);

        rhsVector(2*i  ) += (-mx + (-Kx - Px - Cx) * dTime_ * timeScheme_ 
                             - KLSx * dTime_) * weight_ * djac_;
        rhsVector(2*i+1) += (-my + (-Ky - Py - Cy) * dTime_ * timeScheme_
                             - KLSy * dTime_) * weight_ * djac_;
        rhsVector(12+i) += -Q * dTime_ * weight_ * djac_;
                              
    };

    




    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------

template<>
void Element<2>::getResidualVectorLaplace(){

    rhsVector.clear();
    
    LocalVector U_;
    typename Nodes::VecLocD x_up;

    U_.clear();

    for (int i=0; i<6; i++){

        x_up = nodes_[connect_(i)] -> getUpdatedCoordinates();        
        
        if (nodes_[connect_(i)] -> getConstrainsLaplace(0) == 1){
            U_(2*i  ) = x_up(0) - localNodes_(i,0);
        };
        if (nodes_[connect_(i)] -> getConstrainsLaplace(1) == 1){
            U_(2*i+1) = x_up(1) - localNodes_(i,1);
        };

        //if(connect_(i) == 29) std::cout << "Element 29 " << x_up(0) << " " << x_up(1) <<std::endl;
    };
    
    noalias(rhsVector) += U_;//- prod(laplMatrix,U_); 

    return;
};

//------------------------------------------------------------------------------
//---------------------------ELEMENT LAPLACIAN MATRIX---------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getElemLaplMatrix(){

     for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){        
            laplMatrix(2*i  ,2*j  ) += (dphi_dx(0,i) * dphi_dx(0,j) +
                                        dphi_dx(1,i) * dphi_dx(1,j)) 
                                      * weight_ * djac_ * meshMovingParameter;
            laplMatrix(2*i+1,2*j+1) += (dphi_dx(0,i) * dphi_dx(0,j) +
                                        dphi_dx(1,i) * dphi_dx(1,j)) 
                                      * weight_ * djac_ * meshMovingParameter;
        };
    };
     
    return;
};

//------------------------------------------------------------------------------
//-------------------------COMPUTES NODAL GRADIENT VALUE------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeNodalGradient(){
    
    double dp_dx, dp_dy;
    typename QuadShapeFunction<2>::Coords xsi;
    double potential[6];

    for (int i=0; i<6; i++){
        potential[i] = nodes_[connect_(i)] -> getPotential();
    };
    
    //Node 0
    xsi(0) = 0.;
    xsi(1) = 0.;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(0)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(0)] -> setGradientComponent(-dp_dx, 1);

    gradient_(0,0) = dp_dy;
    gradient_(0,1) =-dp_dx;

    //Node 1
    xsi(0) = 1.;
    xsi(1) = 0.;

    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(1)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(1)] -> setGradientComponent(-dp_dx, 1);      

    gradient_(1,0) = dp_dy;
    gradient_(1,1) =-dp_dx;

    //Node 2
    xsi(0) = 0.;
    xsi(1) = 1.;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(2)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(2)] -> setGradientComponent(-dp_dx, 1);

    gradient_(2,0) = dp_dy;
    gradient_(2,1) =-dp_dx;

    //Node 3
    xsi(0) = .5;
    xsi(1) = 0.;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(3)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(3)] -> setGradientComponent(-dp_dx, 1);

    gradient_(3,0) = dp_dy;
    gradient_(3,1) =-dp_dx;

    //Node 4
    xsi(0) = .5;
    xsi(1) = .5;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(4)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(4)] -> setGradientComponent(-dp_dx, 1);

    gradient_(4,0) = dp_dy;
    gradient_(4,1) =-dp_dx;

    //Node 5
    xsi(0) = 0.;
    xsi(1) = .5;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(5)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(5)] -> setGradientComponent(-dp_dx, 1);
    
    gradient_(5,0) = dp_dy;
    gradient_(5,1) =-dp_dx;

    return;
};

//------------------------------------------------------------------------------
//---------COMPUTES VALUES OF POTENTIAL AND GRADIENT IN A SPECIFIC POINT--------
//------------------------------------------------------------------------------
template<>
std::pair<double,ublas::bounded_vector<double,2> > Element<2>::
                            getPotentialResultsInPoint(VecLoc xsi){
    
    shapeQuad.evaluate(xsi,phi_);
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    double potential[6];
    
    for (int i=0; i<6; i++){
        potential[i] = nodes_[connect_(i)] -> getPressure();
    };
    
    VecLoc dp_dx; 
    double pot = 0.;
    
    dp_dx.clear();
    
    for (int i=0; i<6; i++){
        dp_dx(0) += dphi_dx(1,i) * potential[i];
        dp_dx(1) -= dphi_dx(0,i) * potential[i];

        pot += potential[i] * phi_(i);
    };

    return std::make_pair(pot,dp_dx);

};

//------------------------------------------------------------------------------
//-----------------------COMPUTES THE VELOCITY DIVERGENT------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeVelocityDivergent(){
    
    ublas::bounded_matrix <double, 2, 6>      int_points;
    
    int_points(0,0) = 0.; int_points(1,0) = 0.;
    int_points(0,1) = 1.; int_points(1,1) = 0.;
    int_points(0,2) = 0.; int_points(1,2) = 1.;
    int_points(0,3) = .5; int_points(1,3) = 0.;
    int_points(0,4) = .5; int_points(1,4) = .5;
    int_points(0,5) = 0.; int_points(1,5) = .5;
    

    typename QuadShapeFunction<2>::Coords xsi;
    
    for (int i=0; i<6; i++){
        
        xsi(0) = int_points(0,i);
        xsi(1) = int_points(1,i);
        
        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);
        
        //Computes spatial derivatives
        getSpatialDerivatives(xsi);
        
        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();
        
        if(fabs(nodes_[connect_(i)] -> getVelocityDivergent()) <= 
           fabs(du_dx+dv_dy)){
            nodes_[connect_(i)] -> setVelocityDivergent(du_dx + dv_dy);        
        };
    };
    
    return;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//----------------MOUNTS EACH TYPE OF INCOMPRESSIBLE FLOW PROBEM----------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------STEADY STOKES PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSteadyStokes(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    rhsVector.clear();
    jacobianNRMatrix.clear();

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();
        
        //Computes the element diffusion/viscosity matrix
        getElemMatrix(index);
       
        //Computes the RHS vector
        getResidualVector(index);

        index++;        
    };  
    
    //Apply boundary conditions
    setBoundaryConditions();

    return;
};


//------------------------------------------------------------------------------
//-------------------------STEADY NAVIER-STOKES PROBEM--------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSteadyNavierStokes(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;
    
    rhsVector.clear();
    jacobianNRMatrix.clear();

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();

        //Computes the element diffusion/viscosity matrix
        getElemMatrix(index);

        //Computes the RHS vector
        getResidualVector(index);
        
        index++;        
    };  

    //Apply boundary conditions
    setBoundaryConditions();

    return;
};

//------------------------------------------------------------------------------
//---------------------------TRANSIENT STOKES PROBEM----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getTransientStokes(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    rhsVector.clear();
    jacobianNRMatrix.clear();

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();

        //Computes the element diffusion/viscosity matrix
        getElemMatrix(index);

        //Computes the RHS vector
        getResidualVector(index);
        
        index++;        
    };  
    
    //Apply boundary conditions
    setBoundaryConditions();

    return;
};

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getTransientNavierStokes(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    jacobianNRMatrix.clear();
    rhsVector.clear();
    setLocalNodes();


    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){

        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();

        //Compute Stabilization Parameters
        getParameterSUPG();

        //Computes the element diffusion/viscosity matrix
        getElemMatrix(index);

        //Computes the RHS vector
        getResidualVector(index); 

        index++;        
    };  
    
    //Apply boundary conditions
    setBoundaryConditions();

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSteadyLaplace(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    laplMatrix.clear();
    jacobianNRMatrix.clear();

    setLocalNodes();   
    
    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        getElemLaplMatrix();

        index++;        
    };  
    
    //Computes the Newton's method jacobian
    jacobianNRMatrix = laplMatrix;
    
    //Computes the RHS vector
    getResidualVectorLaplace();

    //Apply boundary conditions
    setBoundaryConditionsLaplace();

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getLagrangeMultipliersSameMesh(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    diffMatrix.clear();
    jacobianNRMatrix.clear();
    rhsVector.clear();
    lagrMultMatrix.clear();
    
    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        getSpatialDerivatives(xsi);
        
        getVelAndDerivatives();     
        
        ublas::bounded_matrix<double, 3,18> Bmatrix, Baux;
        ublas::bounded_matrix<double, 18,3> BmatrixT; 
        ublas::bounded_matrix<double, 3,3>  Maux;

        Bmatrix.clear(); BmatrixT.clear(); Maux.clear(); Baux.clear();

        for (int i = 0; i < 6; i++){
            Bmatrix(0,2*i  ) = dphi_dx(0,i);
            Bmatrix(2,2*i  ) = dphi_dx(1,i);       
            
            BmatrixT(2*i  ,0) = dphi_dx(0,i);
            BmatrixT(2*i  ,2) = dphi_dx(1,i);        
            
            Bmatrix(1,2*i+1) = dphi_dx(1,i);
            Bmatrix(2,2*i+1) = dphi_dx(0,i);
            
            BmatrixT(2*i+1,1) = dphi_dx(1,i);
            BmatrixT(2*i+1,2) = dphi_dx(0,i);
        };
        
        Maux(0,0) =  2.;
        Maux(1,1) =  2.;
        Maux(2,2) =  1.;
        
        Baux = prod(Maux,Bmatrix);
        
        //diffMatrix += prod(BmatrixT,Baux) * weight_ * djac_ * k2;
        
        for (int i = 0; i < 6; i++){
            for (int j = 0; j < 6; j++){        
                lagrMultMatrix(2*i  ,2*j  ) += 
                    (phi_(i) * phi_(j))
                    * weight_ * djac_ * k1;
                lagrMultMatrix(2*i+1,2*j+1) += 
                    (phi_(i) * phi_(j))
                    * weight_ * djac_ * k1;

                diffMatrix(2*i  ,2*j  ) += dphi_dx(0,i) * dphi_dx(0,j) *
                    weight_ * djac_ * k2;
                diffMatrix(2*i+1,2*j  ) += dphi_dx(1,i) * dphi_dx(0,j) *
                    weight_ * djac_ * k2;
                diffMatrix(2*i  ,2*j+1) += dphi_dx(0,i) * dphi_dx(1,j) *
                    weight_ * djac_ * k2;
                diffMatrix(2*i+1,2*j+1) += dphi_dx(1,i) * dphi_dx(1,j) *
                    weight_ * djac_ * k2;

                double LM = 0.;
    
                //Fine
                LM = - ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i))
                    * phi_(j) * tSUPG_;
            
                jacobianNRMatrix(2*i  ,2*j  ) += (timeScheme_ * dTime_ * LM)
                    * weight_ * djac_;
                
                jacobianNRMatrix(2*i+1,2*j+1) += (timeScheme_ * dTime_ * LM)
                    * weight_ * djac_;
        

                double Lx = 0.;
                double Ly = 0.;
                
                Lx = -dphi_dx(0,i) * phi_(j) * tPSPG_ / dens_;
                Ly = -dphi_dx(1,i) * phi_(j) * tPSPG_ / dens_;

                jacobianNRMatrix(2*i  ,12+j) += (Lx * timeScheme_ * dTime_)
                    * weight_ * djac_;
                jacobianNRMatrix(2*i+1,12+j) += (Ly * timeScheme_ * dTime_)
                    * weight_ * djac_;

            };

            double LMx = 0.;
            double LMy = 0.;
            
            LMx = - ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i))
                * (-lagMx_) * tSUPG_
                - dphi_dx(0,i) * (-lagMx_) * tPSPG_ / dens_;
            LMy = - ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i))
                * (-lagMy_) * tSUPG_
                - dphi_dx(0,i) * (-lagMy_) * tPSPG_ / dens_;

            rhsVector(2*i  ) += (-LMx * dTime_ * timeScheme_) * weight_ * djac_;
            rhsVector(2*i+1) += (-LMy * dTime_ * timeScheme_) * weight_ * djac_;


        };         
        
        index++;        
    }; 

    //jacobianNRMatrix.clear();
    // jacobianNRMatrix += diffMatrix;


    lagrMultMatrix += diffMatrix;

    //lagrMultMatrix.clear();


    // ublas::identity_matrix<double>             ident(18);
    // lagrMultMatrix = ident;

    // lagrMultMatrix(12,12) = 0.;
    // lagrMultMatrix(13,13) = 0.;
    // lagrMultMatrix(14,14) = 0.;

    // //Computes the Newton's method jacobian
    //jacobianNRMatrix = lagrMultMatrix;

    // Set boundary conditions to the Lagrange Multipliers
    //setBoundaryConditionsLagrangeMultipliers();

    rhsVectorLM.clear();  
    LocalVector U_;
    U_.clear();

    for (int i=0; i<6; i++){
        U_(2*i  ) = nodes_[connect_(i)] -> getVelocity(0);
        U_(2*i+1) = nodes_[connect_(i)] -> getVelocity(1);
    };

    noalias(rhsVectorLM) = - prod(lagrMultMatrix,U_); 


    //!!!!! O erro é que não pode somar aqui porque senão vai adicionar isso na equação dos multiplicadores também e só deve adicionar na parte da equação do momentum
    // lagrMultMatrix += jacobianNRMatrix;
    // rhsVectorLM += rhsVectorLM;


    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getLagrangeMultipliersDifferentMesh(int ielem){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    typename QuadShapeFunction<2>::Values phiLM_; 
   
    diffMatrix.clear();
    jacobianNRMatrix.clear();
    rhsVector.clear();
    lagrMultMatrix.clear();
    
    for(typename SpecialQuad::QuadratureListIt it = sQuad.begin(); 
        it != sQuad.end(); it++){
        
        if ((intPointCorrespElem[index] == ielem)){
            // std::cout << "Corresp " << intPointCorrespElem[index] << " " << index << " " << index_ << std::endl;
            //Defines the integration points adimentional coordinates
            xsi(0) = sQuad.PointList(index,0);
            xsi(1) = sQuad.PointList(index,1);
            
            //Computes the velocity shape functions
            shapeQuad.evaluate(xsi,phi_);
            
            //Computes the coarse mesh shape functions
            shapeQuad.evaluate(intPointCorrespXsi[index],phiLM_);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightList(index);
            
            //Computes the jacobian matrix
            getJacobianMatrix(xsi);
                        
            getSpatialDerivatives(xsi);

            dphiL_dx = dphi_dx;

            getSpatialDerivatives(intPointCorrespXsi[index]);

            ublas::bounded_matrix<double, 3,18> Bmatrix, Baux;
            ublas::bounded_matrix<double, 18,3> BmatrixT; 
            ublas::bounded_matrix<double, 3,3>  Maux;
            
            Bmatrix.clear(); BmatrixT.clear(); Maux.clear(); Baux.clear();
            
            for (int i = 0; i < 6; i++){
                Bmatrix(0,2*i  ) = dphi_dx(0,i);
                Bmatrix(2,2*i  ) = dphi_dx(1,i);       
                
                BmatrixT(2*i  ,0) = dphiL_dx(0,i);
                BmatrixT(2*i  ,2) = dphiL_dx(1,i);        
                
                Bmatrix(1,2*i+1) = dphi_dx(1,i);
                Bmatrix(2,2*i+1) = dphi_dx(0,i);
                
                BmatrixT(2*i+1,1) = dphiL_dx(1,i);
                BmatrixT(2*i+1,2) = dphiL_dx(0,i);
            };
            
            Maux(0,0) =  2.;
            Maux(1,1) =  2.;
            Maux(2,2) =  1.;
            
            Baux = prod(Maux,Bmatrix);
            
            //diffMatrix += prod(BmatrixT,Baux) * weight_ * djac_ * k2;
            
            for (int i = 0; i < 6; i++){
                for (int j = 0; j < 6; j++){        
                    lagrMultMatrix(2*i  ,2*j  ) += 
                        (phiLM_(i) * phi_(j))
                        * weight_ * djac_ * k1;
                    lagrMultMatrix(2*i+1,2*j+1) += 
                        (phiLM_(i) * phi_(j))
                        * weight_ * djac_ * k1;

                    diffMatrix(2*i  ,2*j  ) += dphiL_dx(0,i) * dphi_dx(0,j) *
                        weight_ * djac_ * k2;
                    diffMatrix(2*i+1,2*j  ) += dphiL_dx(1,i) * dphi_dx(0,j) *
                        weight_ * djac_ * k2;
                    diffMatrix(2*i  ,2*j+1) += dphiL_dx(0,i) * dphi_dx(1,j) *
                        weight_ * djac_ * k2;
                    diffMatrix(2*i+1,2*j+1) += dphiL_dx(1,i) * dphi_dx(1,j) *
                        weight_ * djac_ * k2;

                    //Coarse
                    double LM = 0.;
                    LM = ((u_ - umesh_) * phi_(j) + (v_ - vmesh_) * phi_(j))
                        * phiLM_(i) * tSUPG_;
                    
                    jacobianNRMatrix(2*i  ,2*j  ) += (timeScheme_ * dTime_ * LM)
                        * weight_ * djac_;
                    
                    jacobianNRMatrix(2*i+1,2*j+1) += (timeScheme_ * dTime_ * LM)
                        * weight_ * djac_;
                    
                    
                    double Lx = 0.;
                    double Ly = 0.;
                    
                    Lx = dphi_dx(0,j) * phiLM_(i) * tPSPG_ / dens_;
                    Ly = dphi_dx(1,j) * phiLM_(i) * tPSPG_ / dens_;
                    
                    jacobianNRMatrix(2*i  ,12+j) += (Lx * timeScheme_ * dTime_)
                        * weight_ * djac_;
                    jacobianNRMatrix(2*i+1,12+j) += (Ly * timeScheme_ * dTime_)
                        * weight_ * djac_;    
                    
                };

          
                double LMx = 0.;
                double LMy = 0.;
                
                LMx = - ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i))
                    * (lagMx_) * tSUPG_
                    - dphi_dx(0,i) * (lagMx_) * tPSPG_ / dens_;
                LMy = - ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i))
                    * (lagMy_) * tSUPG_
                    - dphi_dx(0,i) * (lagMy_) * tPSPG_ / dens_;
                
                rhsVector(2*i  ) += (LMx * dTime_ * timeScheme_) * 
                    weight_ * djac_;
                rhsVector(2*i+1) += (LMy * dTime_ * timeScheme_) * 
                    weight_ * djac_;
            };


     
            
            
        };
        index++;        
    };  
    
    

    lagrMultMatrix += diffMatrix;

    //lagrMultMatrix.clear();

    // ublas::identity_matrix<double>             ident(18);
    // lagrMultMatrix = ident;

    // lagrMultMatrix(12,12) = 0.;
    // lagrMultMatrix(13,13) = 0.;
    // lagrMultMatrix(14,14) = 0.;


    //Computes the Newton's method jacobian
    //jacobianNRMatrix = lagrMultMatrix;

    return;
};


//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getLMStabilizationSameMesh(){

    //  typename QuadShapeFunction<2>::Coords xsi;
    // int index = 0;

    // diffMatrix.clear();
    // jacobianNRMatrix.clear();
    // rhsVector.clear();
    // lagrMultMatrix.clear();
    
    // for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
    //     it != nQuad.end(); it++){
        
    //     //Defines the integration points adimentional coordinates
    //     xsi(0) = nQuad.PointList(index,0);
    //     xsi(1) = nQuad.PointList(index,1);

    //     //Computes the velocity shape functions
    //     shapeQuad.evaluate(xsi,phi_);

    //     //Returns the quadrature integration weight
    //     weight_ = nQuad.WeightList(index);

    //     //Computes the jacobian matrix
    //     getJacobianMatrix(xsi);

    //     getSpatialDerivatives(xsi);
        
    //     getVelAndDerivatives();     
        
        
    //     for (int i = 0; i < 6; i++){
    //         for (int j = 0; j < 6; j++){        
                
    //         //Mass matrix (use for both directions)
    //         double mM =  phi_(i) * phi_(j) * dens_ * tARLQ_;
                         
    //         //Convection matrix
    //         double Cxx = phi_(i) * ((u_ - umesh_) * dphi_dx(0,j) + 
    //                                 (v_ - vmesh_) * dphi_dx(1,j)) 
    //             * dens_ * tARLQ_;
            
    //         double Cyy  = Cxx;

            
    //         double Cuu = phi_(i) * (du_dx + dv_dx) * phi_(j) * dens_ * tARLQ_;
                        
    //         double Cvv = phi_(i) * (du_dy + dv_dy) * phi_(j) * dens_ * tARLQ_;
            
    //         // jacobianNRMatrix(2*i  ,2*j  ) += (mM + timeScheme_ * dTime_ * 
    //         //                                   (Cxx + Cuu))
    //         //     * weight_ * djac_;

    //         // jacobianNRMatrix(2*i+1,2*j+1) += (mM + timeScheme_ * dTime_ *
    //         //                                   (Cyy + Cvv))
    //         //     * weight_ * djac_;

    //         };

    //         double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    //         double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;
            
    //         double mx = phi_(i) * (u_ - uPrev_) * tARLQ_ * dens_;
    //         double my = phi_(i) * (v_ - vPrev_) * tARLQ_ * dens_;
        
    //         double Cx = phi_(i) * ((una_ - umesh_) * du_dx + 
    //                                (vna_ - vmesh_) * du_dy) * tARLQ_ *dens_ + 
    //             phi_(i) * (du_dx + dv_dx) * u_ * dens_ * tARLQ_;
    //         double Cy = phi_(i) * ((una_ - umesh_) * dv_dx + 
    //                                (vna_ - vmesh_) * dv_dy) * tARLQ_ *dens_ + 
    //             phi_(i) * (du_dy + dv_dy) * v_ * dens_ * tARLQ_;

            
    //         // rhsVector(2*i  ) += (-mx - Cx * dTime_ * timeScheme_) * 
    //         //     weight_ * djac_;
    //         // rhsVector(2*i+1) += (-my - Cy * dTime_ * timeScheme_) * 
    //         //     weight_ * djac_;
    //     };         
        
    //     index++;        
    // }; 

    return;
};


//------------------------------------------------------------------------------
//----------------------ELEMENT DIFFUSION/VISCOSITY MATRIX----------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getStabCoarse(int ielem){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    typename QuadShapeFunction<2>::Values phiLM_; 
   
    diffMatrix.clear();
    jacobianNRMatrix.clear();
    rhsVector.clear();
    lagrMultMatrix.clear();
    
    for(typename SpecialQuad::QuadratureListIt it = sQuad.begin(); 
        it != sQuad.end(); it++){
        
        if ((intPointCorrespElem[index] == ielem)){

            //Defines the integration points adimentional coordinates
            xsi(0) = sQuad.PointList(index,0);
            xsi(1) = sQuad.PointList(index,1);
            
            //Computes the velocity shape functions
            shapeQuad.evaluate(xsi,phi_);
            
            //Computes the coarse mesh shape functions
            shapeQuad.evaluate(intPointCorrespXsi[index],phiLM_);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightList(index);
            
            //Computes the jacobian matrix
            getJacobianMatrix(xsi);
                        
            getSpatialDerivatives(xsi);

            dphiL_dx = dphi_dx;

            getSpatialDerivatives(intPointCorrespXsi[index]);

            double LM = 0.;

            for (int i = 0; i < 6; i++){
                for (int j = 0; j < 6; j++){
                    LM = ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i))
                        * phiLM_(j) * tARLQ_;
            
                    jacobianNRMatrix(2*i  ,2*j  ) += timeScheme_ * dTime_ * LM
                        * weight_ * djac_ * intPointWeightFunction(index);
                    
                    jacobianNRMatrix(2*i+1,2*j+1) += timeScheme_ * dTime_ * LM
                        * weight_ * djac_ * intPointWeightFunction(index);
                    
          
                    double Lx = 0.;
                    double Ly = 0.;
                    
                    if (glueZone){
                        if (model){
                            //Fine
                            Lx = -dphi_dx(0,i) * phi_(j) * tARLQ_ / dens_;
                            Ly = -dphi_dx(1,i) * phi_(j) * tARLQ_ / dens_;
                        }else{
                            //Coarse
                            Lx = dphi_dx(0,i) * phi_(j) * tARLQ_ / dens_ * 0.;
                            Ly = dphi_dx(1,i) * phi_(j) * tARLQ_ / dens_ * 0.;
                        };
                    };
                    
                    
                    jacobianNRMatrix(12+j,2*i  ) +=  ((Lx) * 
                                                      timeScheme_ * dTime_) * 
                        weight_ * djac_;// * intPointWeightFunction(index);
                    jacobianNRMatrix(12+j,2*i+1) +=  ((Ly) * 
                                                      timeScheme_ * dTime_)
                        * weight_ * djac_;// * intPointWeightFunction(index);
                };

                // double LMx = 0.; double LMy = 0.;
                // LMx = ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i)) 
                //     * lagMx_ * tARLQ_;
                // LMx = ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i)) 
                //     * lagMy_ * tARLQ_;   
                

            };
        };
        index++;
    };
    
    return;
};










#endif

