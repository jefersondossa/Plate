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

    /// Defines the local vector type with dimension 12x12 for DIM=2 and
    /// 30x30 for DIM=3
    typedef ublas::bounded_vector<double, 18*DIM-24>            LocalVector;

    /// Defines the local matrix type with dimension 12x12
    /// for DIM=2 and 30x30 for DIM=3
    typedef ublas::bounded_matrix<double, 18*DIM-24, 18*DIM-24> LocalMatrix;

    /// Defines the normal integration quadrature rule class locally
    typedef IntegQuadrature<DIM>                                NormalQuad;

    /// Defines the boundary integration quadrature rule class locally
    typedef BoundaryIntegQuadrature<DIM>                        BoundaryQuad;
    
    /// Define the type VecLocD from class Node locally
    typedef typename Nodes::VecLocD                             VecLoc;

private:
    QuadShapeFunction<DIM> shapeQuad; //Quadratic shape function
    BoundShapeFunction<DIM>shapeBound;//Boundary shape function
    NormalQuad             nQuad;     //Integration quadrature
    std::vector<Nodes *>   nodes_;    //Nodes
    Connectivity  connect_;           //Element connectivity 
    int           index_;             //Element index
    LocalNodes    localNodes_;        //Nodal coordinates
    LocalNodes    localNodesBoundary_;//Nodal coordinates - boundary
    double        djac_;              //Jacobian determinant
    double        weight_;            //Integration point weight
    DimMatrix     ainv_;             //Inverse of Jacobian transform
    double        poisson_;           //Fluid dynamic viscosity
    double        elastic_;           //Fluid dynamic viscosity
    double        thick_;             //Element thickness
    double        damping_;           //Element damping
    double        dens_;              //Fluid Density
    double        dTime_;             //Time Step
    double        timeBeta_;          //Time Integration scheme
    double        timeGamma_;          //Time Integration scheme
    LocalMatrix   stiffnessMatrix;    //Stiffness Matrix

    DimVector     fieldForce_;        //Element field force
    LocalVector   externalForces;     //Field forces integrated over
                                      //the element domain
    LocalVector   loadVector;         //RHS vector of Newton's method
    int           sideBoundary_;
    bool          loaded_;
    
    std::vector<ublas::bounded_vector<double,DIM> > di;

    //Shape functions second derivatives
    typename QuadShapeFunction<DIM>::ValueDDeriv ddphi_dx; 
    //Shape functions first derivatives
    typename QuadShapeFunction<DIM>::ValueDeriv  dphi_dx, dphiL_dx;
    //Shape functin values
    typename QuadShapeFunction<DIM>::Values      phi_;     
    //Boundary shape functions
    typename BoundShapeFunction<DIM>::Values     phib_;     
    //Boundary shape functions values
    typename BoundShapeFunction<DIM>::ValueDeriv dphib_;     

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
        dens_ = 0.;             dTime_ = 0.;              poisson_ = 0.;
        thick_ = 0.;            elastic_ = 0.;            damping_ = 0.;
        timeBeta_ = 0.;         timeGamma_ = 0.;

        fieldForce_.clear();     externalForces.clear();     loadVector.clear();
        sideBoundary_ = 0;
        loaded_ = false;

        setLocalNodes();
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

    /// Sets the element Poisson coefficient
    /// @param double element Poisson Coeficient
    void setPoissonCoefficient(double& pois){poisson_ = pois;}

    /// Sets the element Elastic modulus
    /// @param double element Elastic Modulus
    void setElasticModulus(double& elast){elastic_ = elast;}

    /// Sets the element Thickness
    /// @param double element thickness
    void setThickness(double& thick){thick_ = thick;}

    /// Sets the element damping
    /// @param double element damping
    void setDamping(double& damp){damping_ = damp;}

    /// Sets the time step size
    /// @param double time step size
    void setTimeStep(double& dt){dTime_ = dt;}

    /// Sets the time integration scheme
    /// @param double time integration scheme: 0.0 - Explicit forward Euler;
    /// 1.0 - Implicit backward Euler;
    /// 0.5 - Implicit Trapezoidal Rule.
    void setTimeIntegrationScheme(double& b, double& g){
        timeBeta_ = b; 
        timeGamma_ = g;};

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

    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian(){return djac_;};

    /// Sets the element side in boundary
    /// @param int side in boundary
    void setElemSideInBoundary(int side){sideBoundary_ = side;};

    /// Gets the element side in boundary
    /// @return side in boundary
    int getElemSideInBoundary(){return sideBoundary_;};

    /// Sets if the element side is loaded
    /// @param int side in boundary
    void setElemLoaded(){loaded_ = true;};

    /// Gets the element side in boundary
    /// @return side in boundary
    bool getElemLoaded(){return loaded_;};

    //.......................Element vectors and matrices.......................
    /// Compute and store the element stiffness matrix
    void getElemStiffMatrix();

    /// Sets the boundary conditions
    void setBoundaryConditions();

    ///Compute and store the residual vector for the incompressible flow problem
    void getElemLoadVector();

    /// Gets the element Stiffness matrix
    LocalMatrix getStiffnessMatrix(){return stiffnessMatrix;};

    /// Gets the load vector
    LocalVector getLoadVector(){return loadVector;};

    //...............................Problem type...............................
    /// Compute the Static problem matrices and vectors
    void getStatic();


};


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

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
//------------------------SEARCH THE ELEMENT LOCAL NODES------------------------
//------------------------------------------------------------------------------
//Creates a vector with the local nodes coordinates
template<>
void Element<2>::setLocalNodes() {
    
    typename Nodes::VecLocD x;

    for (int i=0; i<6; i++){
        x = nodes_[connect_(i)] -> getInitialCoordinates();
        localNodes_(i,0) = x(0);
        localNodes_(i,1) = x(1);
    };

    return;
};

template<>
void Element<3>::setLocalNodes() {
    
    typename Nodes::VecLocD x;

    for (int i=0; i<10; i++){
        x = nodes_[connect_(i)] -> getInitialCoordinates();
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
    
    externalForces.clear();   loadVector.clear();
    localNodes_.clear();     

    ddphi_dx.clear();        dphi_dx.clear();            phi_.clear();     

    loaded_ = false;
    
    setLocalNodes();

    return;
}; 

template<>
void Element<3>::clearVariables(){

    djac_ = 0.;               weight_ = 0.;              ainv_.clear();
    
    externalForces.clear();   loadVector.clear();
    localNodes_.clear();     

    ddphi_dx.clear();        dphi_dx.clear();            phi_.clear();     

    loaded_ = false;
    
    setLocalNodes();

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
//---------------------------ELEMENT STIFFNESS MATRIX---------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getElemStiffMatrix(){

    stiffnessMatrix.clear();

    ublas::bounded_matrix<double, 3, 3 > hooke;

    hooke.clear();
    
    // For EPT
    double k = elastic_ / (1. - poisson_ * poisson_);
    hooke(0,0) = k;
    hooke(0,1) = k * poisson_;
    hooke(1,0) = k * poisson_;
    hooke(1,1) = k;
    hooke(2,2) = k * (1. - poisson_) * 0.5;

    // For EPD
    // double k = elastic_ / ((1. + poisson_) * (1. - 2. * poisson_));
    // hooke(0,0) = k * (1. - poisson_);
    // hooke(0,1) = k * poisson_;
    // hooke(1,0) = k * poisson_;
    // hooke(1,1) = k * (1. - poisson);
    // hooke(2,2) = k * (1. - 2. * poisson_) * 0.5;
    
    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        ublas::bounded_matrix<double, 3,12> bMatrix, aux;
        bMatrix.clear(); aux.clear();
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);
        
        //Computes the shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);
        
        //Computes spatial derivatives
        getSpatialDerivatives(xsi);
        
        for (int i = 0; i < 6; i++){                
            bMatrix(0,2*i  ) = dphi_dx(0,i);
            bMatrix(2,2*i  ) = dphi_dx(1,i);

            bMatrix(1,2*i+1) = dphi_dx(1,i);
            bMatrix(2,2*i+1) = dphi_dx(0,i);
        };
        
        aux = prod(hooke,bMatrix);        
        
        stiffnessMatrix += prod(trans(bMatrix),aux) * djac_ * weight_ * thick_;
        
        index++;
    };

    return;
};


//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setBoundaryConditions(){

    for (int i = 0; i < 6; i++){

        if (nodes_[connect_(i)] -> getConstrains(0) == 1)  {
            for (int j = 0; j < 12; j++){
                stiffnessMatrix(2*i  ,j) = 0.;
                stiffnessMatrix(j,2*i  ) = 0.;
            };
            stiffnessMatrix(2*i  , 2*i  ) = 1.;
            loadVector(2*i  ) = 0.;
        };

        if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
            for (int j = 0; j < 12; j++){
                stiffnessMatrix(2*i+1,j) = 0.;
                stiffnessMatrix(j,2*i+1) = 0.;
            };
            stiffnessMatrix(2*i+1, 2*i+1) = 1.;
            loadVector(2*i+1) =  0.;
        };
    };

    return;
};

//------------------------------------------------------------------------------
//----------------------------------LOAD VECTOR---------------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getElemLoadVector(){
    
    loadVector.clear();

    // Domain Loads
    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);
        
        //Computes the shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);
                
        for (int i = 0; i < 6; i++){                
            loadVector(2*i  ) += fieldForce_(0) * phi_(i) 
                * djac_ * weight_ * thick_;
            loadVector(2*i+1) += fieldForce_(1) * phi_(i) 
                * djac_ * weight_ * thick_;
        };
        
        index++;
    };


    // Surface loads
    if (loaded_) {
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
        
        n_vector.clear();
        gaussQuad.first.clear();
        gaussQuad.second.clear();
        
        gaussQuad = bQuad.GaussQuadrature();
        
        int index = 0;
        for(typename BoundaryQuad::QuadratureListIt it = bQuad.begin(); 
            it != bQuad.end(); it++){
            
            phib_ = shapeBound.getShapeFunction(gaussQuad.first(index));
            dphib_ = shapeBound.
                getShapeFunctionDerivative(gaussQuad.first(index));
            
            double Tx=0.; double Ty = 0.;
            
            for (int i=0; i<3; i++){
                Tx += localNodesBoundary_(i,0) * dphib_(i);
                Ty += localNodesBoundary_(i,1) * dphib_(i);
            };
            
            double jacb_ = sqrt(Tx*Tx + Ty*Ty);
            
            n_vector(0) =  Ty / jacb_;
            n_vector(1) = -Tx / jacb_;
            
            for (int i=0; i<3; i++){
                double qx = nodes_[nodesb_(i)] -> getConstrainValue(0);
                double qy = nodes_[nodesb_(i)] -> getConstrainValue(1);

                if(sideBoundary_ == 0){
                    loadVector(2*1  ) += qx * phib_(i) * weight_ * jacb_;
                    loadVector(2*1+1) += qy * phib_(i) * weight_ * jacb_;
                    
                    loadVector(2*4  ) += qx * phib_(i) * weight_ * jacb_;
                    loadVector(2*4+1) += qy * phib_(i) * weight_ * jacb_; 
                    
                    loadVector(2*2  ) += qx * phib_(i) * weight_ * jacb_;
                    loadVector(2*2+1) += qy * phib_(i) * weight_ * jacb_;
                }else{
                    if(sideBoundary_ == 1){
                        loadVector(2*2  ) += qx * phib_(i) * weight_ * jacb_;
                        loadVector(2*2+1) += qy * phib_(i) * weight_ * jacb_;
                        
                        loadVector(2*5  ) += qx * phib_(i) * weight_ * jacb_;
                        loadVector(2*5+1) += qy * phib_(i) * weight_ * jacb_; 
                        
                        loadVector(2*0  ) += qx * phib_(i) * weight_ * jacb_;
                        loadVector(2*0+1) += qy * phib_(i) * weight_ * jacb_;
                    }else{
                        loadVector(2*0  ) += qx * phib_(i) * weight_ * jacb_;
                        loadVector(2*0+1) += qy * phib_(i) * weight_ * jacb_;
                        
                        loadVector(2*3  ) += qx * phib_(i) * weight_ * jacb_;
                        loadVector(2*3+1) += qy * phib_(i) * weight_ * jacb_; 
                        
                        loadVector(2*1  ) += qx * phib_(i) * weight_ * jacb_;
                        loadVector(2*1+1) += qy * phib_(i) * weight_ * jacb_;
                    };        
                };
            };
        };        
            
        index++;
        
    };

    return;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-----------------------MOUNTS EACH TYPE OF PLATE PROBEM-----------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------STATIC PLATE PROBEM------------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getStatic(){
    
    //Computes the element diffusion/viscosity matrix
    getElemStiffMatrix();
    
    //Computes the RHS vector
    getElemLoadVector();
    
    //Apply boundary conditions
    setBoundaryConditions();

    return;
};

#endif

