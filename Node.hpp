//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------NODES-------------------------------------
//------------------------------------------------------------------------------

#ifndef NODE_H
#define NODE_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric;

/// Defines the node object and stores all nodal variables information

template<int DIM>

class Node{
public:
    /// Defines a type to store double variables with the dimension size
    typedef ublas::bounded_vector<double,DIM>        VecLocD;

    /// Defines a pair with integer and  double variables
    /// which helps to setting the constraints
    typedef std::vector<std::pair<int, double> >     VecConstrD;

    /// Defines the node iterator for helping to build loops in order
    /// to enforce bcs
    typedef VecConstrD::iterator                     VecConstrDIt;

    static const int spaceDim = DIM;

private:
    //Main variables
    VecLocD          coord_;                  //Nodal coordinate vector
    VecLocD          initialCoord_;           //Initial nodal coordinate vector
    VecLocD          coordUpdated_;           //Updated nodal coordinate vector
    int              index_;                  //Node index
    VecLocD          previousCoord_;          //Previous nodal coordinate vector


    //Fluid
    int              constrainType[3];        //Constrain direction
    double           constrainValue[3];       //Nodal prescribed velocity

    int              constrainTypeLaplace[3]; //Constrain direction Laplace
    double           constrainValueLaplace[3];//Nodal prescribed value

    VecLocD          velocity_;               //Nodal velocity
    VecLocD          previousVelocity_;       //Previous time step velocity

    VecLocD          acceleration_;           //Nodal acceleration
    VecLocD          previousAcceleration_;   //Previous time step acceleration

    double           pressure_;               //Nodal pressure 
    double           previousPressure_;       //Previous time step pressure

    double           divergent_;              //Velocity divergent

    VecLocD          meshVelocity_;           //Nodal mesh velocity
    VecLocD          previousMeshVelocity_;   //Previous time step mesh velocity

    //Arlequin
    int              elemCorresp;             //Element correspondence
    VecLocD          xsiCorresp;              //Adimensional coord correspond

    double           presArlequin_;           //Glue zone pressure    
    VecLocD          velArlequin_;            //Glue zone velocity

    VecLocD          lagMultiplier_;          //Nodal Lagrange Multiplier value
    double           weightFunction_;         //Nodal Energy Weight Function
    double           distGlueZone;            //Signaled distance to glue zone

    //Potential problem
    VecLocD          gradient_;               //Potential gradient
    double           potential_;              //Potential value
    
public:
    ///Constructor - Defines a node with index and coordinates
    Node(VecLocD& coor, int index){
        coord_ = coor;
        initialCoord_ = coor;
        coordUpdated_ = coor;
        previousCoord_ = coor;
        index_ = index; 

        constrainType[0] = 0;    constrainType[1] = 0;    constrainType[2] = 0;
        constrainValue[0] = 0;   constrainValue[1] = 0;   constrainValue[2] = 0;
        pressure_ = 0.;          previousPressure_ = 0.;  divergent_ = 0.;
        elemCorresp = 0;         
        velocity_.clear();   previousVelocity_.clear();   acceleration_.clear();
        previousAcceleration_.clear(); xsiCorresp.clear(); gradient_.clear();
        potential_ = 0.;     lagMultiplier_.clear();      weightFunction_ = 0.;
        distGlueZone = 0.;

        constrainTypeLaplace[0] = 0;    constrainTypeLaplace[1] = 0;
        constrainTypeLaplace[2] = 0;    constrainValueLaplace[0] = 0;
        constrainValueLaplace[1] = 0;   constrainValueLaplace[2] = 0;

    };

    /// Clear all node object variables
    void clearVariables();

    /// Returns the node coordinate vector
    /// @return node coordinate vector
    VecLocD getCoordinates() {return coord_;};

    /// Returns the node initial coordinate vector
    /// @return node initial coordinate vector
    VecLocD getInitialCoordinates() {return initialCoord_;};

    /// Returns the node coordinate vector at the previous time step
    /// @return node coordinate vector at the previous time step
    VecLocD getPreviousCoordinates() {return previousCoord_;}

    /// Returns the node updated coordinate vector
    /// @return node coordinate updated vector
    VecLocD getUpdatedCoordinates() {return coordUpdated_;};

    /// Increment the coordinate vector
    /// @param int direction @param double increment value
    void incrementCoordinate(int dir, double u);

    /// Sets the previous coordinate vector
    /// @param int direction @param double value
    void setPreviousCoordinates(int dir, double u);

    /// Sets the node coordinate vector
    /// @param VecLocD Coordinate
    void setCoordinates(VecLocD& coor){coord_ = coor;};

    /// Sets the updated coordinate vector
    /// @param VecLocD Updated Coordinate
    void setUpdatedCoordinates(VecLocD& coor){coordUpdated_ = coor;};

    /// Updates node coordinate vector
    /// @param int direction @param double updated value
    void updateCoordinate(int dir, double val){coord_(dir) = val;};
    
    /// Sets nodal correspondence of overlapped mesh
    /// @param double element @param VecLocD adimensional coordinates
    void setNodalCorrespondence(double elem, const VecLocD& xsi){
        elemCorresp = elem;
        xsiCorresp = xsi;
    };

    /// Gets nodal correspondence of overlapped mesh - element
    /// @return Element correspondence of overlapped mesh
    int getNodalElemCorrespondence() {return elemCorresp;}

    /// Gets nodal correspondence of overlapped mesh - adim. coordinate
    /// @return Adim. coordinate correspondence of overlapped mesh
    VecLocD getNodalXsiCorrespondence() {return xsiCorresp;}
    
    //............................Velocity functions............................
    /// Sets the velocity vector
    /// @param double* velocity vector
    void setVelocity(double *u);

    /// Sets the previous velocity vector
    /// @param double* previous time step velocity vector
    void setPreviousVelocity(double *u);

    /// Increment the velocity vector
    /// @param int direction @param double increment value
    void incrementVelocity(int dir, double u);

    /// Returns the node velocity vector
    /// @return node velocity vector
    double getVelocity(int dir) {return velocity_(dir);}

    /// Returns the node previous time step velocity vector
    /// @return node previous time step velocity vector
    double getPreviousVelocity(int dir) {return previousVelocity_(dir);}

    /// Sets the velocity divergent at the node
    /// @param double velocity divergent
    void setVelocityDivergent(double div) {divergent_ = div;}

    /// Returns the node velocity divergent
    /// @return node velocity divergent
    double getVelocityDivergent() {return divergent_;}

    //..........................Acceleration functions..........................
    /// Sets the acceleration vector
    /// @param double* acceleration vector
    void setAcceleration(double *u);

    /// Sets the previous time step acceleration vector
    /// @param double* previous time step acceleration vector
    void setPreviousAcceleration(double *u);

    /// Gets the acceleration vector
    /// @return acceleration vector
    double getAcceleration(int dir) {return acceleration_(dir);}

    /// Gets the previous time step acceleration vector
    /// @return previous time step acceleration vector
    double getPreviousAcceleration(int dir) {return previousAcceleration_(dir);}

    //............................Pressure functions............................
    /// Sets the nodal pressure
    /// @param double pressure
    void setPressure(double p);

    /// Increments the nodal pressure
    /// @param double pressure increment
    void incrementPressure(double p);

    /// Gets the nodal pressure value
    /// @return nodal pressure value
    double getPressure() {return pressure_;};

    /// Sets the previous time step nodal pressure
    /// @param double previous time step nodal pressure
    void setPreviousPressure(double p);

    /// Gets the previous time step nodal pressure
    /// @return previous time step nodal pressure
    double getPreviousPressure() {return previousPressure_;}

    //.........................Mesh Velocity functions..........................
    /// Sets the node mesh velocity
    /// @param double* mesh velocity
    void setMeshVelocity(double *u);

    /// Sets the previous time step mesh velocity
    /// @param int direction @param double previous time step mesh velocity valu
    void setPreviousMeshVelocity(int dir, double u);

    /// Gets the node mesh velocity
    /// @param int direction @return mesh velocity component
    double getMeshVelocity(int dir) {return meshVelocity_(dir);}

    /// Gets the previous time step mesh velocity
    /// @param int direction @return previous time step mesh velocity component
    double getPreviousMeshVelocity(int dir) {return previousMeshVelocity_(dir);}

    //...........................Constrains functions...........................

    /// Sets all node constrains     
    /// @param int direction 
    /// @param int type: 0 - free, 1 - constrained, 2 - glue zone, 
    /// 3 - fluid-structure interface @param double constrain value
    void setConstrains(int dir, int type, double value){
        constrainType[dir] = type;
        constrainValue[dir] = value;
        velocity_(dir) = value;
        previousVelocity_(dir) = value;
    };

    /// Gets node constrain type
    /// @return constrain type
    int getConstrains(int dir) {return constrainType[dir];}

    /// Gets node constrain value
    /// @return constrain value
    double getConstrainValue(int dir) {return constrainValue[dir];}

    /// Sets constrains for solving the mesh moving problem
    /// @param int direction 
    /// @param int constrain type: 0 - free, 1 - constrained
    /// @param double constrain value
    void setConstrainsLaplace(int dir, int type, double value){
        constrainTypeLaplace[dir] = type;
        constrainValueLaplace[dir] = value;
    };

    /// Gets constrains of mesh moving problem
    /// @return constrain type
    int getConstrainsLaplace(int dir) {return constrainTypeLaplace[dir];}

    /// Gets constrain value of mesh moving problem
    /// return constrain value
    double getConstrainValueLaplace(int dir){return constrainValueLaplace[dir];}

    //............................Arlequin functions............................
    /// Sets Lagrange Multiplier value
    /// @param int direction @param double component value
    void setLagrangeMultiplier(int dir, double lMult){
        lagMultiplier_(dir) = lMult;};

    /// Increment the Lagrange Multiplier vector
    /// @param int direction @param double increment value
    void incrementLagrangeMultiplier(int dir, double lMult){
        lagMultiplier_(dir) += lMult;};

    /// Gets Lagrange Multiplier component value
    /// @param int direction @return component value
    double getLagrangeMultiplier(int dir) {return lagMultiplier_(dir);};

    /// Sets the nodal energy weight function value
    /// @param double weight function value
    void setWeightFunction(double val) {weightFunction_ = val;};

    /// Gets the nodal energy weight function value
    /// @return weight function value
    double getWeightFunction() {return weightFunction_;};

    /// Sets the interpolated Arlequin pressure 
    /// @param double interpolated pressure value
    void setPressureArlequin(double p) {presArlequin_ = p;};

    /// Gets interpolated Arlequin pressure
    /// @return interpolated Arlequin pressure
    double getPressureArlequin() {return presArlequin_;};

    /// Sets the interpolated Arlequin velocity
    /// @param int direction @param double interpolated velocity component value
    void setVelocityArlequin(int dir, double v) {velArlequin_(dir) = v;};

    /// Gets the interpolated Arlequin velocity component
    /// @param int direction @return interpolated Arlequin velocity component
    double getVelocityArlequin(int dir) {return velArlequin_(dir);};

    //.......................Signaled distance functions........................
    /// Sets the Signaled distance function value
    /// @param double signaled distance value
    void setDistFunction(double dist) {distGlueZone = dist;};

    /// Gets the Signaled distance function value
    /// @return signaled distance value
    double getDistFunction(){return distGlueZone;};

    //.......................Potential problem functions........................
    /// Sets the potential gradient component value
    /// @param int direction @param double potential gradient component value
    void setGradientComponent(double val, int dir) {gradient_(dir) = val;};

    /// Gets the potential gradient component value
    /// @param int direction @return potential gradient component value
    double getGradientComponent(int dir) {return gradient_(dir);};

    /// Sets the nodal potential value
    /// @param double potential 
    void setPotential(double p) {potential_ = p;};

    /// Increments the nodal potential value
    /// @param double potential increment value
    void incrementPotential(double p) {potential_ += p;};

    /// Gets the nodal potential value
    /// @return potential value
    double getPotential() {return potential_;};

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------------CLEAR VARIABLES--------------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::clearVariables(){
    pressure_ = 0.;          previousPressure_ = 0.;  divergent_ = 0.;
    elemCorresp = 0;

    if (constrainType[0] != 1){
        velocity_(0) = 0.;   
        previousVelocity_(0) = 0.;   
        acceleration_(0) = 0.;
        previousAcceleration_(0) = 0.;
    };

    if (constrainType[1] != 1){
        velocity_(1) = 0.;   
        previousVelocity_(1) = 0.;   
        acceleration_(1) = 0.;
        previousAcceleration_(1) = 0.;
    };
    
    xsiCorresp.clear();    gradient_.clear();
    potential_ = 0.;       lagMultiplier_.clear();      weightFunction_ = 0.;
    
    return;
};

template<>
void Node<3>::clearVariables(){

    pressure_ = 0.;          previousPressure_ = 0.;  divergent_ = 0.;
    elemCorresp = 0;

    if (constrainType[0] != 1){
        velocity_(0) = 0.;   
        previousVelocity_(0) = 0.;   
        acceleration_(0) = 0.;
        previousAcceleration_(0) = 0.;
    };

    if (constrainType[1] != 1){
        velocity_(1) = 0.;   
        previousVelocity_(1) = 0.;   
        acceleration_(1) = 0.;
        previousAcceleration_(1) = 0.;
    };

    if (constrainType[2] != 1){
        velocity_(2) = 0.;   
        previousVelocity_(2) = 0.;   
        acceleration_(2) = 0.;
        previousAcceleration_(2) = 0.;
    };
    
    xsiCorresp.clear();    gradient_.clear();
    potential_ = 0.;       lagMultiplier_.clear();     weightFunction_ = 0.;
    
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setVelocity(double *u){
    //Sets Velocity value
    velocity_(0) = u[0];
    velocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setVelocity(double *u){
    //Sets Velocity value
    velocity_(0) = u[0];
    velocity_(1) = u[1]; 
    velocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//-----------------------INCREMENT NODAL VELOCITY VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementVelocity(int dir, double u){
    //All element nodes
    velocity_(dir) += u;
    return;
};

template<>
void Node<3>::incrementVelocity(int dir, double u){
    //All element nodes
    velocity_(dir) += u;

    return;
};

//------------------------------------------------------------------------------
//---------------------INCREMENT NODAL COORDINATE VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementCoordinate(int dir, double u){
    //All element nodes
    coord_(dir) += u;
    return;
};

template<>
void Node<3>::incrementCoordinate(int dir, double u){
    //All element nodes
    coord_(dir) += u;

    return;
};

//------------------------------------------------------------------------------
//---------------------SETS PREVIOUS NODAL VELOCITY VALUES----------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousCoordinates(int dir, double u){
    //All element nodes
    previousCoord_(dir) = u;
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL PRESSURE VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPressure(double p){
    //Only element nodes 0, 1 and 2  
    pressure_ = p; 
    return;
};

template<>
void Node<3>::setPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    pressure_ = p; 
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL PRESSURE VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousPressure(double p){
    //Only element nodes 0, 1 and 2  
    previousPressure_ = p; 
    return;
};

template<>
void Node<3>::setPreviousPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    previousPressure_ = p; 
    return;
};

//------------------------------------------------------------------------------
//-----------------------INCREMENT NODAL PRESSURE VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementPressure(double p){
    //Only element nodes 0, 1 and 2  
    pressure_ += p; 
    return;
};

template<>
void Node<3>::incrementPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    pressure_ += p; 
    return;
};

//------------------------------------------------------------------------------
//---------------------SETS PREVIOUS NODAL VELOCITY VALUES----------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousVelocity(double *u){
    //All element nodes
    previousVelocity_(0) = u[0];
    previousVelocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setPreviousVelocity(double *u){
    //All element nodes
    previousVelocity_(0) = u[0];
    previousVelocity_(1) = u[1]; 
    previousVelocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//--------------------------NODAL ACCELERATION VALUES---------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setAcceleration(double *u){
    //All element nodes
    acceleration_(0) = u[0];
    acceleration_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setAcceleration(double *u){
    //All element nodes
    acceleration_(0) = u[0];
    acceleration_(1) = u[1]; 
    acceleration_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//-------------------SETS PREVIOUS NODAL ACCELERATION VALUES--------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousAcceleration(double *u){
    //All element nodes
    previousAcceleration_(0) = u[0];
    previousAcceleration_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setPreviousAcceleration(double *u){
    //All element nodes
    previousAcceleration_(0) = u[0];
    previousAcceleration_(1) = u[1]; 
    previousAcceleration_(2) = u[2]; 
    return;
};


//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setMeshVelocity(double *u){
    //All element nodes
    previousMeshVelocity_(0) = meshVelocity_(0);
    previousMeshVelocity_(1) = meshVelocity_(1); 

    meshVelocity_(0) = u[0];
    meshVelocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setMeshVelocity(double *u){
    //All element nodes
    meshVelocity_(0) = u[0];
    meshVelocity_(1) = u[1]; 
    meshVelocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousMeshVelocity(int dir, double u){
    //All element nodes
    previousMeshVelocity_(dir) = u;

    return;
};


#endif
