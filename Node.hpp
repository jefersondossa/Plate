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
    /// which helps to setting the constrains
    typedef std::vector<std::pair<int, double> >     VecConstrD;

    /// Defines the node iterator for helping to build loops in order
    /// to enforce bcs
    typedef VecConstrD::iterator                     VecConstrDIt;

    static const int spaceDim = DIM;

private:
    //Main variables
    VecLocD          coord_;                  //Nodal coordinate vector
    VecLocD          initialCoord_;           //Initial nodal coordinate vector
    int              index_;                  //Node index
    VecLocD          previousCoord_;          //Previous nodal coordinate vector

    int              constrainType[3];        //Constrain direction
    double           constrainValue[3];       //Nodal prescribed velocity

    int              constrainTypeLaplace[3]; //Constrain direction Laplace
    double           constrainValueLaplace[3];//Nodal prescribed value

    VecLocD          velocity_;               //Nodal velocity
    VecLocD          previousVelocity_;       //Previous time step velocity

    VecLocD          acceleration_;           //Nodal acceleration
    VecLocD          previousAcceleration_;   //Previous time step acceleration
    
public:
    ///Constructor - Defines a node with index and coordinates
    Node(VecLocD& coor, int index){
        coord_ = coor;
        initialCoord_ = coor;
        previousCoord_ = coor;
        index_ = index; 

        constrainType[0] = 0;    constrainType[1] = 0;    constrainType[2] = 0;
        constrainValue[0] = 0;   constrainValue[1] = 0;   constrainValue[2] = 0;
  
        velocity_.clear();   previousVelocity_.clear();   acceleration_.clear();
        previousAcceleration_.clear(); 
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

    /// Increment the coordinate vector
    /// @param int direction @param double increment value
    void incrementCoordinate(int dir, double u);

    /// Sets the previous coordinate vector
    /// @param int direction @param double value
    void setPreviousCoordinates(int dir, double u);

    /// Sets the node coordinate vector
    /// @param VecLocD Coordinate
    void setCoordinates(VecLocD& coor){coord_ = coor;};
    
    
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

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------------CLEAR VARIABLES--------------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::clearVariables(){
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
   
    return;
};

template<>
void Node<3>::clearVariables(){
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

#endif
