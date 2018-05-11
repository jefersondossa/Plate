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

private:
    //Main variables
    VecLocD          coord_;                  //Nodal coordinate vector
    VecLocD          initialCoord_;           //Initial nodal coordinate vector
    int              index_;                  //Node index

    int              constrainType[3];        //Constrain direction
    double           constrainValue[3];       //Nodal prescribed 
                                              // load or displacement
    
public:
    ///Constructor - Defines a node with index and coordinates
    Node(VecLocD& coor, int index){
        coord_ = coor;
        initialCoord_ = coor;
        index_ = index; 

        constrainType[0] = 0;    constrainType[1] = 0;    constrainType[2] = 0;
        constrainValue[0] = 0;   constrainValue[1] = 0;   constrainValue[2] = 0;
    };

    /// Clear all node object variables
    void clearVariables();

    /// Returns the node coordinate vector
    /// @return node coordinate vector
    VecLocD getCoordinates() {return coord_;};

    /// Returns the node initial coordinate vector
    /// @return node initial coordinate vector
    VecLocD getInitialCoordinates() {return initialCoord_;};

    /// Increment the coordinate vector
    /// @param int direction @param double increment value
    void incrementCoordinate(int dir, double u);

    /// Sets the node coordinate vector
    /// @param VecLocD Coordinate
    void setCoordinates(VecLocD& coor){coord_ = coor;};
        
    //...........................Constrains functions...........................
    /// Sets all node constrains     
    /// @param int direction 
    /// @param int type: 0 - free, 1 - constrained, 2 - glue zone, 
    /// 3 - fluid-structure interface @param double constrain value
    void setConstrains(int dir, int type, double value){
        constrainType[dir] = type;
        constrainValue[dir] = value;
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
    coord_ = initialCoord_;
    
    return;
};

template<>
void Node<3>::clearVariables(){
    coord_ = initialCoord_;
    
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

#endif
