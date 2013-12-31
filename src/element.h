#ifndef ELEMENT_H
#define	ELEMENT_H

#include <armadillo>

using namespace arma;

class MechElem {
public:
    MechElem();
    MechElem(int, double, vec*, vec*, int*, mat*, int*);
    MechElem(const MechElem &orig);
    virtual ~MechElem();
    virtual void bodyForce(vec&) =0;
    virtual void traction(vec&) =0;
    virtual void stiffness(mat&, int) =0;
    mat *getGcoords();
    int *getGdofs();
    int *getGnodes();
protected:
    int E;
    double v;
    vec *pb, *pt;
    mat *pgcoords;
    int *gnodes;
    int *gdofs;
    virtual mat N(double, double) =0;
    virtual mat J(double, double) =0;
};


/*
 * Shell elements
 * - Q4:        four node quad
 * - Q4R:       four node quad with reduced integration
 */

enum PlaneState { PSTRESS, PSTRAIN };

#define Q4__NUM_NODES 4
#define Q4__DOF_PER_NODE 2

class Q4 : public MechElem {
private:
    static const unsigned int numPointsX = 2;
    static const unsigned int numPointsY = 2;
    static const double gaussPoints[];
    static const int weights[];
    static const unsigned int totalDOF = Q4__NUM_NODES*Q4__DOF_PER_NODE;
protected:
    double h;
    mat N(double, double);
    mat J(double, double);
public:
    Q4();
    Q4(int, double, double, vec*, vec*, int*, mat*, int*);
    Q4(const Q4 &orig);
    virtual ~Q4();
    void bodyForce(vec&);
    void traction(vec&);
    void stiffness(mat&, int);
};

class Q4R : public MechElem {
private:
    static const unsigned int xi = 0;
    static const unsigned int eta = 0;
    static const int weight = 2;
protected:
    double h;
    // @todo create N1, N2, N3, N4, dN1_dxi, etc., inline functions?
    mat N(double, double);
    mat J(double, double);    
public:
    Q4R();
    Q4R(int, double, double, vec*, vec*, int*, mat*, int*);
    Q4R(const Q4 &orig);
    virtual ~Q4R();
    void bodyForce(vec&);
    void traction(vec&);
    void stiffness(mat&, int);
};

#endif	/* ELEMENT_H */