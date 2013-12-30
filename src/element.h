#ifndef ELEMENT_H
#define	ELEMENT_H

#include <armadillo>

using namespace arma;

class MechElem {
public:
    MechElem();
    MechElem(int, float, vec*, vec*, int*, mat*, int*);
    MechElem(const MechElem& orig);
    virtual ~MechElem();
    virtual vec bodyForce(void) =0;
    virtual vec traction(void) =0;
    virtual mat stiffness(int) =0;
    virtual void bodyForce(vec*) =0;
    virtual void traction(vec*) =0;
    virtual void stiffness(mat*, int) =0;
//protected:
    // @todo rewrite this with accessors and mutators
    int E, v;
    vec *pb, *pt;
    mat *pgcoords;
    int *gnodes;
    int *gdof;
    virtual mat N(double, double) =0;
};


/*
 * Shell elements
 * - Q4:        four node quad
 * - Q4R:       four node quad with reduced integration
 */

enum PlaneState { PSTRESS, PSTRAIN };

#define Q4__NUM_NODES 4
#define Q4__DOF_PER_NODE 2

// @todo rewrite to better encapsulate and protect member data
class Q4 : public MechElem {
private:
    static unsigned int const numPointsX = 2;
    static unsigned int const numPointsY = 2;
    static double const gaussPoints[];
    static int const weights[];
    static unsigned int const totalDOF = Q4__NUM_NODES*Q4__DOF_PER_NODE;
public:
    Q4();
    Q4(int, double, double, vec*, vec*, int*, mat*, int*);
    int E;
    double v, h;
    vec *pb, *pt;
    mat *pgcoords;
    int *gnodes;
    int *gdof;
    mat N(double, double);
    mat J(double, double);
    vec bodyForce(void);
    vec traction(void);
    mat stiffness(int);
    void bodyForce(vec*);
    void traction(vec*);
    void stiffness(mat*, int);
};

class Q4R : public Q4 {
private:
    static unsigned int const xi = 0;
    static unsigned int const eta = 0;
    static int const weight = 2;    
public:
    Q4R(int, float, vec*, vec*, int*, mat, int*);
    int E, v;
    vec *pb, *pt;
    mat *pgcoords;
    int gnodes[Q4__NUM_NODES];
    int gdof[Q4__NUM_NODES*Q4__DOF_PER_NODE];
    mat N(double, double);
    mat J(double, double);
    vec bodyForce(void);
    vec traction(void);
    mat stiffness(bool);
};

#endif	/* ELEMENT_H */