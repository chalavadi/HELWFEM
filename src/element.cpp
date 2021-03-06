#include "element.h"
#include <armadillo>

//#define NDEBUG 1
#include <cassert>

using namespace arma;

// TODO let's store material properties in a struct or an instance and only keep a reference to it
// TODO most matrices we can make a fixed, static size. implement this for better performance
// TODO consider a more flyweight approach... how much can we slim down cost of creating an element?
// TODO convert c style arrays to c++ arrays
// TODO separate material data from mesh data
// TODO convert for loops to ranged loops
// TODO delegate integration to a separate "Integration" class

/**
 * Default constructor for a mechanical finite element
 * 
 * @return MechElem
 */
MechElem::MechElem() {}

/**
 * Constructor for a mechanical finite element
 * 
 * @param [int] E Modulus of elasticity
 * @param [double] v Poissons ratio
 * @param [vec*] b Body forces: x-direction, y-direction
 * @param [vec*] t Traction forces: bottom, right, top, left
 * @param [int*] gnodes Nodal numbers: bot-left, bot-right, top-right, top-left
 * @param [mat*] gcoords Global coordinates corresponding to each node
 * @param [int*] gdof Global degrees of freedom
 * @return MechElem
 */
MechElem::MechElem(long long int E, double v, vec *pb, vec *pt, int *gnodes, 
        mat *pgcoords, int *gdofs) : E(E), v(v), pb(pb), pt(pt), 
        gnodes(gnodes), pgcoords(pgcoords), gdofs(gdofs) {};

/**
 * Clone of a mechanical finite element
 * 
 * @param [int] 
 * @return MechElem
 */
MechElem::MechElem(const MechElem& orig) 
{
}

/**
 * Deallocation of a mechanical finite element
 * 
 * @param void 
 * @return void
 */
MechElem::~MechElem(void) 
{
}

/**
 * Accessor for global nodal coordinates
 * 
 * @return mat *gcoords
 */
mat *MechElem::getGcoords() { return pgcoords; }

/**
 * Accessor for global node numbers
 * 
 * @return int *gnodes
 */
int *MechElem::getGnodes() { return gnodes; }

/**
 * Accessor for global degrees of freedom
 * 
 * @return int *gdofs
 */
int *MechElem::getGdofs() { return gdofs; }

/*
 * Q4 class constants
 */
double const Q4::gaussPoints[] = { -0.5774, 0.5774 };
int const Q4::weights[] = { 1, 1 };

/**
 * Default constructor for a Q4 finite element
 * 
 * @return Q4
 */
Q4::Q4() {}

// TODO: consider separating material data from element object

/**
 * Constructor for a Q4 finite element
 * 
 * @param [int] E Modulus of elasticity
 * @param [double] v Poissons ratio
 * @param [double] h Thickness or 'height' of the element
 * @param [vec] b Body forces: x-direction, y-direction
 * @param [vec] t Traction forces: bottom, right, top, left
 * @param [int*] gnodes Nodal numbers: bot-left, bot-right, top-right, top-left
 * @param [mat*] gcoords Global coordinates corresponding to each node
 * @param [int*] gdof Global degrees of freedom
 * @return Q4
 */
Q4::Q4(long long int E, double v, double h, vec *pb, vec *pt, int *gnodes, 
    mat *pgcoords, int *gdofs) : MechElem(E, v, pb, pt, gnodes, pgcoords, 
    gdofs), h(h) {};


/**
 * Deallocation of a Q4 finite element
 * 
 * @param void 
 * @return void
 */
Q4::~Q4(void) 
{
}

// consider inline functions?
/*
 * Define the Q4 shape functions
 */
#define Q4__N1(xi, eta) (1-xi)*(1-eta)/4
#define Q4__N2(xi, eta) (1+xi)*(1-eta)/4
#define Q4__N3(xi, eta) (1+xi)*(1+eta)/4
#define Q4__N4(xi, eta) (1-xi)*(1+eta)/4
        
/**
 * Calculates the N matrix
 * 
 * @param [double] xi Xi coordinate in parent coordinates
 * @param [double] eta Eta coordinate in parent coordinates
 * @return [mat] N matrix
 */
mat Q4::N(double xi, double eta) 
{
    mat::fixed<Q4__DOF_PER_NODE,Q4__DOF_PER_NODE*Q4__NUM_NODES> N;
    N << Q4__N1(xi, eta) << 0 << Q4__N2(xi, eta) << 0 << Q4__N3(xi, eta) << 0
            << Q4__N4(xi, eta) << 0 << endr
      << 0 << Q4__N1(xi, eta) << 0 << Q4__N2(xi, eta) << 0 << Q4__N3(xi, eta)
            << 0 << Q4__N4(xi, eta) << endr;
    return N;
}

/*
 * Define the partials of the Q4 shape functions with respect to xi and eta
 */
#define Q4__dN1_dxi(eta) -(1-eta)/4
#define Q4__dN2_dxi(eta) (1-eta)/4
#define Q4__dN3_dxi(eta) (1+eta)/4
#define Q4__dN4_dxi(eta) -(1+eta)/4
#define Q4__dN1_deta(xi) -(1-xi)/4
#define Q4__dN2_deta(xi) -(1+xi)/4
#define Q4__dN3_deta(xi) (1+xi)/4
#define Q4__dN4_deta(xi) (1-xi)/4

/**
 * Calculates the jacobian
 * 
 * @param [double] xi Xi coordinate in parent coordinates
 * @param [double] eta Eta coordinate in parent coordinates
 * @return [mat] Jacobian
 */
mat Q4::J(double xi, double eta) 
{
    mat::fixed<Q4__DOF_PER_NODE,Q4__NUM_NODES> dj;
    dj << Q4__dN1_dxi(eta) << Q4__dN2_dxi(eta) << Q4__dN3_dxi(eta)
            << Q4__dN4_dxi(eta) << endr
       << Q4__dN1_deta(xi) << Q4__dN2_deta(xi) << Q4__dN3_deta(xi)
            << Q4__dN4_deta(xi) << endr;
    mat jacobian = dj * *pgcoords;
    return jacobian;
}

// TODO: have Q4 inherit from 'shell element' in order to maximize code reuse

/**
 * Calculate the element body force, mutator
 * 
 * @param [vec&] bodyForce Body force vector
 * @return void
 */
void Q4::mbodyForce(vec &bodyForce)
{
    assert(bodyForce.n_rows == Q4__DOF_PER_NODE*Q4__NUM_NODES);
    // TODO make zeros at or not at the beg of function consistent across functions
    bodyForce.zeros();
    unsigned int i, j;
    double xi, eta, weightX, weightY;
    
    /*
     * calculate element body force using gauss quadrature
     */
    for (i = 0; i < Q4::numPointsX; i++) 
    {
        xi = Q4::gaussPoints[i];
        weightX = Q4::weights[i];
        for (j = 0; j < Q4::numPointsY; j++) 
        {
            eta = Q4::gaussPoints[j];
            weightY = Q4::weights[j];
            bodyForce += weightX * weightY * trans(N(xi, eta)) * *pb 
                    * det(J(xi, eta));
        }
    }
    bodyForce *= h;
}

//TODO: fix this function
/**
 * Calculates the traction for a surface, mutator
 * 
 * @param [vec&] traction Traction force vector
 * @return void
 */
void Q4::mtraction(vec &traction) 
{
    traction.zeros();
    unsigned int i, j;
    double xi, eta, weightX, weightY;
    
    /*
     * calculate element body force using gauss quadrature
     */
    for (i = 0; i < Q4::numPointsX; i++) 
    {
        xi = Q4::gaussPoints[i];
        weightX = Q4::weights[i];
        for (j = 0; j < Q4::numPointsY; j++) 
        {
            eta = Q4::gaussPoints[j];
            weightY = Q4::gaussPoints[j];
            traction += weightX * weightY * trans(N(xi, eta)) * *pt 
                    * det(J(xi, eta));
        }
    }
    traction *= h;
}

#define Q4__STRAIN_COMP 3

/**
 * Calculates the stiffness matrix for the element, mutator
 * 
 * @param [mat*] stiff Element stiffness matrix
 * @param [PlaneState] pState Plane state control
 * @return void
 */
void Q4::mstiffness(mat &stiff, int pState = PSTRESS) 
{
    assert(stiff.n_rows == Q4__DOF_PER_NODE*Q4__NUM_NODES);
    assert(stiff.n_cols == Q4__DOF_PER_NODE*Q4__NUM_NODES);
    
    stiff.zeros();
    Mat<int>::fixed<3,4> e;
    e << 1 << 0 << 0 << 0 << endr
      << 0 << 0 << 0 << 1 << endr
      << 0 << 1 << 1 << 0 << endr;
    unsigned int i, j, k, m;
    double xi, eta, weightX, weightY, temp;
    mat nStar = zeros<mat>(Q4__DOF_PER_NODE*Q4__DOF_PER_NODE,Q4__DOF_PER_NODE*
                Q4__NUM_NODES);
    mat jacobian(Q4__DOF_PER_NODE*Q4__DOF_PER_NODE,
            Q4__DOF_PER_NODE*Q4__DOF_PER_NODE);
    mat invJac(Q4__DOF_PER_NODE*Q4__DOF_PER_NODE,
            Q4__DOF_PER_NODE*Q4__DOF_PER_NODE);
    mat JE = zeros<mat>(Q4__DOF_PER_NODE*Q4__DOF_PER_NODE,
        Q4__DOF_PER_NODE*Q4__DOF_PER_NODE);
    mat B(Q4__STRAIN_COMP,Q4__DOF_PER_NODE*Q4__NUM_NODES);
    mat C(Q4__STRAIN_COMP,Q4__STRAIN_COMP);
    
    /* TODO: separate this logic and material properties from element so that 
     *          the C matrix will not need calculated a billion times */
    // PROB AN EXPENSIVE OPERATION, DO SOMETHING ABOUT THIS
    switch (pState)
    {
        case PSTRESS:
            C << 1 << v << 0 << endr
              << v << 1 << 0 << endr
              << 0 << 0 << (1-v)/2 << endr;
            C *= E/(1-v*v);
            break;
        case PSTRAIN:
            C << 1-v << v << 0 << endr
              << v << 1-v << 0 << endr
              << 0 << 0 << (1-2*v)/2 << endr;
            C *= E/((1-v)*(1-2*v));
            break;
        default:
            // TODO: IDK, FIX THIS I GUESS
            cout << "Cannot understand plane state" << endl;
            assert(0);
    }
    
    /*
     * calculate element stiffness matrix using gauss quadrature
     */
    for (i = 0; i < Q4::numPointsX; i++) 
    {
        xi = Q4::gaussPoints[i];
        weightX = Q4::weights[i];
        /*
         * create "nStar" matrix
         */
        for (j = 1; j < Q4__DOF_PER_NODE*Q4__DOF_PER_NODE; j+=2)
        {
            nStar(j,j/2) = Q4__dN1_deta(xi);
            nStar(j,j/2+2) = Q4__dN2_deta(xi);
            nStar(j,j/2+4) = Q4__dN3_deta(xi);
            nStar(j,j/2+6) = Q4__dN4_deta(xi);
        }
        /*
         * look up 'eta'
         * look up weight in the y-direction
         * calculate jacobian
         * build JE matrix
         * finish building nStar
         */
        for (j = 0; j < Q4::numPointsY; j++) 
        {
            eta = Q4::gaussPoints[j];
            weightY = Q4::weights[j];
            jacobian = J(xi, eta);
            invJac = inv(jacobian);
            for (k = 0; k < Q4__DOF_PER_NODE; k++)
            {
                for (m = 0; m < Q4__DOF_PER_NODE; m++)
                {
                    temp = invJac(k,m);
                    JE(k,m) = temp;
                    JE(k+Q4__DOF_PER_NODE,m+Q4__DOF_PER_NODE) = temp;
                }
            }
            
            for (k = 0; k < Q4__DOF_PER_NODE*Q4__DOF_PER_NODE; k+=2)
            {
                nStar(k,k/2) = Q4__dN1_dxi(eta);
                nStar(k,k/2+2) = Q4__dN2_dxi(eta);
                nStar(k,k/2+4) = Q4__dN3_dxi(eta);
                nStar(k,k/2+6) = Q4__dN4_dxi(eta);
            }
            // this is derived in class notes (Aquino)
            B = e * JE * nStar;
            stiff += weightX * weightY * trans(B) * C * B * det(jacobian);
        }
    }
    stiff *= h;
}

//TODO implement Q4R element
