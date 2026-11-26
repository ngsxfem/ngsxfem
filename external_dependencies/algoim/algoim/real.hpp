#ifndef ALGOIM_REAL_HPP
#define ALGOIM_REAL_HPP

// Header file for defining algoim::real; by default, algoim::real is a typedef for double.
// High precision arithmetic can be used with Algoim by redefining algoim::real as another
// type (e.g., the quadruple-double type qd_real as provided by the QD library). Additional
// modifications may be necessary depending on the end application; e.g., the quadrature
// methods would need high precision versions of eigenvalue and SVD solvers; for further
// details and suggestions, see Algoim's GitHub page.

namespace algoim
{
    // typedef for algoim::real
    using real = double;
} // namespace algoim

#endif
