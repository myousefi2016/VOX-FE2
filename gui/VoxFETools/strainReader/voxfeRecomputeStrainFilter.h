/*
 * Vox-FE (Version 2): Voxel-based Finite Element Modelling
 *
 * Interface developed as a ParaView plugin by
 * Richard Holbrey and Michael J. Fagan
 *
 * Dept of Mechanical Engineering
 * University of Hull
 * Cottingham Road
 * HU6 7RX
 *
 * Vox-FE (Version 1) was created by Andreas Bitternas, George Sisias, Jia Liu and others.
 */

#ifndef VOXFE_RECOMPUTE_STRAIN_FILTER_H
#define VOXFE_RECOMPUTE_STRAIN_FILTER_H

#include <iostream>
#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::RowMajor;
using std::cout;
using std::cerr;

#define BRICK_EL_SHAPE_FUNC_SCALE 0.125
#define NODES_PER_ELEMENT 8
#define NUM_INTEG_POINTS  8
#define N_DOF 3
#define VOXFESIGN(x) ((x > 0) - (x < 0))

/** \class voxfeRecomputeStrainFilter
    \brief Compute strains from displacement data
*/
class voxfeRecomputeStrainFilter {
public:
	typedef Matrix <double, 6,24,RowMajor> BMatrix;
	typedef Matrix <double, 6,6, RowMajor> DMatrix;
	typedef Matrix <double, NODES_PER_ELEMENT,N_DOF, RowMajor> IMatrix;
	typedef Matrix <double, 24,1>          uVector;
	typedef Matrix <double, 6,1>           eVector;

	/**  Build Gradient Matrix (Neelofer Banglawala's code)
	 *
	 * Only need voxel size A, B and C for this (assumed same here)
	 * NB this produces matrix which is a factor 8 x method of Grandin, so need to do 0.125*B
	 */
	int computeGradientMatrix( const double& voxel_size, BMatrix& m, int node );

	/** Build Property Matrix (aka 'E') */
	int computePropertyMatrix(const double& ym, const double& pr, DMatrix& DM);

	//================================================================================================

	/** Compute the gradient matrix (after Hartley Grandin) */
	void computeVoxelGradientMatrix( const double& voxel_size, BMatrix& B );

        /** Compute principal strains from the strain tensor */
	void  computePrincipalStrain(double StrainTensor[6], double pstrain[3]);

	/** strain (tensor) = B*u */
	inline void getStrainTensor( const BMatrix& B, const uVector& u, eVector& strain ) { strain = B * u; }

	//================================================================================================
	/** Compute G lame coeff */
	double getModulusOfRigidity( const double& ym, const double& pr );

	/** Compute lambda lame coeff */
	double getLambda( const double& ym, const double& pr );

	//================================================================================================

	/** Compute the SED of an element using the normal strain components.
	 */
	void computeStrainEnergyDensity(const double* strainTensor, const double& ym, const double& pr, double& sed);

	/** Compute SED as per Eq 7.8, Elasticity, Chou & Pagano  */
	void computeStrainEnergyDensity2(const double* strainTensor, const double& rigidity_modulus,
		                               const double& lambda, double& sed);


protected:

	/** Find cubic roots
	 *  \brief Solve: coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
	 *  \return false if no real roots found
	 */
	bool FindCubicRoots(double coeff[4], double x[3]);



	//================================================================================================

	/** Compute Voxel element Integration matrix
	 * from Hartley Grandin jr, "Fundamentals of the Finite Element Method" 1991 p388
	 * using s,t,r
	 */
	void CreateIntegrationMatrix( IMatrix& Ig );

	/** Create a matrix of reference brick element coords */
	void CreateRefElement( IMatrix& RefElement );

	inline double dN1ds(const double& t, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(t-1.)*(1.+r); }   ///< Shape func (after Grandin)
	inline double dN1dt(const double& s, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(s-1.)*(1.+r); }   ///< Shape func (after Grandin)
	inline double dN1dr(const double& s, const double& t) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.-s)*(1.-t); }   ///< Shape func (after Grandin)

	inline double dN2ds(const double& t, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.-t)*(1.+r); }   ///< Shape func (after Grandin)
	inline double dN2dt(const double& s, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(-s-1.)*(1.+r); }   ///< Shape func (after Grandin)
	inline double dN2dr(const double& s, const double& t) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.+s)*(1.-t); }   ///< Shape func (after Grandin)

	inline double dN3ds(const double& t, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.+t)*(1.+r); }   ///< Shape func (after Grandin)
	inline double dN3dt(const double& s, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(s+1.)*(1.+r); }   ///< Shape func (after Grandin)
	inline double dN3dr(const double& s, const double& t) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.+s)*(1.+t); }   ///< Shape func (after Grandin)

	inline double dN4ds(const double& t, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(-t-1.)*(1.+r); }   ///< Shape func (after Grandin)
	inline double dN4dt(const double& s, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.-s)*(1.+r); }   ///< Shape func (after Grandin)
	inline double dN4dr(const double& s, const double& t) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.-s)*(1.+t); }   ///< Shape func (after Grandin)

	inline double dN5ds(const double& t, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(t-1.)*(1.-r); }   ///< Shape func (after Grandin)
	inline double dN5dt(const double& s, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(s-1.)*(1.-r); }   ///< Shape func (after Grandin)
	inline double dN5dr(const double& s, const double& t) { return BRICK_EL_SHAPE_FUNC_SCALE*(s-1.)*(1.-t); }   ///< Shape func (after Grandin)

	inline double dN6ds(const double& t, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.-t)*(1.-r); }   ///< Shape func (after Grandin)
	inline double dN6dt(const double& s, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(s+1.)*(r-1.); }   ///< Shape func (after Grandin)
	inline double dN6dr(const double& s, const double& t) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.+s)*(t-1.); }   ///< Shape func (after Grandin)

	inline double dN7ds(const double& t, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(t+1.)*(1.-r); }   ///< Shape func (after Grandin)
	inline double dN7dt(const double& s, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(s+1.)*(1.-r); }   ///< Shape func (after Grandin)
	inline double dN7dr(const double& s, const double& t) { return BRICK_EL_SHAPE_FUNC_SCALE*(-1.-s)*(1.+t); }   ///< Shape func (after Grandin)

	inline double dN8ds(const double& t, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(t+1.)*(r-1.); }   ///< Shape func (after Grandin)
	inline double dN8dt(const double& s, const double& r) { return BRICK_EL_SHAPE_FUNC_SCALE*(1.-s)*(1.-r); }   ///< Shape func (after Grandin)
	inline double dN8dr(const double& s, const double& t) { return BRICK_EL_SHAPE_FUNC_SCALE*(s-1.)*(1.+t); }   ///< Shape func (after Grandin)


	/** Sum across the integration points to get each dNx/ds term for index x */
	double sum_dNds(const IMatrix& Ig, const int& index);

	/** Sum across the integration points to get each dNx/dt term for index x */
	double sum_dNdt(const IMatrix& Ig, const int& index);

	/** Sum across the integration points to get each dNx/dr term for index x */
	double sum_dNdr(const IMatrix& Ig, const int& index);

	/** Sum the product of dN[]/ds and the relevant x,y,z ordinates */
	double sum_dNds(const IMatrix& NiSums, double ord[8]);

	/** Sum the product of dN[]/dt and the relevant x,y,z ordinates */
	double sum_dNdt(const IMatrix& NiSums, double ord[8]);

	/** Sum the product of dN[]/dr and the relevant x,y,z ordinates */
	double sum_dNdr(const IMatrix& NiSums, double ord[8]);

	/** Compute the 'a' coeffs of Hartley Grandin Jacobian matrix */
	void getVoxelShapeFuncCoeffs( const double& voxel_size, Matrix3d& coeff, double& detJ, IMatrix& NiSums);


};


#endif
