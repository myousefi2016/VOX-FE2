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

#include "voxfeRecomputeStrainFilter.h"
  
int voxfeRecomputeStrainFilter::
computeGradientMatrix( const double& voxel_size, BMatrix& m, int node ) {
  
  const double GP = 1.0/(2.0*sqrt(3.0)); // gauss points/2
  const double VR = 1.0/voxel_size;      // voxel size ratio
  
  MatrixXd integMtx(NODES_PER_ELEMENT,N_DOF);
  integMtx(0,0) = -GP; integMtx(0,1) = -GP; integMtx(0,2) = -GP;
  integMtx(1,0) = +GP; integMtx(1,1) = -GP; integMtx(1,2) = -GP;
  integMtx(2,0) = -GP; integMtx(2,1) = +GP; integMtx(2,2) = -GP;
  integMtx(3,0) = +GP; integMtx(3,1) = +GP; integMtx(3,2) = -GP;
  integMtx(4,0) = -GP; integMtx(4,1) = -GP; integMtx(4,2) = +GP;
  integMtx(5,0) = +GP; integMtx(5,1) = -GP; integMtx(5,2) = +GP;
  integMtx(6,0) = -GP; integMtx(6,1) = +GP; integMtx(6,2) = +GP;
  integMtx(7,0) = +GP; integMtx(7,1) = +GP; integMtx(7,2) = +GP;  
  
  
  double temp;
  for( int i=0; i<NUM_INTEG_POINTS; i++){ // 8 integration points
    
    temp = VR * VOXFESIGN(integMtx(i,0)) *                            //A
           ( 0.5 + integMtx(node,1) * VOXFESIGN(integMtx(i,1)) ) *
           ( 0.5 + integMtx(node,2) * VOXFESIGN(integMtx(i,2)) );
    m(0, 0+3*i) = temp;
    m(3, 1+3*i) = temp;
    m(5, 2+3*i) = temp;
    
    temp = VR * VOXFESIGN(integMtx(i,1)) *                            //B
           ( 0.5 + integMtx(node,0) * VOXFESIGN(integMtx(i,0)) ) *
           ( 0.5 + integMtx(node,2) * VOXFESIGN(integMtx(i,2)) );
    m(1, 1+3*i) = temp;
    m(3, 0+3*i) = temp;
    m(4, 2+3*i) = temp;
    
    
    temp = VR * VOXFESIGN(integMtx(i,2)) *                            //C
           ( 0.5 + integMtx(node,0) * VOXFESIGN(integMtx(i,0)) ) *
           ( 0.5 + integMtx(node,1) * VOXFESIGN(integMtx(i,1)) );
    m(2, 2+3*i) = temp;
    m(4, 1+3*i) = temp;
    m(5, 0+3*i) = temp;
    
  }
  
  //================ **NB** Factored down by 8 here ==============================
  m*= 0.125;
  
  return 0;
}




//================================================================================================

/*
 * Build Property Matrix
 *
 */
int voxfeRecomputeStrainFilter::
computePropertyMatrix(const double& ym, const double& pr, DMatrix& DM){

  DM.fill(0);

  DM(0, 0) = 1.0;
  DM(0, 1) = pr/(1.0-pr);
  DM(0, 2) = pr/(1.0-pr);
  DM(1, 0) = pr/(1.0-pr);
  DM(1, 1) = 1.0;
  DM(1, 2) = pr/(1.0-pr);
  DM(2, 0) = pr/(1.0-pr);
  DM(2, 1) = pr/(1.0-pr);
  DM(2, 2) = 1.0;
  DM(3, 3) = (1.0-2.0*pr)/(2.0*(1.0-pr));
  DM(4, 4) = (1.0-2.0*pr)/(2.0*(1.0-pr));
  DM(5, 5) = (1.0-2.0*pr)/(2.0*(1.0-pr));

  DM *= ym*(1.0-pr)/((1.0+pr)*(1.0-2.0*pr));

  return 0;
}


/*
 *      Solve:
 *              coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
 * 
 *   return false if no real roots found
 */
bool voxfeRecomputeStrainFilter::
FindCubicRoots(double coeff[4], double x[3])
{
  long double a1 = (long double)(coeff[2]) / (long double)(coeff[3]);
  long double a2 = (long double)(coeff[1]) / (long double)(coeff[3]);
  long double a3 = (long double)(coeff[0]) / (long double)(coeff[3]);

  long double Q = (a1 * a1 - 3 * a2) / 9;
  long double R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;

  long double Q3 = Q * Q * Q;
  long double R2 = R*R;
  const double pi = atan(1.)*4.; //3.141592653589;

  // Three real roots
  if (R2 <= Q3) {
    double temp = sqrt(Q3);
    double alpha = 0;
    if(temp>0)
      alpha = (R / sqrt(Q3));
    else alpha = 0;
    if(alpha < -1) alpha = -1;
    if(alpha > 1) alpha = 1;
    double theta = acos(alpha);
    double sqrtQ = sqrt(Q);
    x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
    x[1] = -2 * sqrtQ * cos((theta + 2 * pi) / 3) - a1 / 3;
    x[2] = -2 * sqrtQ * cos((theta - 2 * pi) / 3) - a1 / 3;
  }
  else {
    cerr << "\n\n** Roots of principal strains may be non-real... ***\n";
    cerr << "(Perhaps the strain tensor is very close to the solution?)\n\n"; 
    return false;
  }
  
  return true;
}


void  voxfeRecomputeStrainFilter::
computePrincipalStrain(double StrainTensor[6], double pstrain[3])
{
  // calculate the principal strain
  double I1 = StrainTensor[0] + StrainTensor[1] + StrainTensor[2];
  double factor = 4.0; // divide shear by 2
  double I2 = StrainTensor[0]*StrainTensor[1] +
    StrainTensor[1]*StrainTensor[2] +
    StrainTensor[2]*StrainTensor[0] -
    StrainTensor[3]*StrainTensor[3]/factor -
    StrainTensor[4]*StrainTensor[4]/factor -
    StrainTensor[5]*StrainTensor[5]/factor;
  double I3 = StrainTensor[0]*StrainTensor[1]*StrainTensor[2] +
    StrainTensor[3]*StrainTensor[4]*StrainTensor[5]/factor -
    StrainTensor[0]*(StrainTensor[3]*StrainTensor[3]/factor) -
    StrainTensor[1]*(StrainTensor[4]*StrainTensor[4]/factor) -
    StrainTensor[2]*(StrainTensor[5]*StrainTensor[5]/factor);
  double coeff[4] = {-I3,I2,-I1,1};
  double roots[3] = {0,0,0};
  
  if( ! FindCubicRoots(coeff,roots) ) {
    
    //fixme: if roots not found, set roots to normal strains
    // (ie presume only normal strain is present)
    for( int p=0; p<3; p++ ) roots[p] = StrainTensor[p];
  }
  
  // Ep1, find maximum value from cubic roots
  double Ep1 = roots[0]>roots[1]?roots[0]:roots[1];
  Ep1 = Ep1>roots[2]?Ep1:roots[2];

  // sort root values into Ep1 > Ep2 > Ep3
  if(Ep1 == roots[0])
    {
      if(roots[1]>roots[2]) {
        pstrain[0] = Ep1;
        pstrain[1] = roots[1];
        pstrain[2] = roots[2];
      }
      else {
        pstrain[0] = Ep1;
        pstrain[1] = roots[2];
        pstrain[2] = roots[1];
      }
    }
  else if (Ep1 == roots[1]) {
    if(roots[0] > roots[2]) {
      pstrain[0] = Ep1;
      pstrain[1] = roots[0];
      pstrain[2] = roots[2];
    }
    else {
      pstrain[0] = Ep1;
      pstrain[1] = roots[2];
      pstrain[2] = roots[0];
    }
  }
  else if(Ep1 == roots[2]) {
    if(roots[0]>roots[1]) {
      pstrain[0] = Ep1;
      pstrain[1] = roots[0];
      pstrain[2] = roots[1];
    }
    else {
      pstrain[0] = Ep1;
      pstrain[1] = roots[1];
      pstrain[2] = roots[0];
    }
  }
}

//================================================================================================

double voxfeRecomputeStrainFilter::
getModulusOfRigidity( const double& ym, const double& pr ) {
 
  return ym/( 2.*(1+pr));
}

double voxfeRecomputeStrainFilter::
getLambda( const double& ym, const double& pr ) {
  
 return ( ym * pr )/( (1.+pr) * (1.-(2.*pr)) ); 
}


//================================================================================================

/*
 * Compute the SED of an element using the normal strain components.
 */
void voxfeRecomputeStrainFilter::
computeStrainEnergyDensity(const double* strainTensor, const double& ym, const double& pr, double& sed) {
  
  double k = ym/(2.0*(1.0 + pr)*(1 - (2.0*pr)));
  
  double s0 = (strainTensor[0]*strainTensor[0] + strainTensor[1]*strainTensor[1] + strainTensor[2]*strainTensor[2]);
  double s1 = (strainTensor[0]*strainTensor[1] + strainTensor[1]*strainTensor[2] + strainTensor[2]*strainTensor[0]);
  
  sed = k * ( ((1.0 - pr) * s0) + (2.0 * pr * s1) );
}

/*
 * Compute SED as per Eq 7.8, Elasticity, Chou & Pagano
 * 
 */
void voxfeRecomputeStrainFilter::
computeStrainEnergyDensity2(const double* strainTensor, const double& rigidity_modulus, const double& lambda, double& sed) {
  
  double s0 = (strainTensor[0] + strainTensor[1] + strainTensor[2]);
  double s1 = (strainTensor[0]*strainTensor[0] + strainTensor[1]*strainTensor[1] + strainTensor[2]*strainTensor[2]);
  double s2 = (strainTensor[3]*strainTensor[3] + strainTensor[4]*strainTensor[4] + strainTensor[5]*strainTensor[5]);
  
  sed = 0.5 * ( (lambda*s0*s0) + (2.* rigidity_modulus * s1) + (rigidity_modulus * s2 ) );
}


//================================================================================================
//================================================================================================

//from Hartley Grandin jr, "Fundamentals of the Finite Element Method" 1991 p388
//using s,t,r

void voxfeRecomputeStrainFilter::
CreateIntegrationMatrix( IMatrix& Ig ) {
  
  const double GP = 1.0/sqrt(3.0); // gauss points
  
  Ig(0,0) = -GP; Ig(0,1) = -GP; Ig(0,2) = +GP;
  Ig(1,0) = +GP; Ig(1,1) = -GP; Ig(1,2) = +GP;
  Ig(2,0) = +GP; Ig(2,1) = +GP; Ig(2,2) = +GP;
  Ig(3,0) = -GP; Ig(3,1) = +GP; Ig(3,2) = +GP;
  Ig(4,0) = -GP; Ig(4,1) = -GP; Ig(4,2) = -GP;
  Ig(5,0) = +GP; Ig(5,1) = -GP; Ig(5,2) = -GP;
  Ig(6,0) = +GP; Ig(6,1) = +GP; Ig(6,2) = -GP;
  Ig(7,0) = -GP; Ig(7,1) = +GP; Ig(7,2) = -GP;  
}

void voxfeRecomputeStrainFilter::
CreateRefElement( IMatrix& RefElement ) {
  
  RefElement(0,0) = -1.; RefElement(0,1) = -1.; RefElement(0,2) = 1.;
  RefElement(1,0) =  1.; RefElement(1,1) = -1.; RefElement(1,2) = 1.;
  RefElement(2,0) =  1.; RefElement(2,1) =  1.; RefElement(2,2) = 1.;
  RefElement(3,0) = -1.; RefElement(3,1) =  1.; RefElement(3,2) = 1.;
  RefElement(4,0) = -1.; RefElement(4,1) = -1.; RefElement(4,2) = -1.;
  RefElement(5,0) =  1.; RefElement(5,1) = -1.; RefElement(5,2) = -1.;
  RefElement(6,0) =  1.; RefElement(6,1) =  1.; RefElement(6,2) = -1.;
  RefElement(7,0) = -1.; RefElement(7,1) =  1.; RefElement(7,2) = -1.;
  
}



double voxfeRecomputeStrainFilter::
sum_dNds(const IMatrix& Ig, const int& index) {
  double sum = 0; 
  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 
    const double & t = Ig(p,1);
    const double & r = Ig(p,2);
    switch(index) {
      case 1: sum += dN1ds(t,r); break;
      case 2: sum += dN2ds(t,r); break;
      case 3: sum += dN3ds(t,r); break;
      case 4: sum += dN4ds(t,r); break;
      case 5: sum += dN5ds(t,r); break;
      case 6: sum += dN6ds(t,r); break;
      case 7: sum += dN7ds(t,r); break;
      case 8: sum += dN8ds(t,r); break;
    }
  }
  return sum;
}

double voxfeRecomputeStrainFilter::
sum_dNdt(const IMatrix& Ig, const int& index) {
  double sum = 0; 
  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 
    const double & s = Ig(p,0);
    const double & r = Ig(p,2);
    switch(index) {
      case 1: sum += dN1dt(s,r); break;
      case 2: sum += dN2dt(s,r); break;
      case 3: sum += dN3dt(s,r); break;
      case 4: sum += dN4dt(s,r); break;
      case 5: sum += dN5dt(s,r); break;
      case 6: sum += dN6dt(s,r); break;
      case 7: sum += dN7dt(s,r); break;
      case 8: sum += dN8dt(s,r); break;
    }
  }
  return sum;
}

double voxfeRecomputeStrainFilter::
sum_dNdr(const IMatrix& Ig, const int& index) {
  double sum = 0; 
  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 
    const double & s = Ig(p,0);
    const double & t = Ig(p,1);
    switch(index) {
      case 1: sum += dN1dr(s,t); break;
      case 2: sum += dN2dr(s,t); break;
      case 3: sum += dN3dr(s,t); break;
      case 4: sum += dN4dr(s,t); break;
      case 5: sum += dN5dr(s,t); break;
      case 6: sum += dN6dr(s,t); break;
      case 7: sum += dN7dr(s,t); break;
      case 8: sum += dN8dr(s,t); break;
    }
  }
  return sum;
}


double voxfeRecomputeStrainFilter::
sum_dNds(const IMatrix& NiSums, double ord[8]) {
  
  double sum = 0; 
  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 
 
    sum += NiSums(p,0)*ord[p];
  }
  return sum;
}

double voxfeRecomputeStrainFilter::
sum_dNdt(const IMatrix& NiSums, double ord[8]) {
  double sum = 0;
  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 
    
    sum += NiSums(p,1)*ord[p];
  }
  return sum;
}

double voxfeRecomputeStrainFilter::
sum_dNdr(const IMatrix& NiSums, double ord[8]) {
  double sum = 0;
  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 
   
    sum += NiSums(p,2)*ord[p];
  }
  return sum;
}

void voxfeRecomputeStrainFilter::
getVoxelShapeFuncCoeffs( const double& voxel_size, Matrix3d& coeff, double& detJ, IMatrix& NiSums) {
  
  IMatrix Ig, refEl, unitEl;
  CreateIntegrationMatrix( Ig );
  CreateRefElement( refEl );

  //Here we use the voxel size to define a single element assumed to be used throughout.
  //If arbitrary brick elements were to be used 'unitEl' (and el_x etc) would have to be set
  //for each element and this function would need to be recomputed every time.
  unitEl = (voxel_size/2.) * refEl; // 2. is the ref element size, so we always need to scale
  //cout << "\nuser defined element: \n" << unitEl << "\n";
  
  double el_x[8], el_y[8], el_z[8];
  //=================================================
  //NiSums [8x3]... sum each individual dNi/d[s,t,r] term across all the 
  //integration points to arrive at dNi/ds, dNi/dt and dNi/dr for i=1,2,..8
  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 
    
    el_x[p] = unitEl(p,0);
    el_y[p] = unitEl(p,1);
    el_z[p] = unitEl(p,2);
    
    /*//debug
    cout << "\nx[" << p << "] = " << el_x[p] << "  " 
         << "  y[" << p << "] = " << el_y[p] << "  "
	 << "  z[" << p << "] = " << el_z[p] << "\n";*/
 
    NiSums(p,0) = sum_dNds(Ig, p+1);
    NiSums(p,1) = sum_dNdt(Ig, p+1);
    NiSums(p,2) = sum_dNdr(Ig, p+1);
  }
  
  //========================================================
  //Get the nine (a11 -- a33) coefficients to form the Jacobian
  coeff.fill(0);

  //a11
  double sum0 = 0., sum1 = 0., sum2 = 0.;
  sum0 = sum_dNdt( NiSums, el_y ) * sum_dNdr( NiSums, el_z );
  sum1 = sum_dNdt( NiSums, el_z ) * sum_dNdr( NiSums, el_y );
  coeff(0,0) = sum0 - sum1;
  
  //a12
  sum0 = sum_dNdt( NiSums, el_z ) * sum_dNdr( NiSums, el_x );
  sum1 = sum_dNdt( NiSums, el_x ) * sum_dNdr( NiSums, el_z );
  coeff(0,1) = sum0 - sum1;
  
  //a13
  sum0 = sum_dNdt( NiSums, el_x ) * sum_dNdr( NiSums, el_y );
  sum1 = sum_dNdt( NiSums, el_y ) * sum_dNdr( NiSums, el_x );
  coeff(0,2) = sum0 - sum1;
  
  //a21
  sum0 = sum_dNds( NiSums, el_z ) * sum_dNdr( NiSums, el_y );
  sum1 = sum_dNds( NiSums, el_y ) * sum_dNdr( NiSums, el_z );
  coeff(1,0) = sum0 - sum1;  

  //a22
  sum0 = sum_dNds( NiSums, el_x ) * sum_dNdr( NiSums, el_z );
  sum1 = sum_dNds( NiSums, el_z ) * sum_dNdr( NiSums, el_x );
  coeff(1,1) = sum0 - sum1;  

  //a23
  sum0 = sum_dNds( NiSums, el_y ) * sum_dNdr( NiSums, el_x );
  sum1 = sum_dNds( NiSums, el_x ) * sum_dNdr( NiSums, el_y );
  coeff(1,2) = sum0 - sum1;    

  //a31
  sum0 = sum_dNds( NiSums, el_y ) * sum_dNdt( NiSums, el_z );
  sum1 = sum_dNds( NiSums, el_z ) * sum_dNdt( NiSums, el_y );
  coeff(2,0) = sum0 - sum1;  

  //a32
  sum0 = sum_dNds( NiSums, el_z ) * sum_dNdt( NiSums, el_x );
  sum1 = sum_dNds( NiSums, el_x ) * sum_dNdt( NiSums, el_z );
  coeff(2,1) = sum0 - sum1;  

  //a33
  sum0 = sum_dNds( NiSums, el_x ) * sum_dNdt( NiSums, el_y );
  sum1 = sum_dNds( NiSums, el_y ) * sum_dNdt( NiSums, el_x );
  coeff(2,2) = sum0 - sum1;
  
  //cout << "\n\nCoeffs matrix: " << coeff << "\n\n";
  
  //=================================================
  //compute detJ
  sum0 = sum_dNds( NiSums, el_x );
  sum1 = sum_dNds( NiSums, el_y );
  sum2 = sum_dNds( NiSums, el_z );
  detJ = coeff(0,0)*sum0 + coeff(0,1)*sum1 + coeff(0,2)*sum2;

}

void voxfeRecomputeStrainFilter::
computeVoxelGradientMatrix( const double& voxel_size, BMatrix& B ) {
  
  double detJ;
  Matrix3d coeffs;
  IMatrix NiSums;
  getVoxelShapeFuncCoeffs( voxel_size, coeffs, detJ, NiSums );
  
  //debug
  //cout << "\n coeffs = \n" << coeffs << "\n";
  //cout << "\n detJ = "     << detJ << "\n";
  //cout << "\n NiSums = "   << NiSums << "\n";
  
  B.fill(0);    //construct B
#ifdef USE_HARTLEY_GRANDIN_NODE_ORDER

  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 
    //                                           's'                       't'                       'r'
    B(0,3*p  ) = (1./detJ)*( coeffs(0,0)*NiSums(p,0) + coeffs(0,1)*NiSums(p,1) + coeffs(0,2)*NiSums(p,2) );
    B(1,3*p+1) = (1./detJ)*( coeffs(1,0)*NiSums(p,0) + coeffs(1,1)*NiSums(p,1) + coeffs(1,2)*NiSums(p,2) );
    B(2,3*p+2) = (1./detJ)*( coeffs(2,0)*NiSums(p,0) + coeffs(2,1)*NiSums(p,1) + coeffs(2,2)*NiSums(p,2) );
    
    B(3,3*p  ) = B(1,3*p+1);
    B(3,3*p+1) = B(0,3*p  );
    
    B(4,3*p+1) = B(2,3*p+2);
    B(4,3*p+2) = B(1,3*p+1);
    
    B(5,3*p  ) = B(2,3*p+2);
    B(5,3*p+2) = B(0,3*p  );
  } 
#else
  // use VoxFE node order  --  NB to change node order we have to change 
  // where things get stored( eg. B(x,y) ), not what gets computed   
  int q[] = { 4,5,7,6,0,1,3,2 };
  double b11,b22,b33; //temp vars
  for( int p=0; p<NODES_PER_ELEMENT; p++ ) { 

    //                                    's'                       't'                       'r'
    b11 = (1./detJ)*( coeffs(0,0)*NiSums(p,0) + coeffs(0,1)*NiSums(p,1) + coeffs(0,2)*NiSums(p,2) );
    b22 = (1./detJ)*( coeffs(1,0)*NiSums(p,0) + coeffs(1,1)*NiSums(p,1) + coeffs(1,2)*NiSums(p,2) );
    b33 = (1./detJ)*( coeffs(2,0)*NiSums(p,0) + coeffs(2,1)*NiSums(p,1) + coeffs(2,2)*NiSums(p,2) );
    
    B(0,3*q[p]  ) = b11;
    B(1,3*q[p]+1) = b22;
    B(2,3*q[p]+2) = b33;
    
    B(3,3*q[p]  ) = b22;
    B(3,3*q[p]+1) = b11;
    
    B(4,3*q[p]+1) = b33;
    B(4,3*q[p]+2) = b22;
    
    B(5,3*q[p]  ) = b33;
    B(5,3*q[p]+2) = b11;
  } 
#endif

}



//================================================================================================
#ifdef VOXFE_TEST_MAIN

int main()
{
  VectorXd s0(6); s0.fill(0);
  VectorXd s1(6); s1.fill(0);
  VectorXd disp0(24), disp1(24);
  
  double ym = 1e6;
  double pr = 0.3;

//#define HALF_CUBE  1
#ifdef HALF_CUBE
  
  const double A_SIZE = 0.5;
  const double B_SIZE = 0.5;
  const double C_SIZE = 0.5;
 
  disp0 << 0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,
          -1.49999753853241166556e-07,0.00000000000000000000e+00,0.00000000000000000000e+00,
           0.00000000000000000000e+00,-1.49999753853040975613e-07,0.00000000000000000000e+00,
          -1.50000282096109980103e-07,-1.50000282095894198459e-07,0.00000000000000000000e+00,
          -1.65930904772066341399e-13,-1.65931025539529025862e-13,4.99999773391260189412e-07,
          -1.50000115041848496650e-07,9.51130335688584924537e-14,5.00000059163966182348e-07,
           9.51134165929765961194e-14,-1.50000115041425747800e-07,5.00000059164154011904e-07,
          -1.49999943515780303234e-07,-1.49999943515310438176e-07,4.99999951588959489666e-07;

  disp1 << -1.50000282096109980103e-07,-1.50000282095894198459e-07,0.00000000000000000000e+00,
           -3.00000343653272903065e-07,-1.49999780086494282210e-07,0.00000000000000000000e+00,
           -1.49999780086673350269e-07,-3.00000343653513989818e-07,0.00000000000000000000e+00,
           -2.99999776412594851704e-07,-2.99999776412727200602e-07,0.00000000000000000000e+00,
           -1.49999943515780303234e-07,-1.49999943515310438176e-07,4.99999951588959489666e-07,
           -2.99999882482599385875e-07,-1.49999914491547855965e-07,5.00000031273934780133e-07,
           -1.49999914491961102164e-07,-2.99999882482982297707e-07,5.00000031274148338315e-07,
           -2.99999918433791148703e-07,-2.99999918433907351036e-07,4.99999835114742273521e-07;
	  
  
#else  
	  
  const double A_SIZE = 1.;
  const double B_SIZE = 1.;
  const double C_SIZE = 1.;	  

#ifdef USE_HARTLEY_GRANDIN_NODE_ORDER  
  //As above but node order changed to that used by Hartley Grandin
  disp0 <<  0,0,1e-06,
	    -3e-07,0,1e-06, 
	    -3e-07,-3e-07,1e-06,
	    0,-3e-07,1e-06,   
	    0, 0, 0, 
	    -3e-07, 0, 0,  
	    -3e-07,-3e-07, 0,
	    0,-3e-07, 0 ;
	    
  disp1 << 0,0,2e-06,
	    -3e-07,0,2e-06,
	    -3e-07,-3e-07,2e-06,
	    0,-3e-07,2e-06,
            0,0,1e-06,
	    -3e-07,0,1e-06,
	    -3e-07,-3e-07,1e-06,
	    0,-3e-07,1e-06;
#else
  disp0 <<  0, 0, 0, 
	    -3e-07, 0, 0,  
	    0,-3e-07, 0, 
	    -3e-07,-3e-07, 0, 
	    0,0,1e-06,
	    -3e-07,0,1e-06, 
	    0,-3e-07,1e-06,   
	    -3e-07,-3e-07,1e-06 ;
	    
  disp1 << 0,0,1e-06,
	    -3e-07,0,1e-06, 
	    0,-3e-07,1e-06,   
	    -3e-07,-3e-07,1e-06,
	    0,0,2e-06,
	    -3e-07,0,2e-06,
	    0,-3e-07,2e-06,
	    -3e-07,-3e-07,2e-06;
#endif	    

  

#endif	    
	    
  std::cout << "\n==============================================\n";
  voxfeRecomputeStrainFilter recomp;
  

  voxfeRecomputeStrainFilter::BMatrix m0, B0;
  m0.setZero();
  for(int n=0; n<NODES_PER_ELEMENT; ++n) {
  
    recomp.computeGradientMatrix(A_SIZE, m0, n);
    //std::cout << n << " m0 :\n" << m0.transpose() << "\n\n";
    
    s0 += (m0 * disp0);
  }
  std::cout << " strain tensor (el 0) :\n" << s0.transpose() << "\n\n";
  
  
  std::cout << "\n==============================================\n";  

  //second element   
  m0.setZero();
  s0.setZero();
  B0.setZero();
  for(int n=0; n<NODES_PER_ELEMENT; ++n) {

    recomp.computeGradientMatrix(A_SIZE, m0, n);
    //std::cout << n << " m0 :\n" << m0.transpose() << "\n\n";
    
    s0 += (m0 * disp1);
    B0 += m0;
  }
  std::cout << "B1 matrix: \n" << B0.transpose() << "\n\n";
  std::cout << " strain tensor (el 1) :\n" << s0.transpose() << "\n\n";
  
  std::cout << "\n==============================================\n";

  voxfeRecomputeStrainFilter::DMatrix dm;
  recomp.computePropertyMatrix(ym, pr, dm);
  
  //std::cout << " d matrix :\n" << dm << "\n\n";
  
  std::cout << "\n==============================================\n";
  
  double pstrain[3];
  recomp.computePrincipalStrain(s0.data(), pstrain);
  
  std::cout << " p strains :\n" << pstrain[0] << "  " << pstrain[1] << "  " << pstrain[2] << "\n\n";
  
  std::cout << "\n==============================================\n";
  
  double sed;
  recomp.computeStrainEnergyDensity(s0.data(), ym, pr, sed);
  
  std::cout << " SED :\n" << sed << "\n\n";
  
  
  std::cout << "\n=================================================================\n";
  std::cout << "\n=================================================================\n";
  

  voxfeRecomputeStrainFilter::BMatrix B5;
  recomp.computeVoxelGradientMatrix( A_SIZE, B5 );

  cout << "B5 matrix: \n" << B5.transpose() << "\n\n";
  
  cout << "\n strain tensor (B5 * disp0):\n" << (B5 * disp0).transpose() << "\n\n";
  cout << "\n strain tensor (B5 * disp1):\n" << (B5 * disp1).transpose() << "\n\n";  
  
  s1 = B5 * disp0;
  cout << " HG strain tensor:\n" << s1.transpose() << "\n\n";
  
  
  recomp.computePrincipalStrain(s1.data(), pstrain);
  cout << " HG p-strains :\n" << pstrain[0] << "  " << pstrain[1] << "  " << pstrain[2] << "\n\n";  
  
  recomp.computeStrainEnergyDensity(s1.data(), ym, pr, sed);
  cout << " HG SED :\n" << sed << "\n\n";  
  
  double rigidity_modulus = recomp.getModulusOfRigidity( ym, pr );
  double lambda           = recomp.getLambda( ym, pr );
  recomp.computeStrainEnergyDensity2( s1.data(), rigidity_modulus, lambda, sed);
  cout << " HG SED (2):\n" << sed << "\n\n";
  
  return 0;
}

#endif

