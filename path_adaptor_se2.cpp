/* CHOMP class implementation
 *
 * Copyright (C) 2016 Rafael Valencia. All rights reserved.
 * License (3-Cluase BSD): https://github.com/rafaelvalencia
 * 
 * This code uses and is based on code from:
 *   Project: trychomp https://github.com/poftwaresatent/trychomp
 *   Copyright (C) 2014 Roland Philippsen. All rights reserved.
 *   License (3-Clause BSD) : https://github.com/poftwaresatent/trychomp
 * **
 * \file chomp.cpp
 *
 * CHOMP for point vehicles (x,y) moving holonomously in the plane. It will
 * plan a trajectory (xi) connecting start point (qs_) to end point (qe) while
 * avoiding obstacles (obs)
 */
#include "path_adaptor.hpp"
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdlib.h>
#include <sys/time.h>
#include <err.h>

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::Isometry3d Transform;
static size_t const obs_dim(3); 	// obstacle dimensions (x,y,Radius)

using namespace std;
#define PI (3.141592653589793)

////////////////////////////////////////////////////////////////////////////
// Auxiliar functions

void pi2pi(double &angle)
{
	double ang;
    /* Process angle*/     
    if( (angle < -2*PI) ||  (angle > 2* PI) ) 
		ang = fmod(angle, 2*PI); 
    else 
        ang = angle;
    if (ang > PI) ang = ang - 2*PI;
    if (ang < -PI) ang = ang + 2*PI;    
    angle = ang;    
}

////////////////////////////////////////////////////////////////////////////

CHOMP_SE2::CHOMP_SE2(double dt_input, double eta_input, double lambda_input, size_t nq_input, size_t numIt_input, double gain) : 
CHOMP(dt_input, eta_input, lambda_input, nq_input, 3, numIt_input, gain) 
{	
	cout << "-----------------------------------------------------:"<< endl;  		  	
	cout << "CHOMP SE2 started "<< endl;  
	cout << "-----------------------------------------------------:"<< endl;  	  	
}


void CHOMP_SE2::boundTrajectoryAngles (void)
{
	for (size_t iq (0); iq < nq_; ++iq) 
	{	  	   
		//Makes the angles to be between Pi to -Pi	   
		Vector q (xi_.block  (iq * cdim_, 0, cdim_, 1) ); 
		pi2pi ( q(2) );
		xi_.block  (iq * cdim_, 0, cdim_, 1) = q;
	}
}

void CHOMP_SE2::boundVectorAngles (Vector& VV)
{	 
    static size_t const VVsize = 	(VV.size()/cdim_) - 1;	
	for (size_t iq (0); iq < VVsize; ++iq) 
	{
		//Makes the angles to be between Pi to -Pi	   
		Vector q (VV.block  (iq * cdim_, 0, cdim_, 1) ); 
		pi2pi ( q(2));
		VV.block  (iq * cdim_, 0, cdim_, 1) = q;
	}
}

double CHOMP_SE2::chompIteration(Vector  &xi)
{  	
	
	// Before performing the iteration check if a path has been given	
	if (PATH_INIT_==false)
	{
		cout << "A path was not initialized. Leaving CHOMP iteration! " << endl;
		return NAN;	
	}
	
	
	
	//////////////////////////////////////////////////
	// beginning of "the" constrained CHOMP iteration
	 
	Vector nabla_smooth (AA_ * xi_ + bb_);   
	Vector const & xidd (nabla_smooth); // indeed, it is the same in this formulation...
  
	// Constrained CHOMP. 
	// Impose nonholonmic (NH) restrictions with the rolling constraint. 
	// Next we evaluate the constraint functional 
	// and its Jacobian b and C, respectively, 
	// as it appears in CHOMP's IJRR paper.
	//
	Matrix CC ( Matrix::Zero (xidim_, 1) );  
	double b = 0;  
	Vector c1 (Vector::Zero (3)); 
	Vector c2 (Vector::Zero (3));	 
	Vector cf1 (Vector::Zero (3)); 
	Vector cf2 (Vector::Zero (3));	   

	for (size_t iq (0); iq < (nq_-1); ++iq) 
	{
		   		    
		//Evaluate the constraint functional. It is defined by a sum of auxiliar functions
		//that depend on only two consecutive robot poses.
		Vector const q1 (xi_.block  (iq * cdim_, 0, cdim_, 1) ); 
		Vector const q2 (xi_.block  ((iq+1) * cdim_, 0, cdim_, 1) ); 

		double nhc =  ( q2(0)  -  q1(0) )* sin(q1(2)) - ( q2(1)  -  q1(1) )* cos(q1(2)); //nonholonomic constraint
		double fmc =   cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1)) + sqrt(pow(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1)),2)); //forward motion constraint
	  
		b +=  nhc * nhc + fmc; //square of the rolling constraint + forward motion constrain
	  
		//Computation of the Jacobian of the NH constraint. 
		//  Jacobian of the auxiliar functions  
		c1(0) = -2*sin(q1(2))*(cos(q1(2))*( q1(1)  -  q2(1) ) - sin(q1(2))*( q1(0)  -  q2(0) )); 
		c1(1) =  2*cos(q1(2))*(cos(q1(2))*( q1(1)  -  q2(1) ) - sin(q1(2))*( q1(0)  -  q2(0) )); 
		c1(2) =  -2*(cos(q1(2))*( q1(0)  -  q2(0) ) + sin(q1(2))*( q1(1)  -  q2(1) ))*(cos(q1(2))*( q1(1)  -  q2(1) ) - sin(q1(2))*( q1(0)  -  q2(0) ));  
		c2(0) =  2*sin(q1(2))*(cos(q1(2))*( q1(1)  -  q2(1) ) - sin(q1(2))*( q1(0)  -  q2(0) ));
		c2(1) = -2*cos(q1(2))*(cos(q1(2))*( q1(1)  -  q2(1) ) - sin(q1(2))*( q1(0)  -  q2(0) ));
		c2(2) =  0;
		//end of NHC Jacobian
	  
		//Computation of the Jacobian of the NH constraint. 
		//  Jacobian of the auxiliar functions 
		cf1(0)=     cos(q1(2)) + (cos(q1(2))*(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1))))/sqrt(pow(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1)),2)); 
		cf1(1)=     sin(q1(2)) + (sin(q1(2))*(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1))))/sqrt(pow(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1)),2)); 
		cf1(2)=     cos(q1(2))*(q1(1) - q2(1)) - sin(q1(2))*(q1(0) - q2(0)) + ((cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1)))*(cos(q1(2))*(q1(1) - q2(1)) - sin(q1(2))*(q1(0) - q2(0))))/sqrt(pow(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1)),2)); 
		cf2(0)=     -cos(q1(2)) - (cos(q1(2))*(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1))))/sqrt(pow(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1)),2)); 
		cf2(1) =    - sin(q1(2)) - (sin(q1(2))*(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1))))/sqrt(pow(cos(q1(2))*(q1(0) - q2(0)) + sin(q1(2))*(q1(1) - q2(1)),2)); 
		cf2(2)=     0;             
		//end of FMC Jacobian
	  
		if (iq == 0 ) 
		{
			  c2(0) +=  2*sin(qs_(2))*(cos(qs_(2))*( qs_(1)  -  q1(1) ) - sin(qs_(2))*( qs_(0)  -  q1(0) ));
			  c2(1) += -2*cos(qs_(2))*(cos(qs_(2))*( qs_(1)  -  q1(1) ) - sin(qs_(2))*( qs_(0)  -  q1(0) ));
			  c2(2) +=  0;
			  
			  cf2(0) +=     -cos(qs_(2)) - (cos(qs_(2))*(cos(qs_(2))*(qs_(0) - q1(0)) + sin(qs_(2))*(qs_(1) - q1(1))))/sqrt(pow(cos(qs_(2))*(qs_(0) - q1(0)) + sin(qs_(2))*(qs_(1) - q1(1)),2)); 
			  cf2(1) +=    - sin(qs_(2)) - (sin(qs_(2))*(cos(qs_(2))*(qs_(0) - q1(0)) + sin(qs_(2))*(qs_(1) - q1(1))))/sqrt(pow(cos(qs_(2))*(qs_(0) - q1(0)) + sin(qs_(2))*(qs_(1) - q1(1)),2)); 
			  cf2(2) +=     0;      
				   
			  double nhstart  =  ( q1(0)  -  qs_(0) )* sin(qs_(2)) - ( q1(1)  -  qs_(0) )* cos(qs_(2)); 
			  double fmstart =   cos(qs_(2))*(qs_(0) - q1(0)) + sin(qs_(2))*(qs_(1) - q1(1)) + sqrt(pow(cos(qs_(2))*(qs_(0) - q1(0)) + sin(qs_(2))*(qs_(1) - q1(1)),2)); //forward motion constraint  
			  b +=  nhstart * nhstart + fmstart; //square of the rolling constraint	          	          
		} 
	 
		//Update Jacobian with the contributions from the Jacobian of the auxiliar functions
		CC.block  (iq * cdim_, 0, cdim_, 1) = CC.block  (iq * cdim_, 0, cdim_, 1) + c1 + cf1;
		CC.block  ((iq+1) * cdim_, 0, cdim_, 1) = CC.block  ((iq+1) * cdim_, 0, cdim_, 1) + c2 + cf2;
	}
	 
	  
	CC /=  dt_ * dt_ * (nq_ + 1) ; 
	b /=  dt_ * dt_ * (nq_ + 1);
	
	Matrix CCtrans = CC.transpose(); 
    // end of the computation of C and b 
  
  
	Vector nabla_obs (Vector::Zero (xidim_));
  
    
	for (size_t iq (0); iq < nq_; ++iq) 
	{
		Vector const qq (xi_.block (iq * cdim_, 0, cdim_, 1));
		Vector qd;
		if (0 == iq) 
		{
			qd = 0.5 * (xi_.block ((iq+1) * cdim_, 0, cdim_, 1) - qs_);
		}
		else if (iq == nq_ - 1) 
		{
			qd = 0.5 * (qe_ - xi_.block ((iq-1) * cdim_, 0, cdim_, 1));
		}
		else 
		{
			qd = 0.5 * (xi_.block ((iq+1) * cdim_, 0, cdim_, 1) - xi_.block ((iq-1) * cdim_, 0, cdim_, 1));;
		}

		Vector const & xx (qq.block (0,0,2,1));
		Vector const & xd (qd.block (0,0,2,1));

		// In this case, C and W are NOT the same
		Matrix JJ(Matrix::Zero (2, 3));
		JJ(0,0)=1; 
		JJ(1,1)=1; 
	      	
		double const vel (xd.norm());
		if (vel < 1.0e-3) 
		{	
			// avoid div by zero further down
			continue;
		}
		Vector const xdn (xd / vel);

		Vector const xdd (JJ * xidd.block (iq * cdim_, 0, cdim_ , 1));

		Matrix const prj (Matrix::Identity (2, 2) - xdn * xdn.transpose()); // hardcoded planar case
		Vector const kappa (prj * xdd / pow (vel, 2.0));
		
		//Add obstacles			 	 
		for (int ii = 0; ii < OBS_.cols(); ii++) 
		{
			Vector delta(xx - OBS_.block(0, ii, 2, 1));
			double const dist(delta.norm());
			if ((dist >= OBS_(2, ii)) || (dist < 1e-9))
				continue;
			double const cost(costGain_ * OBS_(2, ii) * pow(1.0 - dist / OBS_(2, ii), 3.0) / 3.0);  
			delta *= - costGain_ *pow(1.0 - dist / OBS_(2, ii), 2.0) / dist;                        
			nabla_obs.block(iq * cdim_, 0, cdim_, 1) += JJ.transpose() * vel * (prj * delta - cost * kappa);
		} 
		 
	}
    double residual;
  
  
	Vector dxi (Ainv_ * (nabla_obs + lambda_ * nabla_smooth)); //unconstrained step

	if (b < 1.0e-10) 
	{ 
		//unconstrained update to initialize trajectory (it starts with b aprox to zero)
		//cout << " One unconstrained update to initialize trajectory  " << endl;
		xi_ -= dxi / eta_; 
		
		Vector dxi (Ainv_ * (nabla_obs + lambda_ * nabla_smooth));
		residual =  dxi.norm() / eta_;
	}
	else
	{
		// Constrained optimization update	
		Vector CAC (CCtrans * Ainv_ * CC);
	 
		Matrix CACinv (CAC.inverse());
	 
		Vector Proj (Ainv_ * CC * CACinv); //auxiliar matrix for the update equation
 
		//cout << "b =" << b << "\n";
		Vector cdxi ( - dxi/eta_  +  Proj*CCtrans*dxi / eta_ -  Proj * b  );

		xi_ += cdxi; //constrained update
		
		residual = cdxi.norm() / eta_;
	}
	// end of "the" constrainedCHOMP iteration
	//////////////////////////////////////////////////
	boundTrajectoryAngles();
	
	res_ = residual;
	xi = xi_; //updated path
	
	return res_;
 
}
