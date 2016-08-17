/*CHOMP class definition 
 *
 * Copyright (C) 2016 Rafael Valencia. All rights reserved.
 * License (3-Cluase BSD): https://github.com/rafaelvalencia
 * 
 * This code uses and is based on code from:
 *   Project: trychomp https://github.com/poftwaresatent/trychomp
 *   Copyright (C) 2014 Roland Philippsen. All rights reserved.
 *   License (3-Clause BSD) : https://github.com/poftwaresatent/trychomp
 * **
 * \file path_adaptor.hpp
 *
 * CHOMP for point vehicles (x,y) moving holonomously in the plane. It will
 * plan a trajectory (xi) connecting start point (qs) to end point (qe) while
 * avoiding obstacles (obs)
 * 
 * 
 * USAGE:
 * 			1. Declare a CHOMP object with your selected parameters (passed as arguments in contructor).
 * 			2. Initialize a path, starting pose and ending pose. E.g. use  'setPath', 'initStraightLinePath' or 'initStackedPath'.
 * 			3. Set obstacles with 'setObstacles' (add more with 'addObstacle').
 * 			4. Run CHOMP with 'generatePath' or run a single iteration with 'chompIteration' (you need an initialized path for this step).
 * 
 */
#ifndef CHOMP_H
#define CHOMP_H

#include <Eigen/Dense>
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::Isometry3d Transform;

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Base CHOMP class
class CHOMP
{
	protected:		
		// trajectory 
		Vector xi_;			// the trajectory (q_1, q_2, ...q_n)
		Vector qs_;			// the start config a.k.a. q_0
		Vector qe_;			// the end config a.k.a. q_(n+1)
		size_t nq_;			// number of robot poses q in xi
		size_t cdim_;		// dimension of config space
		size_t xidim_;		// dimension of trajectory, xidim = nq * cdim
		double dt_;	       	// time step
		double eta_; 		// >= 1, regularization factor for gradient descent
		double lambda_; 	// weight of smoothness objective (usually 10)
		
		//obstacles
		Matrix OBS_;		// Obstacles, each column contains (x,y,radius) of the obstacle
	
		// gradient descent 
		Matrix AA_;			// metric
		Vector bb_;			// acceleration bias for start and end config
		Matrix Ainv_;		// inverse of AA
		size_t  numIt_;		// Number of iterations
		double costGain_;   // Gain inside cost function (usually 10) 
		double res_;		// Residual
		size_t cter_;		// Count of optimization iterations
		bool PATH_INIT_; //Indicates whether a path has been initialized 
			
	public:
		CHOMP(double dt_input, double eta_input, double lambda_input, size_t nq_input, size_t cdim_input, size_t numIt_, double gain);
		void initCHOMP(void); 			//Initializes gradient descent vectors and matrices
		void generatePath(Vector  &xi);	//optimize path using with numIt_ iterations
		double chompIteration(Vector  &xi);	//performs a single gradient descent iteration
		double chompUpdate(Vector  &xi, double &U_cost, double &curv); //performs a single gradient descent iteration and gives more data
		double chompUpdateWithSearchRegion(Vector  &xi, double orientation, double &U_cost, double &curv);	//performs a single gradient descent iteration withing predefined search region	
		void setObstacles(Matrix obs);	//sets the obstacle matrix OBS_		
		void addObstacle(double px, double py, double radius);		//Appends a column in matrix OBS_ with the (x,y,radius) of the new obstacle
		void initStraightLinePath(Vector  &qs, Vector &qe);//initializes path xi as a straight line, with qs and qe, starting and ending point
		void initStackedPath(Vector  &qs, Vector &qe);//initializes path with qs and qe, starting and ending point, where xi_ is stacked to qs point
		void setPath(Vector  &qs, Vector &qe, Vector &xi); // sets an aribitrary value for xi_, qs_, qe_
		void getPath(Vector &xi); // gets xi_ 	
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Derived CHOMP SE2 class
/*
 * CHOMP class to work with planar differential-drive robot, i.e. it can translate and rotate 
 * and its workspace is still 2-dimensional, 
 * but Configuration Space is the special Euclidean group SE(2) = R2 × SO(2) 
 * (where SO(2) is the special orthogonal group of 2D rotations). 
 * Its configuration is represented using 3 parameters (x, y, orientation).
*/
class CHOMP_SE2: public CHOMP
{
   public:
   	CHOMP_SE2(double dt_input, double eta_input, double lambda_input, size_t nq_input, size_t numIt_, double gain);
	double chompIteration(Vector  &xi);	//performs a single gradient descent iteration considering nonholonomic constraints for a diff. drive robot.
	void boundTrajectoryAngles (void); // keep angles in trajectory inside pi to pi range
	void boundVectorAngles (Vector& VV); // keep angles in any trajectory vector inside pi to pi range
};

#endif
