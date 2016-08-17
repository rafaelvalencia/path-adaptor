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

/**
   \file chomp-class.cpp
   
   Interactive demo with CHOMP class for point vehicles moving
   holonomously in the plane.  There is a fixed start and goal
   configuration, and you can drag a circular obstacle around to see
   how the CHOMP algorithm reacts to that.  
*/

#include "gfx.hpp"
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdlib.h>
#include <sys/time.h>
#include <err.h>
#include "path_adaptor.hpp"

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::Isometry3d Transform;

static size_t const obs_dim(3); 	// obstacle dimensions (x,y,Radius)

using namespace std;


//////////////////////////////////////////////////
// trajectory etc
Vector trajectory;				// the trajectory (q_1, q_2, ...q_n)
Vector trajectory2;				// the trajectory (q_1, q_2, ...q_n)
Matrix  obs(obs_dim,2);   		// Matrix containninging all obstacles, each column is (x,y,R) of the obstacle
Vector qs;						// the start config a.k.a. q_0
Vector qe;						// the end config a.k.a. q_(n+1)
static size_t const cdim (2);	// dimension of config space
static size_t const nq (20);	// number of q stacked into xi
static size_t const xidim (nq * cdim); // dimension of trajectory, xidim = nq * cdim

double cost1, cost2, curv1, curv2;
//////////////////////////////////////////////////
// CHOMP
CHOMP chomp(1.0, 100.0, 1.0, nq, cdim, 1000, 10.0);

/* CHOMP parameters
dt     = 1.0 	(time step)
eta    = 100.0 	(regularization factor for gradient descent)
lambda = 1.0	(weight of smoothness objective)
nq     = 20		(number of poses q in xi)
cdim   = 2		(dimension of config space)
numIt  = 1000	(number of iterations)
gain   = 10.0	(Gain inside cost function)
*/
//////////////////////////////////////////////////
// Auxiliar CHOMP
CHOMP aux_chomp(1.0, 100.0, 1.0, nq, cdim, 1000, 10.0);

/* CHOMP parameters
dt     = 1.0 	(time step)
eta    = 100.0 	(regularization factor for gradient descent)
lambda = 1.0	(weight of smoothness objective)
nq     = 20		(number of poses q in xi)
cdim   = 2		(dimension of config space)
numIt  = 1000	(number of iterations)
gain   = 10.0	(Gain inside cost function)
*/


//////////////////////////////////////////////////
// gui stuff

enum { PAUSE, STEP, RUN } state;

struct handle_s {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  handle_s (double px, double py, double radius, double red, double green, double blue, double alpha)
    : point_(2),
      radius_(radius),
      red_(red),
      green_(green),
      blue_(blue),
      alpha_(alpha)
  {
    point_ << px, py;
  }
  
  Vector point_;
  double radius_, red_, green_, blue_, alpha_;
};

// Set two obstacles
static handle_s rep1 (3.0, 0.0,   2.0, 0.0, 0.0, 1.0, 0.2);
static handle_s rep2 (0.0, 3.0,   2.0, 0.0, 0.5, 1.0, 0.2);

static handle_s * handle[] = { &rep1, &rep2, 0 };
static handle_s * grabbed (0);
static Vector grab_offset (3);



//////////////////////////////////////////////////
// robot (one per waypoint)

class Robot
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  Robot ()
    : position_ (Vector::Zero(2))
  {
  }
  
  
  void update (Vector const & position)
  {
    if (position.size() != 2) {
      errx (EXIT_FAILURE, "Robot::update(): position has %zu DOF (but needs 2)",
	    (size_t) position.size());
    }
    position_ = position;
  }
  
  
  void draw () const
  {
    // translucent disk for base
    gfx::set_pen (1.0, 0.7, 0.7, 0.7, 0.5);
    gfx::fill_arc (position_[0], position_[1], radius_, 0.0, 2.0 * M_PI);
    
    // thick circle outline for base
    gfx::set_pen (3.0, 0.2, 0.2, 0.2, 1.0);
    gfx::draw_arc (position_[0], position_[1], radius_, 0.0, 2.0 * M_PI);
  }
  
  static double const radius_;
  
  Vector position_;
};

double const Robot::radius_ (0.5);

Robot rstart;
Robot rend;
vector <Robot> robots;


//////////////////////////////////////////////////

static void update_robots ()
{	
  rstart.update (qs);
  rend.update (qe);
  if (nq != robots.size()) 
  {
    robots.resize (nq);
  }  
  for (size_t ii (0); ii < nq; ++ii) 
  {
    robots[ii].update (trajectory.block (ii * cdim, 0, cdim, 1));
  }
}


static void cb_step ()
{
  state = STEP;
}


static void cb_run ()
{
  if (RUN == state) {
    state = PAUSE;
  }
  else {
    state = RUN;
  }
}


static void cb_jumble ()
{
  for (size_t ii (0); ii < xidim; ++ii) 
  {
    trajectory[ii] = double (rand()) / (0.1 * numeric_limits<int>::max()) - 5.0;
    trajectory2[ii] = double (rand()) / (0.1 * numeric_limits<int>::max()) - 5.0;
  }
  //update new initial path
  chomp.setPath(qs, qe, trajectory);
  aux_chomp.setPath(qs, qe, trajectory2);
  update_robots();
}


static void cb_idle ()
{
  if (PAUSE == state) 
  {
    return;
  }
  if (STEP == state) 
  {
    state = PAUSE;
  }
  
  //static size_t stepcounter (0);
  //cout << "step " << stepcounter++ << "\n";
   
  /////////////////////////////////////////////////
  // Get updated obstacles
    //cout << "Getting updated obstacles ..." << endl;
    int ii = 0;
	for (handle_s ** hh (handle); *hh != 0; ++hh) 
	{
		if (ii < obs.cols())
		{
			Vector position ( (*hh)->point_ );		
			obs(0,ii) = position(0);
			obs(1,ii) = position(1);
			obs(2,ii) = (*hh)->radius_;
			++ii;
		}	 
	}
  //cout << "OBSTACLES : \n" << obs << endl;
  
  //////////////////////////////////////////////////
  // One CHOMP iteration  
  chomp.setObstacles(obs); //update obstacles  
  aux_chomp.setObstacles(obs); //update obstacles  
    
  
  //chomp.chompIteration(trajectory); //update path
  
  //update path considering a selected search region. Assuming robot's orientation = 45deg
  
  chomp.chompUpdateWithSearchRegion(trajectory, 0.79, cost1, curv1); //update path 
  aux_chomp.chompUpdateWithSearchRegion(trajectory2, 0.79, cost2, curv2); //update path 
  
  // if axiliar chomp gives lower cost use its result
  if(cost2 < cost1)
  {
	  cout << "  Using second CHOMP for this iteration!!! " << endl;
	  trajectory = trajectory2;
  }
  
  update_robots ();
}


static void cb_draw ()
{
  //////////////////////////////////////////////////
  // set bounds
  
  Vector bmin (qs);
  Vector bmax (qs);
  for (size_t ii (0); ii < 2; ++ii) {
    if (qe[ii] < bmin[ii]) {
      bmin[ii] = qe[ii];
    }
    if (qe[ii] > bmax[ii]) {
      bmax[ii] = qe[ii];
    }
    for (size_t jj (0); jj < nq; ++jj) {
      if (trajectory[ii + cdim * jj] < bmin[ii]) {
	bmin[ii] = trajectory[ii + cdim * jj];
      }
      if (trajectory[ii + cdim * jj] > bmax[ii]) {
	bmax[ii] = trajectory[ii + cdim * jj];
      }
    }
  }
  
  gfx::set_view (bmin[0] - 2.0, bmin[1] - 2.0, bmax[0] + 2.0, bmax[1] + 2.0);
  
  //////////////////////////////////////////////////
  // robots
  
  rstart.draw();
  for (size_t ii (0); ii < robots.size(); ++ii) {
    robots[ii].draw();
  }
  rend.draw();
  
  //////////////////////////////////////////////////
  // trj
  
  gfx::set_pen (1.0, 0.2, 0.2, 0.2, 1.0);
  gfx::draw_line (qs[0], qs[1], trajectory[0], trajectory[1]);
  for (size_t ii (1); ii < nq; ++ii) {
    gfx::draw_line (trajectory[(ii-1) * cdim], trajectory[(ii-1) * cdim + 1], trajectory[ii * cdim], trajectory[ii * cdim + 1]);
  }
  gfx::draw_line (trajectory[(nq-1) * cdim], trajectory[(nq-1) * cdim + 1], qe[0], qe[1]);
  
  gfx::set_pen (5.0, 0.8, 0.2, 0.2, 1.0);
  gfx::draw_point (qs[0], qs[1]);
  gfx::set_pen (5.0, 0.5, 0.5, 0.5, 1.0);
  for (size_t ii (0); ii < nq; ++ii) {
    gfx::draw_point (trajectory[ii * cdim], trajectory[ii * cdim + 1]);
  }
  gfx::set_pen (5.0, 0.2, 0.8, 0.2, 1.0);
  gfx::draw_point (qe[0], qe[1]);
  
  //////////////////////////////////////////////////
  // handles
  
  for (handle_s ** hh (handle); *hh != 0; ++hh) {
    gfx::set_pen (1.0, (*hh)->red_, (*hh)->green_, (*hh)->blue_, (*hh)->alpha_);
    gfx::fill_arc ((*hh)->point_[0], (*hh)->point_[1], (*hh)->radius_, 0.0, 2.0 * M_PI);
  }
}


static void cb_mouse (double px, double py, int flags)
{
  if (flags & gfx::MOUSE_PRESS) {
    for (handle_s ** hh (handle); *hh != 0; ++hh) {
      Vector offset ((*hh)->point_);
      offset[0] -= px;
      offset[1] -= py;
      if (offset.norm() <= (*hh)->radius_) {
    	grab_offset = offset;
    	grabbed = *hh;
    	break;
      }
    }
  }
  else if (flags & gfx::MOUSE_DRAG) {
    if (0 != grabbed) {
      grabbed->point_[0] = px;
      grabbed->point_[1] = py;
      grabbed->point_ += grab_offset;
    }
  }
  else if (flags & gfx::MOUSE_RELEASE) {
    grabbed = 0;
  }
}


int main()
{
  struct timeval tt;
  gettimeofday (&tt, NULL);
  srand (tt.tv_usec);
  
 
  
  //Starting and ending points
  qs.resize (cdim);
  qs << -5.0, -5.0;
  qe.resize (cdim);
  qe << 7.0, 7.0;
  
  //Initialize path with all poses in Xi stacked in qs
  chomp.initStackedPath(qs, qe);  
  chomp.getPath(trajectory);

  aux_chomp.initStackedPath(qs, qe);  
  aux_chomp.getPath(trajectory2);
  
  
  update_robots();   
  state = PAUSE;
  
  gfx::add_button ("jumble", cb_jumble);
  gfx::add_button ("step", cb_step);
  gfx::add_button ("run", cb_run);
  gfx::main ("chomp", cb_idle, cb_draw, cb_mouse);
}
