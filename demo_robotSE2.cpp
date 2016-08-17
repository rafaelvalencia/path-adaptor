/* demo_robotSE2
 *
 * Copyright (C) 2016 Rafael Valencia. All rights reserved.
 * License (3-Cluase BSD):  
 *
 * This code uses and is based on code from:
 *   Project: trychomp https://github.com/poftwaresatent/trychomp
 *   Copyright (C) 2014 Roland Philippsen. All rights reserved.
 *   License (3-Clause BSD) : https://github.com/poftwaresatent/trychomp
 * **
 * \file demo_robotSE2.cpp
 *
 * Shows an extension to CHOMP to account for a planar differential drive mobile robot, 
 * including its nonholonomic restrictions. It will
 * plan a trajectory (xi) connecting start point (qs) to end point (qe) while
 * avoiding obstacles (obs)
 */
 


#include "gfx.hpp"
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdlib.h>
#include <sys/time.h>
#include <err.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "path_adaptor.hpp"

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::Isometry3d Transform;

static size_t const obs_dim(3); 	// obstacle dimensions (x,y,Radius)

using namespace std;

//ofstream myfile("chompResults.txt");

//////////////////////////////////////////////////
// trajectory etc
Vector trajectory;				// the trajectory (q_1, q_2, ...q_n)
Matrix  obs(obs_dim,2);   		// Matrix containninging all obstacles, each column is (x,y,R) of the obstacle
Vector qs;						// the start config a.k.a. q_0
Vector qe;						// the end config a.k.a. q_(n+1)
static size_t const cdim (3);	// dimension of config space
static size_t const nq (20);	// number of q stacked into xi
static size_t const xidim (nq * cdim); // dimension of trajectory, xidim = nq * cdim

//////////////////////////////////////////////////
// CHOMP
CHOMP_SE2 chomp(1.0, 100.0, 1.0, 20, 1000, 10.0);

/* CHOMP parameters
dt     = 1.0 	(time step)
eta    = 100.0 	(regularization factor for gradient descent)
lambda = 1.0	(weight of smoothness objective)
nq     = 20		(number of poses q in xi)
cdim   = 3		(dimension of config space)
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
    : position_ (Vector::Zero(2)),
    heading_(0)
  {
  }
  
  
  void update (Vector const & pose)
  {
    if (pose.size() != 3) {
      errx (EXIT_FAILURE, "Robot::update(): position has %zu DOF (but needs 3)",
	    (size_t) pose.size());
    }
    position_ = pose.block(0,0,2,1);
    heading_ = pose(2);
  }
  
  
  void draw () const
  {
    // translucent disk for base
    gfx::set_pen (1.0, 0.7, 0.7, 0.7, 0.5);
    gfx::fill_arc (position_[0], position_[1], radius_, 0.0, 2.0 * M_PI);
    
    // thick circle outline for base
    gfx::set_pen (3.0, 0.2, 0.2, 0.2, 1.0);
    gfx::draw_arc (position_[0], position_[1], radius_, 0.0, 2.0 * M_PI);
    
    // thick line for heading
    gfx::set_pen (3.0, 0.2, 0.2, 0.2, 1.0);
    gfx::draw_line (position_[0], position_[1], position_[0] + 0.5 * cos(heading_), position_[1] + 0.5* sin(heading_));
  }
  
  
 
  
  static double const radius_;
  
  Vector position_;
  double heading_;
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
  if (nq != robots.size()) {
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
	}
	//update new initial path
	chomp.setPath(qs, qe, trajectory);
	update_robots();
}


static void cb_idle ()
{
  if (PAUSE == state) {
    return;
  }
  if (STEP == state) {
    state = PAUSE;
  }

  static size_t stepcounter (0);
  //cout << "step " << stepcounter++ << "\n";
  
  /////////////////////////////////////////////////
  // Get updated obstacles  
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
  //////////////////////////////////////////////////
  // One CHOMP iteration  
  chomp.setObstacles(obs); //update obstacles  
  chomp.chompIteration(trajectory); //update path
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
  //Vector qqs, qqe;
  qs.resize (3);
  qs << -5.0, -5.0, -1.57;
  qe.resize (3);
  qe << 7.0, 7.0, 1.57;
 
  //Initialize path with all poses in Xi stacked in qs
  chomp.initStackedPath(qs, qe);   
  chomp.getPath(trajectory);
  update_robots(); 

  state = PAUSE;
  
  gfx::add_button ("jumble", cb_jumble);
  gfx::add_button ("step", cb_step);
  gfx::add_button ("run", cb_run);
  gfx::main ("chomp", cb_idle, cb_draw, cb_mouse);   
}
