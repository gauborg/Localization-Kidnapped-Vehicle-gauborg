/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Gaurav Borgaonkar Implementation of Udacity Kidnapped Vehicle project
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

// define random engine
std::default_random_engine gen;

// include normal distribution
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles

  // std is the uncertainities in GPS data
  // get uncertainties in x, y and theta
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x_init(x, std_x);
  
  // TODO: Create normal distributions for y and theta
  normal_distribution<double> dist_y_init(y, std_y);
  normal_distribution<double> dist_theta_init(theta, std_theta);

  // initialize all particles to first position (based on GPS)
  for (int i = 0; i <= num_particles; i++)
  {
    // define a particle object
    Particle p;

    // assign properties
    p.id = i;
    p.x = x;
    p.y = y;
    p.theta = theta;

    // initialize the weights to 1.0
    p.weight = 1.0;

    // adding noise

    p.x = p.x + dist_x_init(gen);
    p.y = p.y + dist_y_init(gen);
    p.theta = p.theta + dist_theta_init(gen);

    // add created particle to the set of particles
    particles.push_back(p);

    // std::cout<<p.id<<", "<<p.x<<", "<<p.y<<std::endl;

  }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  // std_pos is GPS measurement uncertainty

  // adding noise
  normal_distribution<double> dist_x(0, std_pos[0]);
  
  // TODO: Create normal distributions for y and theta
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  // iterate through all the particles
  // for (auto it = particles.begin(); particles.end(); particles++)
  for (int i = 0; i < num_particles; i++)
  {
    // apply motion models here to calculate new position after time delta_t
    particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));

    particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));

    particles[i].theta = particles[i].theta + (yaw_rate*delta_t);

    // add random noise
    particles[i].x = particles[i].x + dist_x(gen);
    particles[i].y = particles[i].y + dist_y(gen);
    particles[i].theta = particles[i].theta + dist_theta(gen);

  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  // we compare the predictions with our observations

  for (unsigned int i = 0; i < observations.size(); i++)
  {
    // current observation
    LandmarkObs obs = observations[i];

    // init distance
    double min_dist = std::numeric_limits<double>::max();

    // init id of landmark from map placeholder to be associated with the observation
    int map_id = -1;

    for (unsigned int j = 0; j < predicted.size(); j++)
    {
      LandmarkObs pred = predicted[j];

      // use helper function to get distance between observation and prediction
      double cur_dist = dist(obs.x, obs.y, pred.x, pred.y);

      // find nearest landmark to predicted landmark
      if (cur_dist < min_dist)
      {
        min_dist = cur_dist;
        map_id = pred.id;
      }

    }
    // set the observation's id to the nearest predicted landmark's id
    observations[i].id = map_id;

  }

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  // for every particle

  for (int i = 0; i < num_particles; i++)
  {

    // get the properties of particles
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;

    vector<LandmarkObs> predictions;

    // first isolate landmarks within range of sensor(vehicle) location
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      
      double x_map = map_landmarks.landmark_list[j].x_f;
      double y_map = map_landmarks.landmark_list[j].y_f;
      // double id_map = map_landmarks.landmark_list[j].id_i;

      // first check if the landmark is within the sensor range
      if(fabs(x_map - p_x) <= sensor_range && fabs(y_map - p_y))
        predictions.push_back(LandmarkObs{x_map, y_map});
        
    } // end of j loop for map landmarks

    // vector to store transformed observed positions
    vector<LandmarkObs> transformed_obs;

    // iterate through the observations to transform them to map co-ordinates
    for (unsigned int j = 0; j < observations.size(); j++)
    {
      // transformed coordinates
      // perform transformation from vehicle co-ordinates to map co-ordinates
      double trans_x, trans_y;
      trans_x = p_x + (cos(p_theta) * observations[j].y) - (sin(p_theta) * observations[j].y);
      trans_y = p_y + (sin(p_theta) * observations[j].x) + (cos(p_theta) * observations[j].y);

      // add id and transformed co-ordinates to vector
      transformed_obs.push_back(LandmarkObs{observations[j].id, trans_x, trans_y });

    }

    // associate transformed map landmarks with predictions
    dataAssociation(predictions, transformed_obs);

    // we re-initialize weights
    particles[i].weight = 1.0;

    // iterate through the transformed observations and calculate weights
    for (unsigned int j = 0; j < transformed_obs.size(); j++)
    {
      double prediction_x, prediction_y;
      
      double mu_x = transformed_obs[j].x;
      double mu_y = transformed_obs[j].y;

      // get the x,y coordinates of the prediction associated with the current observation
      for (unsigned int k = 0; k < predictions.size(); k++)
      {
        // if prediction matches the transformed observations
        if (predictions[k].id == transformed_obs[j].id)
        {
          prediction_x = predictions[k].x;
          prediction_y = predictions[k].y;
        }
      }

      /* --- calculate weight for this observation with multivariate Gaussian --- */
      // landmark measurement uncertainty (standard deviations)
      double sigma_x = std_landmark[0];
      double sigma_y = std_landmark[1];

      // compute multi-variate Gaussian
      double exponent = (pow((prediction_x-mu_x),2)/(2*pow(sigma_x, 2))) + (pow((prediction_y-mu_y),2)/(2*pow(sigma_y, 2)));
      double particle_weight = (1/(2*M_PI*sigma_x*sigma_y)) * exp(-exponent);

      // product of this obersvation weight with total observations weight
      // calculate total probability
      particles[i].weight = particles[i].weight * particle_weight;

    }   // end of transformed observations loop

  }     // end of particle loop

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}