/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"
#define USE_MATH_DEFINES
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

// random number generator to use in various methods
std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  if(is_initialized) {
    return;
  }

  num_particles = 100;  // TODO: Set the number of particles

  // Initialise particles
  // Create normal distr's for x,y,theta
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++) {
    Particle p_temp;

    // Set ID to i
    p_temp.id = i;

    // Sample from the above distr's to set x,y,theta
    p_temp.x = dist_x(gen);
    p_temp.y = dist_y(gen);
    p_temp.theta = dist_theta(gen);

    // Set weight to 1
    p_temp.weight = 1.0;

    // Append the particle to the object's particles vector
    particles.push_back(p_temp);
  }
  is_initialized = true;
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

  // Create normal distributions
  std::normal_distribution<double> noise_x(0, std_pos[0]);
  std::normal_distribution<double> noise_y(0, std_pos[1]);
  std::normal_distribution<double> noise_theta(0, std_pos[2]);

  // Loop through all particles
  for(int i = 0; i < num_particles; i++) {
    Particle p = particles[i];

    double vel_yr = velocity / yaw_rate;

    // Calculate new position
    if(std::isinf(vel_yr)) { 
      // Velocity/Yaw rate is infinity
      p.x += velocity*delta_t*cos(p.theta);
      p.y += velocity*delta_t*sin(p.theta);

    } else {
      p.x += vel_yr*( sin(p.theta + yaw_rate*delta_t) - sin(p.theta) );
      p.y += vel_yr*( -cos(p.theta + yaw_rate*delta_t) + cos(p.theta) );
      p.theta += yaw_rate * delta_t;
    }    

    // Add Gaussian noise
    p.x += noise_x(gen);
    p.y += noise_y(gen);
    p.theta += noise_theta(gen);

    particles[i] = p;
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

  // Loop through observations
  for (int i = 0; i < observations.size(); i++) {
    LandmarkObs obs = observations[i];

    // Initialise matching landmark distance and ID
    double min_distance = std::numeric_limits<double>::max();

    // Loop through predictions
    for (int j = 0; j < predicted.size(); j++) {
      LandmarkObs pred = predicted[j];

      // Find the distance
      double cur_distance = dist(obs.x, obs.y, pred.x, pred.y);

      if(cur_distance < min_distance) {
        // Set the obs ID to the pred id and update min dist
        observations[i].id = predicted[j].id;
        min_distance = cur_distance;
      }
    }
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
  
  // Loop through each particle
  for (int i = 0; i < num_particles; i++) {
    double x_part = particles[i].x;
    double y_part = particles[i].y;
    double theta_part = particles[i].theta;

    vector<LandmarkObs> predictions;

    // Only consider landmarks within sensor range for predictions
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      float x_lm = map_landmarks.landmark_list[j].x_f;
      float y_lm = map_landmarks.landmark_list[j].y_f;
      int id_lm = map_landmarks.landmark_list[j].id_i;

      if (dist(x_lm, y_lm, x_part, y_part) <= sensor_range) {
        predictions.push_back(LandmarkObs{id_lm, x_lm, y_lm});
      }
    }

    // Convert observations to map coordinates using 2D transform
    vector<LandmarkObs> transformed_obs;

    for (int j = 0; j < observations.size(); j++) {
      // Convert vehicle coordinates to map coordinates
      double x_obs = observations[j].x;
      double y_obs = observations[j].y;
      
      // Apply transform
      double x_transf = x_part + x_obs*cos(theta_part) - y_obs*sin(theta_part);
      double y_transf = y_part + x_obs*sin(theta_part) + y_obs*cos(theta_part);

      // Append to transformed_obs
      transformed_obs.push_back(LandmarkObs{observations[j].id, x_transf, y_transf});
    }

    // Perform data association - sets the obs ID to the closest pred landmark ID
    dataAssociation(predictions, transformed_obs);

    // Re-initialise weight
    particles[i].weight = 1.0;

    for (int j = 0; j < transformed_obs.size(); j++) {     
      // Placeholder positions for matches
      double x_obs, y_obs, x_pred, y_pred; 
      x_obs = transformed_obs[j].x;
      y_obs = transformed_obs[j].y;
      int associated_prediction = transformed_obs[j].id;

      // Assign position values if id's match
      for (int k = 0; k < predictions.size(); k++) {
        if(predictions[k].id == associated_prediction) {
          x_pred = predictions[k].x;
          y_pred = predictions[k].y;
        }
      }

      // Standard deviations
      double std_x = std_landmark[0];
      double std_y = std_landmark[1];

      // Calculate Gaussian pdf
      double norm_term = 1/(2*M_PI*std_x*std_y);
      double exponent_term1 = -pow(x_obs - x_pred,2)/(2*pow(std_x, 2));
      double exponent_term2 = -pow(y_obs - y_pred,2)/(2*pow(std_y, 2));
      double pdf_obs =  norm_term * exp(exponent_term1 + exponent_term2);

      // Multiply to current weight
      particles[i].weight *= pdf_obs;
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  // Vector to hold the result
  vector<Particle> resampled_particles;

  // Vector of weights
  vector<double> weights; 

  for (int i = 0; i < particles.size(); i++) {
    weights.push_back(particles[i].weight);
  }

  // Find the max weight
  double max_weight = *max_element(weights.begin(), weights.end());;

  // Set a random initial index
  std::uniform_int_distribution<int> uni_int_dist(0, num_particles-1);
  int index = uni_int_dist(gen);

  // Create a random distr for beta calculation
  std::uniform_real_distribution<double> uni_real_dist(0.0, max_weight);

  double beta = 0.0;

  // Use a resampling wheel to pick particles
  for(int i = 0; i < num_particles; i++) {
    beta += uni_real_dist(gen) * 2.0;
    while(beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampled_particles.push_back(particles[index]);
  }

  particles = resampled_particles;

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