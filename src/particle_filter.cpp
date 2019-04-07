/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <random>
#include <iostream>

using std::string;
using std::vector;
using std::normal_distribution;

/**
 * Set the number of particles. Initialize all particles to
 *   first position (based on estimates of x, y, theta and their uncertainties
 *   from GPS) and all weights to 1. And add random Gaussian noise to each particle.
 * NOTE: Consult particle_filter.h for more information about this method
 *   (and others in this file).
 * @param x  sense_x
 * @param y  sense_y
 * @param theta theta
 * @param std std array
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  std::default_random_engine r;
  num_particles = 99;  // TODO: Set the number of particles


  // Create a normal distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i) {
    Particle p;

    p.id = i;
    p.x = dist_x(r);
    p.y = dist_y(r);
    p.theta = dist_theta(r);
    p.weight = 1.0;

    particles.push_back(p);
  }

  is_initialized = true;

}

/**
 * Add measurements to each particle and add random Gaussian noise.
 * @param delta_t
 * @param std_pos
 * @param velocity
 * @param yaw_rate
 *
 * When adding noise you may find std::normal_distribution
 *   and std::default_random_engine useful.
 *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
 *  http://www.cplusplus.com/reference/random/default_random_engine/
 */
void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  std::default_random_engine r;

  // Create a normal (Gaussian) distribution for sensor noise
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; i++) {
    // predict next position x and sy
    if (fabs(yaw_rate) > 0.0001) {
      particles[i].x += (velocity / yaw_rate)
          * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      particles[i].y += (velocity / yaw_rate)
          * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      particles[i].theta += yaw_rate * delta_t;
    } else {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }

    // Add random Gaussian noise
    particles[i].x += dist_x(r);
    particles[i].y += dist_y(r);
    particles[i].theta += dist_theta(r);
  }

}

/**
 * Find the predicted measurement that is closest to each
 *   observed measurement and assign the observed measurement to this
 *   particular landmark.
 * NOTE: this method will NOT be called by the grading code. But you will
 *   probably find it useful to implement this method and use it as a helper
 *   during the updateWeights phase.
 * @param predicted
 * @param observations
 */
void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs> &observations) {
  for (unsigned int i = 0; i < observations.size(); i++) {
    int closest_landmark = 0;
    double min_dist = std::numeric_limits<float>::max();

    // Iterate through all predicted landmarks to check which is closest
    for (unsigned int j = 0; j < predicted.size(); j++) {
      // Calculate EuclideanDistance
      double curr_dist = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
      // Compare to min_dist and update if closest
      if (curr_dist < min_dist) {
        min_dist = curr_dist;
        closest_landmark = predicted[j].id;
      }
    } // end for loop of predicted
    observations[i].id = closest_landmark;
  }// end for loop of observations
}

/**
 * Update the weights of each particle using a mult-variate Gaussian
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
 * @param sensor_range range of sensor
 * @param std_landmark landmark std
 * @param observations observations
 * @param map_landmarks landmarks of map
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  // Update all particle weights
  for (int i = 0; i < num_particles; i++) {

    // Find landmark
    vector<LandmarkObs> landmarks;
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      float landmark_x = map_landmarks.landmark_list[j].x_f;
      float landmark_y = map_landmarks.landmark_list[j].y_f;

      // Only update dist small than sensor range
      if (dist(landmark_x, landmark_y, particles[i].x, particles[i].y) <= sensor_range) {
        landmarks.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i, landmark_x,
                                        landmark_y});
      }
    }

    // Transform observations from vehicle coordinates to map coordinates
    vector<LandmarkObs> mapObservations;
    for (unsigned int j = 0; j < observations.size(); j++) {
      LandmarkObs observation;
      observation.id = observations[j].id;

      // Transform to map coordinate
      observation.x = particles[i].x + (cos(particles[i].theta) * observations[j].x)
          - (sin(particles[i].theta) * observations[j].y);
      observation.y = particles[i].y + (sin(particles[i].theta) * observations[j].x)
          + (cos(particles[i].theta) * observations[j].y);

      mapObservations.push_back(observation);
    }

    dataAssociation(landmarks, mapObservations);

    // Reset weight to 1
    particles[i].weight = 1.0;

    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];

    // Calculate normalization term
    double gaussian_norm = 1 / (2 * M_PI * sig_x * sig_y);

    double x, y;
    for (unsigned int j = 0; j < mapObservations.size(); j++) {
      unsigned int k = 0;
      do {
        if (landmarks[k].id == mapObservations[j].id) {
          x = landmarks[k].x;
          y = landmarks[k].y;
          break;
        }
        k++;
      } while (k < landmarks.size());

      double obs_x = mapObservations[j].x;
      double obs_y = mapObservations[j].y;

      // exponent
      double exponent =
          (pow(obs_x - x, 2) / (2 * pow(sig_x, 2))) + (pow(obs_y - y, 2) / (2 * pow(sig_y, 2)));

      // weight using normalization terms and exponent
      double weight = gaussian_norm * exp(-exponent);

      // update weight
      particles[i].weight *= weight;
    }

  }

}

/**
 * Resample particles with replacement with probability proportional
 *   to their weight.
 * NOTE: You may find std::discrete_distribution helpful here.
 *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
 */
void ParticleFilter::resample() {

  std::default_random_engine r;

  vector<Particle> resample_particles;

  vector<double> weights;
  double max_weight = std::numeric_limits<double>::min();
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if (particles[i].weight > max_weight) {
      max_weight = particles[i].weight;
    }
  }

  std::uniform_real_distribution<double> uniformRealDist(0, max_weight);
  std::uniform_int_distribution<int> uniformIntDist(0, num_particles - 1);
  int index = uniformIntDist(r);

  double beta = 0.0;
  double mw = *std::max_element(std::begin(weights), std::end(weights));

  for (int i = 0; i < num_particles; i++) {
    beta += uniformRealDist(r) * 2.0 * mw;

    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }

    resample_particles.push_back(particles[index]);
  }

  particles = resample_particles;

  if (IS_DEBUG) {
    std::cout << "resampled particle number: " << resample_particles.size() << std::endl;
  }

}

/**
 * SetAssociations
 * @param particle the particle to which assign each listed association, and association's (x,y) world coordinates mapping
 * @param associations The landmark id that goes along with each listed association
 * @param sense_x the associations x mapping already converted to world coordinates
 * @param sense_y the associations y mapping already converted to world coordinates
 */
void ParticleFilter::SetAssociations(Particle &particle,
                                     const vector<int> &associations,
                                     const vector<double> &sense_x,
                                     const vector<double> &sense_y) {
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

/**
 * Get current associations
 * @param best
 * @return
 */
string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

/**
 * Get sensor coordnate
 * @param best
 * @param coord
 * @return
 */
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
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}