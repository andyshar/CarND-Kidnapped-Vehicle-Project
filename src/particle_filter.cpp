/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_t(theta, std[2]);

    num_particles = 100;

    for (int i = 0; i < num_particles; ++i) {

        Particle p = Particle();

        p.id = i;
        p.x = dist_x(dre);
        p.y = dist_y(dre);
        p.theta = dist_t(dre);
        p.weight = 1.0;

        particles.push_back(p);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    for (int i = 0; i < num_particles; ++i) {

        std::normal_distribution<double> dist_x(0, std_pos[0]);
        std::normal_distribution<double> dist_y(0, std_pos[1]);
        std::normal_distribution<double> dist_t(0, std_pos[2]);
        if (abs(yaw_rate) < 0.000001) {
            particles[i].x = particles[i].x + velocity * sin(particles[i].theta) * delta_t + dist_x(dre);
            particles[i].y = particles[i].y + velocity * cos(particles[i].theta) * delta_t + dist_y(dre);
        } else {
            particles[i].x = particles[i].x + velocity * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) / yaw_rate + dist_x(dre);
            particles[i].y = particles[i].y + velocity * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)) / yaw_rate + dist_y(dre;
            particles[i].theta = particles[i].theta + yaw_rate * delta_t + dist_t(dre);
        }

    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    for (size_t i = 0; i < predicted.size(); ++i) {

        double min_dist = numeric_limits<double>::max();
        int id = -1;

        for (size_t j = 0; j < observations.size(); ++j) {
            double distance = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
            if (distance < min_dist) {
                min_dist = distance;
                id = j;
            }
        }

        predicted[i].id = id;

    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    double gauss_norm = 1/(2 * M_PI * std_landmark[0] * std_landmark[1]);

    for (int i = 0; i < num_particles; ++i) {

        vector<LandmarkObs> predicted;
        for (size_t j = 0; j < observations.size(); ++j) {

            LandmarkObs l;
            l.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
            l.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
            predicted.push_back(l);
        }

        vector<LandmarkObs> possible;
        for (size_t j = 0; j < map_landmarks.landmark_list.size(); ++j) {
            if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) < sensor_range) {
                LandmarkObs p = LandmarkObs();
                p.x = map_landmarks.landmark_list[j].x_f;
                p.y = map_landmarks.landmark_list[j].y_f;
                p.id = map_landmarks.landmark_list[j].id_i;
                possible.push_back(p);
            }
        }

        dataAssociation(predicted, possible);

        double weight = 1.0;
        for (size_t j = 0; j < predicted.size(); ++j) {
            double exponent = pow(predicted[j].x - possible[predicted[j].id].x, 2) / (2 * std_landmark[0] * std_landmark[0]) + pow(predicted[j].y - possible[predicted[j].id].y, 2) / (2 * std_landmark[1] * std_landmark[1]);
            weight *= gauss_norm * exp(-exponent);
        }

        particles[i].weight = weight;

    }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    vector<double> weights;
    for (int i = 0; i < num_particles; ++i) {
        weights.push_back(particles[i].weight);
    }

    discrete_distribution<> dist(weights.begin(), weights.end());
    std::vector<Particle> particles_resample;
    for (size_t i = 0; i < particles.size(); ++i) {
        particles_resample.push_back(particles[dist(dre)]);
    }

    particles = particles_resample;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
