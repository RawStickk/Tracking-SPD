#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

using namespace std;

double const pi = acos(-1);
double const e = 1.602 * 1e-19;
double const c = 2.998 * 1e8;

///all magnitudes of geometrical parameters are signed in mm

struct Point {
    double x, y, z;
};

struct Wire {
    Point point = {0, 0, 0};
    double alpha = 0;
};

struct Particle {
    Point p {0, 0, 0}; /// in [MeV/c]
    double charge = 1; /// in [e]
};

struct Hit {
    Wire pipe;
    double R;
    unsigned ntrk;
    Point p;
};

void construct(vector<Wire> *pipes, double const wire_rad, double const wire_angle,
               double const R0, double const d,
               unsigned const N) {

    double R = R0, angle = 0;

    for (unsigned i = 0; i < N; i++) {
        double l = 0;
        while (2 * pi * R - l >= wire_rad * 2) {
            pipes[i].push_back({{R * cos(l/R), -R * sin(l/R), 0}, angle});
            l += 2 * wire_rad;
        }
        if (i % 2 == 0) {
            R += 2 * wire_rad;
            angle = wire_angle;
        }
        else {
            R += d;
            angle = 0;
        }
    }
}

void pipes_output(vector<Wire> const *pipes, unsigned const N) {
    fstream fout;
    fout.open("wires.csv", ios::out | ios::trunc);

    fout << "x" << "," << "y" << "," << "alpha" << "\n";
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < pipes[i].size(); j++) {
            fout << pipes[i][j].point.x << "," << pipes[i][j].point.y << "," << pipes[i][j].alpha << "\n";
        }
    }

    fout.close();
}

void particle_parameters(Particle &particle, unsigned const num_particles, unsigned const seed) {
    default_random_engine rng(seed);
    uniform_real_distribution<double> dstr(-500, 500);
    particle.p.x = dstr(rng);
    particle.p.y = dstr(rng);
    particle.p.z = dstr(rng);
}

Point rotate_coord_system(Point const &point, double const angle) {
    return {point.x * cos(angle) + point.y * sin(angle), -point.x * sin(angle) + point.y * cos(angle), point.z};
}

double modulus_(Point const &point) {
    return point.x * point.x + point.y * point.y + point.z * point.z;
}

Point vect_prod(Point const &a, Point const &b) {
    return {a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y};
}

Point dot(double const c, Point const p) {
    return {c * p.x, c * p.y, c * p.z};
}

double calc_angle(double const sine, double const cosine) {
    double angle;
    if (abs(cosine) < 1e-9) {
        angle = asin(sine);
    }
    else {
        angle = atan(sine / cosine);
        if (cosine < 0) angle += pi;
    }
    return angle;
}

double dist(Point const &track_point, Wire const &pipe) {
    double teta = calc_angle(pipe.point.y / sqrt(pipe.point.x * pipe.point.x + pipe.point.y * pipe.point.y),
                             pipe.point.x / sqrt(pipe.point.x * pipe.point.x + pipe.point.y * pipe.point.y));

    Point l = {sin(pipe.alpha * pi / 180), 0, cos(pipe.alpha * pi / 180)},
    a = {track_point.x - pipe.point.x, track_point.y - pipe.point.y, track_point.z - pipe.point.z};
    a = rotate_coord_system(a, teta + pi/2);

    return sqrt(modulus_(vect_prod(l, a)));
}

void hits_output(vector<Hit> const &hits) {
    fstream fout;
    fout.open("hits.csv", ios::out | ios::trunc);
    fout << "x,y,R,alpha,NTRK,p_x,p_y,p_z" << "\n";
    for (unsigned i = 0; i < hits.size(); i++) {
        fout << hits[i].pipe.point.x << "," << hits[i].pipe.point.y << "," << hits[i].R << "," << hits[i].pipe.alpha << ",";
        fout << hits[i].ntrk << "," << hits[i].p.x << "," << hits[i].p.y << "," << hits[i].p.z << "\n";
    }
    fout.close();
}

void hit_detection(vector<Wire> const *pipes, unsigned const N, double const len, double const pipe_rad,
        Point const &B, Particle const &particle, unsigned const ntrk, vector<Hit> &hits) {
    Point p_xy = {particle.p.x, particle.p.y, 0};
    Point rho = dot(c * 1e-3 * modulus_(particle.p) / (particle.charge * modulus_(B) * modulus_(p_xy)), vect_prod(particle.p, B));
    double beta = calc_angle(particle.p.y / sqrt(modulus_(p_xy)), particle.p.x / sqrt(modulus_(p_xy)));

    for (unsigned i = 0; i < N; i++) {
        double R_i = sqrt(modulus_({pipes[i][0].point.x, pipes[i][0].point.y, 0})) - pipe_rad;

        if (R_i * R_i / modulus_(rho) > 4) break;

        double phi = acos(1 - R_i * R_i / (2 * modulus_(rho))), phi0;
        double *gap = new double;
        *gap = sqrt((len * pipes[i][0].alpha * pi / 180) * (len * pipes[i][0].alpha * pi / 180) + (R_i + pipe_rad) * (R_i + pipe_rad))
        - R_i + pipe_rad;

        if ((R_i + *gap) * (R_i + *gap) / modulus_(rho) > 4)
            phi0 = pi;
        else
            phi0 = acos(1 - (R_i + *gap) * (R_i + *gap) / (2 * modulus_(rho)));
        delete gap;

        double const d_phi = (phi0 - phi) / 100;

        bool hit = false, out = false;
        Wire wire;
        double min_dist;

        while (phi < phi0) {
            Point r = {sqrt(modulus_(rho)) * sin(phi - beta) + rho.x, sqrt(modulus_(rho)) * cos(phi - beta) + rho.y,
            sqrt(modulus_(rho)) * particle.p.z * phi / sqrt(modulus_(p_xy))};
            fout << r.x << "," << r.y << "," << r.z << "\n";
            if (r.z > len) {
                cout << i << " out \n";
                out = true;
                break;
            }

            for (unsigned j = 0; j < pipes[i].size() && !hit; j++) {
                double s = dist(r, pipes[i][j]);
                if (dist(r, pipes[i][j]) < pipe_rad) {
                    hit = true;
                    cout << i << " hit \n";
                    min_dist = dist(r, pipes[i][j]);
                    wire = pipes[i][j];
                }
            }
            if (dist(r, wire) < min_dist)
                min_dist = dist(r, wire);

            phi += d_phi;
        }
        if (out) break;
        hits.push_back({wire, min_dist, ntrk, particle.p});
    }
}

int main() {
    double const R0 = 100, d = 50,
    pipe_rad = 5, len = 1000, wire_angle = 5;
    unsigned const N = 40;
    Point const B = {0, 0, 1}; /// in [T]

    vector<Wire> pipes[N];
    construct(pipes, pipe_rad, wire_angle, R0, d, N);
    pipes_output(pipes, N);

    unsigned seed = 1111;
    default_random_engine rng(seed);
    uniform_int_distribution<unsigned> dstr(1, 10);
    uniform_int_distribution<unsigned> seed_dstr(1, 1000);

    unsigned num_particles = dstr(rng);
    cout << num_particles << endl;
    Particle *particles = new Particle[num_particles];
    vector<Hit> hits;
    for (unsigned i = 0; i < num_particles; i++) {
        particle_parameters(particles[i], num_particles, seed_dstr(rng));
        cout << particles[i].p.x << " " << particles[i].p.y << " " << particles[i].p.z << "\n";
        hit_detection(pipes, N, len, pipe_rad, B, particles[i], i + 1, hits);
    }
    //shuffle(hits.begin(), hits.end(), default_random_engine(seed));
    hits_output(hits);
    delete[] particles;
    return 0;
}
