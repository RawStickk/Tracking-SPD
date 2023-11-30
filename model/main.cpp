#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>

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
    double R; ///min dist from particle to wire
    unsigned ntrk; /// number of track
    Point p; ///particle momentum
};

void construct(vector<Wire> *pipes, double const wire_rad, double const wire_angle,
               double const R0, double const d,
               unsigned const N) {
    /// Description: constructs the detector by placing wires in layers
    /// Variables:  pipes - array of wires
    ///             wire_rad - radius of a pipe
    ///             wire_angle - incline angle
    ///             R0 - radius of the inner layer (of smallest radius)
    ///             d - distance between inclined layer and next parallel one)
    ///             N - number of layers

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
    /// Description: outputs wires coordinates into csv file

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
    /// Description: generates randomly particle momentum components

    default_random_engine rng(seed);
    uniform_real_distribution<double> dstr(-500, 500);
    particle.p.x = -2;
    particle.p.y = 0;
    particle.p.z = 0;
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
    ///Description: multiplies vector by scalar

    return {c * p.x, c * p.y, c * p.z};
}

double calc_angle(double const sine, double const cosine) {
    ///Description: finds angle value using it's sine and cosine values

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
    ///Description: finds distance from a point to a line
    ///by rotating coord system and calculating area of parallelogram

    double teta = calc_angle(pipe.point.y / sqrt(pipe.point.x * pipe.point.x + pipe.point.y * pipe.point.y),
                             pipe.point.x / sqrt(pipe.point.x * pipe.point.x + pipe.point.y * pipe.point.y));

    Point l = {sin(pipe.alpha * pi / 180), 0, cos(pipe.alpha * pi / 180)},
    a = {track_point.x - pipe.point.x, track_point.y - pipe.point.y, track_point.z - pipe.point.z};
    a = rotate_coord_system(a, teta + pi/2);

    return sqrt(modulus_(vect_prod(l, a)));
}

void hits_output(vector<Hit> const &hits) {
    ///Description: outputs all the data characterizing hits
    ///(num of track, wire coordinates and inclination angle, min dist between particle and wire, particle momentum)

    fstream fout;
    fout.open("hits.csv", ios::out | ios::trunc);
    fout << "NTRK,x,y,R,alpha,p_x,p_y,p_z" << "\n";
    for (unsigned i = 0; i < hits.size(); i++) {
        fout << hits[i].ntrk << "," << hits[i].pipe.point.x << "," << hits[i].pipe.point.y << "," << hits[i].R << "," << hits[i].pipe.alpha << ",";
        fout << hits[i].p.x << "," << hits[i].p.y << "," << hits[i].p.z << "\n";
    }
    fout.close();
}

void get_track(vector<Wire> const *pipes, unsigned const N, double const len, double const pipe_rad,
        Point const &B, Particle const &particle, unsigned const ntrk) {

    fstream fout;
    string file = "track" + to_string(ntrk) + ".csv";
    fout.open(file, ios::out | ios::trunc);

    Point p_xy = {particle.p.x, particle.p.y, 0};
    Point rho = dot(1e12/ (c * particle.charge * modulus_(B)), vect_prod(particle.p, B));
    cout << "rho: " << sqrt(modulus_(rho)) << endl;
    double beta = calc_angle(particle.p.y / sqrt(modulus_(p_xy)), particle.p.x / sqrt(modulus_(p_xy)));
    double R = sqrt(modulus_({pipes[N- 2][0].point.x, pipes[ N -2][0].point.y, 0})) - pipe_rad;

    double phi = 0, phi0 = acos(1 - R * R / (2 * modulus_(rho))), d_phi = phi0 / 500;

    while (phi <= phi0) {
        Point r = {sqrt(modulus_(rho)) * sin(phi - beta) + rho.x, sqrt(modulus_(rho)) * cos(phi - beta) + rho.y,
        sqrt(modulus_(rho)) * particle.p.z * phi / sqrt(modulus_(p_xy))}; ///current position of the particle in XYZ

        if (abs(r.z) > len) break;

        fout << "NTRK,x,y,z" << "\n";
        fout << ntrk << "," << r.x << "," << r.y << "," << r.z << "\n";

        phi += d_phi;
    }
    fout.close();
}

void hit_detection(vector<Wire> const *pipes, unsigned const N, double const len, double const pipe_rad,
        Point const &B, Particle const &particle, unsigned const ntrk, vector<Hit> &hits) {
    ///Description: detects hits for a particle
    ///             1.iterates through layers and extrapolates particle track in their vicinity
    ///             2.particle track is parametrized by azimuthal angle phi
    ///             3.iterates through this angle, particle coordinates transfers to XYZ system
    ///             and checks for hits by calc of dist between particle and wire
    ///             (if it is less than pipe_rad -> hit)

    get_track(pipes, N, len, pipe_rad, B, particle, ntrk);

    Point p_xy = {particle.p.x, particle.p.y, 0}; ///perpendicular component of particle momentum to Z axis
    Point rho = dot(1e12/ (c * particle.charge * modulus_(B)), vect_prod(particle.p, B)); ///helix radius
    double beta = calc_angle(particle.p.y / sqrt(modulus_(p_xy)), particle.p.x / sqrt(modulus_(p_xy)));
    ///iteration through layers
    for (unsigned i = 0; i < N; i++) {
        double R_i = sqrt(modulus_({pipes[i][0].point.x, pipes[i][0].point.y, 0})) - pipe_rad; ///radius of "i" layer in detector

        if (R_i * R_i / modulus_(rho) > 4) break; /// check whether a particle can reach this layer
                                                  ///(in case of strong magnetic field or small momentum its track can be curved)

        double phi = acos(1 - R_i * R_i / (2 * modulus_(rho))), phi0;
        ///inclined layers of small radius are not cylindrical near their ends even approximately
        ///[double gap] is the gap between the end of a pipe and the inner layer
        ///estimation of the gap and calculation of the ending value(phi0) of phi
        double *gap = new double;
        *gap = sqrt((len * pipes[i][0].alpha * pi / 180) * (len * pipes[i][0].alpha * pi / 180) + (R_i + pipe_rad) * (R_i + pipe_rad))
        - R_i + pipe_rad;

        if ((R_i + *gap) * (R_i + *gap) / modulus_(rho) > 4)
            phi0 = pi;
        else
            phi0 = acos(1 - (R_i + *gap) * (R_i + *gap) / (2 * modulus_(rho)));
        delete gap;

        double const d_phi = (phi0 - phi) / 100; ///step
        ///in this part extrapolation process starts
        bool hit = false, out = false;
        Wire wire;
        double min_dist;

        while (phi < phi0) {
            Point r = {sqrt(modulus_(rho)) * sin(phi - beta) + rho.x, sqrt(modulus_(rho)) * cos(phi - beta) + rho.y,
            sqrt(modulus_(rho)) * particle.p.z * phi / sqrt(modulus_(p_xy))}; ///current position of the particle in XYZ
            if (abs(r.z) > len) { ///check whether particle has flown out through the end of detector
                cout << i << " out \n";
                out = true;
                break;
            }
            ///iteration through all the pipes in the layer till hit is not detected
            ///by checking if the distance between particle ad wire is less than the pipe radius
            for (unsigned j = 0; j < pipes[i].size() && !hit; j++) {
                if (dist(r, pipes[i][j]) < pipe_rad) {
                    hit = true;
                    cout << i << " hit \n";
                    min_dist = dist(r, pipes[i][j]);
                    wire = pipes[i][j];
                }
            }
            if (dist(r, wire) < min_dist) ///calculation of min dist between particle and wire
                min_dist = dist(r, wire);

            phi += d_phi;
        }
        if (out) break;
        hits.push_back({wire, min_dist, ntrk, particle.p});
    }
}

int main() {
    double const R0 = 100, d = 50,
    pipe_rad = 5, len = 1000, wire_angle = 5; ///angle is i deg
    unsigned const N = 40; ///number of layers
    Point const B = {0, 0, 1}; /// in [T]

    vector<Wire> pipes[N];
    construct(pipes, pipe_rad, wire_angle, R0, d, N);
    pipes_output(pipes, N);

    unsigned seed = 1111;
    default_random_engine rng(seed);
    uniform_int_distribution<unsigned> dstr(1, 10);
    uniform_int_distribution<unsigned> seed_dstr(1, 1000);

    unsigned num_particles = 1;
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
