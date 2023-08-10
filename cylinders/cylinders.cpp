#include <iostream>
#include <cmath>

using namespace std;

double const pi = acos(-1);

struct Point {
    double x, y, z;
};

struct Pipe {
    Point point;
    double rad;
    double alpha = 0;
};

bool input(Pipe &pipe1, Pipe &pipe2) {
    cout << "Input x0, y0, z0: ";
    cin >> pipe1.point.x >> pipe1.point.y >> pipe1.point.z;
    cout << "Input x1, y1, z1: ";
    cin >> pipe2.point.x >> pipe2.point.y >> pipe2.point.z;
    cout << "Input alpha(deg): ";
    cin >> pipe2.alpha;
    pipe2.alpha /= 180, pipe2.alpha *= pi;
    cout << "Input radii R0, R1: ";
    cin >> pipe1.rad >> pipe2.rad;
    if (pipe2.alpha == 0) return false;
    else return true;
}

Point rotate_coord_system(Point const &point, double const sine, double const cosine) {
    return {point.x * cosine + point.y * sine, -point.x * sine + point.y * cosine, point.z};
}

int sgn(double const y0, double const y1) {
    if (y1 - y0 >= 0) return 1;
    else return -1;
}

Point calc(Pipe const &pipe1, Pipe const &pipe2) {
    /// Description:
    /// if XZ plane is parallel to both of the pipes then
    /// coordinates X Z of equidistant point are the same of the intersection point of pipes axes
    /// if the axes of pipes lied in the same plane
    /// Y coord can be calculated easily as arithmetic average of  y coordinates of axes (taking into account radii of pipes)

    /// this function rotates coord system so that the condition above is valid
    /// then calculates the intersection point and converts its coords back to XYZ
    Point coord;
    double sin_a = pipe2.point.x / sqrt(pipe2.point.x * pipe2.point.x + pipe2.point.y * pipe2.point.y),
    cos_a = -pipe2.point.y  / sqrt(pipe2.point.x * pipe2.point.x + pipe2.point.y * pipe2.point.y);

    Point p_wire1 = rotate_coord_system(pipe1.point, sin_a, cos_a),
    p_wire2 = rotate_coord_system(pipe2.point, sin_a, cos_a);
    int sign = sgn(p_wire1.y, p_wire2.y);

    coord.x = p_wire1.x;
    coord.y = ((p_wire1.y + sign * pipe1.rad) + (p_wire2.y - sign * pipe2.rad)) / 2;
    coord.z = (p_wire1.x - p_wire2.x) * cos(pipe2.alpha) / sin(pipe2.alpha) + p_wire2.z;

    return rotate_coord_system(coord, -sin_a, cos_a);
}

int main()
{
    Pipe pipe1, pipe2;
    if (!input(pipe1, pipe2)) {
        cout << "The problem does not have an unambiguous solution.";
        return 0;
    }
    Point approx_p = calc(pipe1, pipe2);
    cout << "Coordinates of the equidistant point X, Y, Z: ";
    cout << approx_p.x << " " << approx_p.y << " " << approx_p.z;
    return 0;
}
