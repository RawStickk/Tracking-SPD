#include <iostream>
#include <cmath>

using namespace std;

double const pi = acos(-1);

struct Point {
    double x, y, z;
};

struct Line {
    Point point;
    double alpha = 0;
};

bool input(Line &wire1, Line &wire2) {
    cout << "Input x0, y0, z0: ";
    cin >> wire1.point.x >> wire1.point.y >> wire1.point.z;
    cout << "Input x1, y1, z1: ";
    cin >> wire2.point.x >> wire2.point.y >> wire2.point.z;
    cout << "Input alpha(deg): ";
    cin >> wire2.alpha;
    wire2.alpha /= 180;
    wire2.alpha *= pi;
    if (wire2.alpha == 0) return false;
    else return true;
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

Point rotate_coord_system(Point const &point, double const angle) {
    return {point.x * cos(angle) + point.y * sin(angle), -point.x * sin(angle) + point.y * cos(angle), point.z};
}

Point calc(Line const &wire1, Line const &wire2) {
    Point coord;
    double teta = calc_angle(wire2.point.y / sqrt(wire2.point.x * wire2.point.x + wire2.point.y * wire2.point.y),
                             wire2.point.x / sqrt(wire2.point.x * wire2.point.x + wire2.point.y * wire2.point.y));

    Point p_wire1 = rotate_coord_system(wire1.point, teta + pi/2), p_wire2 = rotate_coord_system(wire2.point, teta + pi/2);

    coord.x = p_wire1.x;
    coord.y = (p_wire1.y + p_wire2.y) / 2;
    coord.z = (p_wire1.x - p_wire2.x) * cos(wire2.alpha) / sin(wire2.alpha) + p_wire2.z;

    return rotate_coord_system(coord, -(teta + pi/2));
}

int main()
{
    Line wire1, wire2;
    if (!input(wire1, wire2)) {
        cout << "The problem does not have an unambiguous solution.";
        return 0;
    }
    Point intersection_point = calc(wire1, wire2);
    cout << "Coordinates of the equidistant point X, Y, Z: ";
    cout << intersection_point.x << " " << intersection_point.y << " " << intersection_point.z;
    return 0;
}
