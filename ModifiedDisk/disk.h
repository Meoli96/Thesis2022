#include "pch.h"

#include "utils.h"

typedef enum
{
    Outside = 0,
    Inside,
    InsideInfluence
} IsInRange;

class Point
{
private:
    // Stored in radians
    double lat;
    double lon;

public:
    Point(double a, double b) : lat(degree2radians(a)), lon(degree2radians(b)) {}
    ~Point() {}
    double getLat() { return this->lat; }
    double getLon() { return this->lon; }
};

class Disk
{
private:
    Point c;
    double radius; // km
    double roi;    // km

public:
    Disk(Point p, double r) : c(p), radius(r), roi(1.2 * r) {}
    ~Disk() {}
    IsInRange isInRange(double lat, double lon, double r);
};

IsInRange mergeResults(IsInRange *, size_t);