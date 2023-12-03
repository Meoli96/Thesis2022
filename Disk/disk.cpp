#include "../include/disk.h"


bool Disk::isInRange(double tlat, double tlon, double r_b)
{

    
    bool res = 0;
    double dlat = this->c.getLat();
    double dlon = this->c.getLon();
    double r_tlat = degree2radians(tlat);
    double r_tlon = degree2radians(tlon);

    double sinlat = sin(r_tlat - dlat);
    double sinlon = sin(r_tlon - dlon);

    double dist = 2 * r_b * asin(sqrt((sinlat * sinlat / 2 + cos(dlat) * cos(r_tlat) * (sinlon * sinlon / 2))));
  

    if (dist <= this->radius)
    {
        // Body is in range
        return true;
    }
    else
    {
        return false;
    }
}