#include "../include/disk.h"

IsInRange Disk::isInRange(double tlat, double tlon, double r_b)
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
        return IsInRange::Inside;
    }
    else if (dist <= this->roi)
    {
        return IsInRange::InsideInfluence;
    }
    else
    {
        return IsInRange::Outside;
    }
}

IsInRange mergeResults(IsInRange *enums, size_t N)
{
    bool isInsideInfluence = 0;
    for (size_t i = 0; i < N; i++)
    {
        IsInRange temp = enums[i];
        if (temp == 1)
        {
            /*
                If some enum reports Inside, return Inside regardless of
                the other results
            */
            return IsInRange::Inside;
        }
        else if (temp == 2)
        {
            // if some enum reports InsideInfluence switch flag
            isInsideInfluence = 1;
        }
    }
    if (isInsideInfluence)
    {
        return IsInRange::InsideInfluence;
    }
    else
    {
        return IsInRange::Outside;
    }
}