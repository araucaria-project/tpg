"""Dome azimuth for a telescope on a German equatorial mount (GEM).

TEMPORARY copy of the geometry in `pyaraucaria.dome_eq.dome_eq_azimuth`. The pyaraucaria
version computes for "now" only (it calls `ra_dec_2_az_alt` / `site_sidereal_time` without a
time), while tpg needs it for every minute of a future night. Once pyaraucaria accepts a time
(or az/alt + hour angle directly), switch MaskWind/CalcObject to it and delete this file.

The math is kept identical to pyaraucaria (tests/test_dome_eq.py checks that), with two
differences that don't change results:
  * inputs are the az/alt and hour angle tpg already has from ephem, so no coordinate
    conversion happens here,
  * vectorized with numpy, and the quadrant ladder is replaced by atan2 (same angle, but no
    ZeroDivisionError when the final vector points exactly E or W).
"""

import numpy


def dome_eq_azimuth(az, alt, ha, r_dome, spx, spy, gem, side_of_pier, latitude):
    """
    @param az: telescope (target) azimuth [deg], scalar or array
    @param alt: telescope (target) altitude [deg], scalar or array
    @param ha: hour angle, sidereal time - ra [deg], scalar or array
    @param r_dome: radius of dome [mm]
    @param spx: X distance from dome center to intersection of ra - dec axes: (+) E, (-) W [mm]
    @param spy: Y distance from dome center to intersection of ra - dec axes: (+) N, (-) S [mm]
    @param gem: distance from intersection of ra - dec axes to telescope optical axis [mm]
    @param side_of_pier: on east = 0, on west = 1 (anything else is treated as 0, as in pyaraucaria)
    @param latitude: site latitude [deg]
    @return: dome azimuth [deg], 0-360, same shape as the inputs
    """
    az = numpy.radians(numpy.asarray(az, dtype=float))
    alt = numpy.radians(numpy.asarray(alt, dtype=float))
    ha = numpy.asarray(ha, dtype=float)

    if int(side_of_pier) not in [0, 1]:
        side_of_pier = 0

    # target point projected on the dome xy plane
    target_x = r_dome * numpy.cos(numpy.pi / 2 - az) * numpy.cos(alt)
    target_y = r_dome * numpy.sin(numpy.pi / 2 - az) * numpy.cos(alt)

    # gem arm, rotated with the hour angle, flipped by 180 deg on the west side of pier
    angle = numpy.radians(ha + int(side_of_pier) * 180)
    n_s = -numpy.sign(latitude)
    gem_x = numpy.cos(angle) * gem
    gem_y = n_s * numpy.sin(angle) * gem * numpy.cos(numpy.radians(90 - abs(latitude)))

    x = target_x + spx + gem_x
    y = target_y + spy + gem_y

    return numpy.degrees(numpy.arctan2(x, y)) % 360
