from astropy.coordinates import SkyCoord, AltAz, get_body
from astropy.time import Time
from astropy import units as u
from astroquery.jplhorizons import Horizons


def compute_object_attributes_old(dt, loc, ra, dec):
    observer_alt_az = AltAz(obstime=dt, location=loc)
    target_coord = SkyCoord(ra=ra, dec=dec, frame='icrs')
    target_alt_az = target_coord.transform_to(observer_alt_az)

    moon_coord = get_body("moon", Time(dt), location=loc)
    moon_alt_az = moon_coord.transform_to(observer_alt_az)

    # ID 301 is Moon
    moon_phase = Horizons(id='301', location=f"geo", epochs=Time(dt).jd).ephemerides()["illumination"][0]

    azimuth_diff = abs(moon_alt_az.az - target_alt_az.az)

    is_target_visible = target_alt_az.alt > 0 * u.deg
    is_target_above_25_degree = target_alt_az.alt > 25 * u.deg
    is_target_in_10_degrees_of_zenith = target_close_to_zenith(target_alt_az, 10 * u.deg)
    is_target_in_30_degrees_of_zenith = target_close_to_zenith(target_alt_az, 30 * u.deg)
    is_target_in_45_degrees_of_zenith = target_close_to_zenith(target_alt_az, 45 * u.deg)

    is_moon_visible = moon_alt_az.alt > 0 * u.deg
    is_moon_phase_acceptable = moon_phase < 30 or azimuth_diff > 90 * u.deg
    is_moon_in_180_degrees = azimuth_diff <= 180 * u.deg

    attributes = {
        "target_alt": target_alt_az.alt,
        "target_az": target_alt_az.az,
        "target_visible": is_target_visible,
        "is_target_above_25_degree": is_target_above_25_degree,
        "is_target_in_10_degrees_of_zenith": is_target_in_10_degrees_of_zenith,
        "is_target_in_30_degrees_of_zenith": is_target_in_30_degrees_of_zenith,
        "is_target_in_45_degrees_of_zenith": is_target_in_45_degrees_of_zenith,
        "moon_alt": moon_alt_az.alt,
        "moon_az": moon_alt_az.az,
        "moon_phase": moon_phase,
        "is_moon_visible": is_moon_visible,
        "is_moon_phase_acceptable": is_moon_phase_acceptable,
        "is_moon_in_180_degrees": is_moon_in_180_degrees
    }

    return attributes


def target_close_to_zenith(target_alt_az, threshold):
    return abs(target_alt_az.alt - 90 * u.deg) <= threshold
