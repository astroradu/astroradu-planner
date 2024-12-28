from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_moon, get_body
from astropy.time import Time
from astropy import units as u
from datetime import datetime, time, timedelta
import ephem
import pytz
import json
import os

from astroquery.jplhorizons import Horizons
from timezonefinder import TimezoneFinder

from utils.interval_type import IntervalType

lat = 44.4268
lng = 26.1025
elev = 85
location = EarthLocation(lat=lat * u.deg, lon=lng * u.deg, height=elev * u.m)

# For later
startDate = datetime(2024, 12, 26)
endDate = datetime(2025, 3, 25)

currentDay = datetime(2024, 12, 26)
currentDateTime = Time(currentDay)

datetime_format = '%Y-%m-%d %H:%M:%S'


def read_constellations(file_name):
    file_path = os.path.join(os.path.dirname(__file__), file_name)
    with open(file_path, 'r') as file:
        json_data = json.load(file)
    return tuple(json_data.items())


def get_location_timezone(latitude, longitude):
    tz_finder = TimezoneFinder()
    timezone_str = tz_finder.timezone_at(lat=latitude, lng=longitude)
    if timezone_str:
        return pytz.timezone(timezone_str)
    else:
        raise ValueError(f"Timezone could not be determined for coordinates: ({latitude}, {longitude})")


def convert_to_utc(date_time, timezone):
    current_time = date_time.astimezone(pytz.timezone(timezone.zone)).strftime(datetime_format)
    utc_time = date_time.astimezone(pytz.utc).strftime(datetime_format)
    print("===============================")
    print("Converting to UTC:")
    print("Input time: " + current_time)
    print("UTC time: " + utc_time)
    return utc_time


# CHECKED
def get_utc_midnight_for_tonight(timezone):
    tz = pytz.timezone(timezone.zone)

    today = datetime.now(tz)
    midnight = tz.localize(datetime.combine(today, time(0, 0)), is_dst=None)
    midnight_utc = convert_to_utc(midnight, timezone)
    print("===============================")
    print("Computed UTC Midnight: " + midnight_utc)
    return midnight_utc


# noinspection PyUnresolvedReferences
def calculate_astronomical_night(latitude, longitude, elevation, timezone):
    midnight = get_utc_midnight_for_tonight(timezone)
    observer = ephem.Observer()
    observer.date = ephem.Date(midnight)
    observer.lat = str(latitude)
    observer.lon = str(longitude)
    observer.elev = elevation
    observer.horizon = '0'
    observer.pressure = 0

    sunset = observer.previous_setting(ephem.Sun())  # Sunset
    sunrise = observer.next_rising(ephem.Sun())  # Sunrise

    print("===============================")
    print("Computed Upcoming Sunset: ")
    print(ephem.localtime(sunset))
    print("Computed Upcoming Sunrise: ")
    print(ephem.localtime(sunrise))

    # -6=civil twilight, -12=nautical, -18=astronomical
    observer.horizon = '-18'

    astronomical_sunset = observer.previous_setting(ephem.Sun(), use_center=True)
    astronomical_sunrise = observer.next_rising(ephem.Sun(), use_center=True)

    astronomical_sunset_utc = ephem.to_timezone(astronomical_sunset, pytz.utc)
    astronomical_sunrise_utc = ephem.to_timezone(astronomical_sunrise, pytz.utc)

    astronomical_sunset_local = ephem.localtime(astronomical_sunset)
    astronomical_sunrise_local = ephem.localtime(astronomical_sunrise)

    print("===============================")
    print("Computed Upcoming Astronomical Sunset: ")
    print(astronomical_sunset_local)
    print("Computed Upcoming Astronomical Sunrise: ")
    print(astronomical_sunrise_local)

    return astronomical_sunset_utc, astronomical_sunrise_utc


def round_to_nearest(time, mins, round_up=True):
    computed_mins = min(mins, 60) if round_up else mins
    total_seconds = (time - time.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
    remainder = total_seconds % (computed_mins * 60)

    if round_up:
        rounded_seconds = total_seconds + (computed_mins * 60 - remainder) if remainder != 0 else total_seconds
    else:
        rounded_seconds = total_seconds - remainder

    rounded_time = time.replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(seconds=rounded_seconds)
    return rounded_time


def generate_time_interval(start_date, end_date, mins):
    rounded_start_date = round_to_nearest(start_date, mins, round_up=True)
    rounded_end_date = round_to_nearest(end_date, mins, round_up=False)

    interval_in_minutes = (rounded_end_date - rounded_start_date).total_seconds() / 60
    if interval_in_minutes < mins:
        raise ValueError("The interval between rounded start and end dates is less than the specified minutes.")

    if interval_in_minutes / mins > 48:
        raise ValueError("The interval is too large to compute.")

    print("Rounded start date: ", rounded_start_date)
    print("Rounded end date: ", rounded_end_date)

    interval_list = []
    current_time = rounded_start_date
    while current_time <= end_date:
        interval_list.append(current_time)
        current_time += timedelta(minutes=mins)

    return interval_list


def check_star_visibility(dt, star_coord):
    observer_alt_az = AltAz(obstime=dt,
                            location=EarthLocation(lat=44.4268 * u.deg, lon=26.1025 * u.deg,
                                                   height=85 * u.m))
    star_alt_az = star_coord.transform_to(observer_alt_az)
    star_visible = star_alt_az.alt > 0 * u.deg
    return star_visible


def check_star_visibility_in_time_intervals(brightest_star, datetime_list, timezone):
    star_coord = SkyCoord.from_name(brightest_star)
    results = []
    for dt in datetime_list:
        star_visible = check_star_visibility(dt, star_coord)
        result = f"At datetime: {dt.astimezone(tz=timezone)} - " \
                 f"Star: {brightest_star} was {'visible' if star_visible else 'NOT visible'}."
        results.append(result)
    return results


def clear_file(f):
    with open(f, "w") as file:
        file.write("")
    return file

def target_close_to_zenith(target_alt_az, threshold):
    return abs(target_alt_az.alt - 90 * u.deg) <= threshold


def check_ra_dec_visibility(dt, loc, ra, dec):
    observer_alt_az = AltAz(obstime=dt, location=loc)
    target_coord = SkyCoord(ra=ra, dec=dec, frame='icrs')
    target_alt_az = target_coord.transform_to(observer_alt_az)
    zenith_threshold = 30 * u.deg

    moon_coord = get_body("moon", Time(dt), location=loc)
    moon_alt_az = moon_coord.transform_to(observer_alt_az)

    # ID 301 is Moon
    moon_phase = Horizons(id='301', location=f"geo", epochs=Time(dt).jd).ephemerides()["illumination"][0]

    azimuth_diff = abs(moon_alt_az.az - target_alt_az.az)

    is_target_visible = target_alt_az.alt > 0 * u.deg
    is_target_above_25_degree = target_alt_az.alt > 25 * u.deg
    is_moon_dim = moon_phase < 30
    is_moon_far = azimuth_diff > 90 * u.deg
    is_close_to_zenith = target_close_to_zenith(target_alt_az, zenith_threshold)

    return {
        "target_visible": is_target_visible,
        "is_target_above_25_degree": is_target_above_25_degree,
        "is_target_close_to_zenith": is_close_to_zenith,
        "target_alt": target_alt_az.alt,
        "target_az": target_alt_az.az,
        "moon_alt": moon_alt_az.alt,
        "moon_az": moon_alt_az.az,
        "moon_phase": moon_phase,
        "moon_brightness_below_30": is_moon_dim,
        "moon_farther_than_90_deg_az": is_moon_far,
    }


def compute_bright_star_visibility():
    timezone = get_location_timezone(lat, lng)
    astronomical_sunset, astronomical_sunrise = calculate_astronomical_night(lat, lng, elev, timezone)
    constellations = read_constellations("constellations.json")
    output_file = "visibility_results.txt"

    try:

        print("----")

        mins = IntervalType.QUARTER

        datetime_list = generate_time_interval(astronomical_sunset, astronomical_sunrise, mins.value)

        with open(output_file, "w") as file:
            for constellation, brightest_star in constellations[:1]:
                visible = check_star_visibility_in_time_intervals(brightest_star, datetime_list, timezone)
                for line in visible:
                    file.write(line + "\n")

    except ValueError as e:
        print("Error:", e)


def compute_ra_dec_visibility(ra, dec):
    timezone = get_location_timezone(lat, lng)
    bucharest_location = EarthLocation(lat=lat * u.deg, lon=lng * u.deg, height=elev * u.m)

    astronomical_sunset, astronomical_sunrise = calculate_astronomical_night(lat, lng, elev, timezone)

    output_file = "visibility_results.txt"

    mins = IntervalType.THREE_HOUR
    datetime_list = generate_time_interval(astronomical_sunset, astronomical_sunrise, mins.value)

    with open(output_file, "w") as file:
        for dt in datetime_list:
            visibility_result = check_ra_dec_visibility(dt, bucharest_location, ra, dec)
            s_date_time = f"At datetime: {str(dt.astimezone(tz=timezone))[:-6]}"
            s_target = " - Andromeda Galaxy has the following attributes: \n"
            file.write(s_date_time)
            file.write(s_target)
            for key, value in visibility_result.items():
                file.write(f"{key}: {value}\n")
            file.write("\n")


def main():
    andromeda_ra, andromeda_dec = "0h42m44s", "+41d16m9s"
    compute_ra_dec_visibility(andromeda_ra, andromeda_dec)


main()
