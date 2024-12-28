from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy.time import Time
from astropy import units as u
from datetime import datetime, time, timedelta
import ephem
import pytz
import json
import os

from timezonefinder import TimezoneFinder

from utils.intervalType import IntervalType

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
def calculate_astronomical_night(utc_midnight, latitude, longitude, elevation):
    observer = ephem.Observer()
    observer.date = ephem.Date(utc_midnight)
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


def check_star_visibility(brightest_star, datetime_list, timezone):
    star_coord = SkyCoord.from_name(brightest_star)
    results = []
    for dt in datetime_list:
        observer_alt_az = AltAz(obstime=dt,
                                location=EarthLocation(lat=44.4268 * u.deg, lon=26.1025 * u.deg,
                                                       height=85 * u.m))
        star_alt_az = star_coord.transform_to(observer_alt_az)
        star_visible = star_alt_az.alt > 0 * u.deg
        result = f"At datetime: {dt.astimezone(tz=timezone)} - " \
                 f"Star: {brightest_star} was {'visible' if star_visible else 'NOT visible'}."
        results.append(result)
    return results


def clear_file(f):
    with open(f, "w") as file:
        file.write("")
    return file


def main():
    timezone = get_location_timezone(lat, lng)
    midnight = get_utc_midnight_for_tonight(timezone)
    astronomical_sunset, astronomical_sunrise = calculate_astronomical_night(midnight, lat, lng, elev)
    constellations = read_constellations("constellations.json")
    output_file = "visibility_results.txt"

    try:

        print("----")

        mins = IntervalType.QUARTER

        datetime_list = generate_time_interval(astronomical_sunset, astronomical_sunrise, mins.value)

        with open(output_file, "w") as file:
            for constellation, brightest_star in constellations[:1]:
                visible = check_star_visibility(brightest_star, datetime_list, timezone)
                for line in visible:
                    file.write(line + "\n")

    except ValueError as e:
        print("Error:", e)


main()
