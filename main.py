
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy import units as u
from datetime import datetime, time
import ephem
import pytz
from timezonefinder import TimezoneFinder

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


def calculate_astronomical_night(utc_midnight, latitude, longitude, elevation):
    fred = ephem.Observer()
    fred.date = ephem.Date(utc_midnight)
    fred.lat = str(latitude)
    fred.lon = str(longitude)
    fred.elev = elevation
    fred.horizon = '0'
    fred.pressure = 0

    sunset = fred.previous_setting(ephem.Sun())  # Sunset
    sunrise = fred.next_rising(ephem.Sun())  # Sunrise

    print("===============================")
    print("Computed Upcoming Sunset: ")
    print(ephem.localtime(sunset))
    print("Computed Upcoming Sunrise: ")
    print(ephem.localtime(sunrise))

    # -6=civil twilight, -12=nautical, -18=astronomical
    fred.horizon = '-18'

    astronomical_sunset = fred.previous_setting(ephem.Sun(), use_center=True)  # Begin civil twilight
    astronomical_sunrise = fred.next_rising(ephem.Sun(), use_center=True)  # End civil twilight

    print("===============================")
    print("Computed Upcoming Astronomical Sunset: ")
    print(ephem.localtime(astronomical_sunset))
    print("Computed Upcoming Astronomical Sunrise: ")
    print(ephem.localtime(astronomical_sunrise))


def main():
    timezone = get_location_timezone(lat, lng)
    midnight = get_utc_midnight_for_tonight(timezone)
    calculate_astronomical_night(midnight, lat, lng, elev)


main()
