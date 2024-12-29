from collections import OrderedDict
from operator import itemgetter

import pandas as pd
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
from utils.time_utils import current_milli_time, current_milli_time_formatted, format_mills

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


def compute_ra_dec_attributes(dt, loc, ra, dec):
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


def compute_target_potential_score(attributes):
    """
    To help compare multiple targets for astrophotography potential for a night or even period , we must calculate a
    score for each target and for each time segment we observe (see IntervalType.py for the options). If for example
    we are measuring the potential for a 5-hour astronomical night, and we fragment into half-hour intervals,
    the algorithm will compute 10 different scores for the target in that night and will add them up to create a
    master score for the night. When comparing to other targets' master score, we can easily decide which target is
    worth the effort and which is not.

    If a target is not visible or is below 25 degrees Altitude (too low for quality imagery), we break the
    computation with a 0 score.

    If the moon is visible but the conditions are not acceptable (moon phase too advanced or moon closer than 90
    degrees Azimuth to our target), then again we return a score of 0. Astr

    The closer the target is to Zenith at the observed time, the higher score we add up. If the target is within 10
    degreez Altitude of Zenith, we add 3 score points for exceptional imagery potential.

    Astrophotography done under an acceptable moon should also be rewarded, so if the moon is visible and the moon
    conditions are good, the target adds a higher score if its position is antipodal (other side of the sky dome) to
    the moon position. Narrowband signal acquisition is recommended during full moon to lower the effects of moon light.
    """
    if (
            not attributes["target_visible"] or
            not attributes["is_target_above_25_degree"] or
            (
                    attributes["is_moon_visible"] and not attributes["is_moon_phase_acceptable"]
            )
    ):
        return 0

    # Initialize score
    score = 1  # Default score if the target is visible and above 25 degrees

    # Scoring based on target's position relative to zenith
    if attributes["is_target_in_45_degrees_of_zenith"]:
        if not attributes["is_target_in_30_degrees_of_zenith"]:
            score += 1
        elif not attributes["is_target_in_10_degrees_of_zenith"]:
            score += 2
        else:
            score += 3

    # Adjust score based on moon visibility and position
    if attributes["is_moon_visible"] and attributes["is_moon_phase_acceptable"]:
        if not attributes["is_moon_in_180_degrees"]:
            score += 2
    else:
        score += 1

    return score


def test_attribute_computation():
    attributes = {
        "target_alt": None,
        "target_az": None,
        "target_visible": True,
        "is_target_above_25_degree": True,
        "is_target_in_10_degrees_of_zenith": False,
        "is_target_in_30_degrees_of_zenith": True,
        "is_target_in_45_degrees_of_zenith": True,
        "moon_alt": None,
        "moon_az": None,
        "moon_phase": 30,
        "is_moon_visible": True,
        "is_moon_phase_acceptable": True,
        "is_moon_in_90_degrees": False,
        "is_moon_in_180_degrees": True
    }

    score = compute_target_potential_score(attributes)
    return score


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


def compute_ra_dec_score(ra, dec):
    timezone = get_location_timezone(lat, lng)
    bucharest_location = EarthLocation(lat=lat * u.deg, lon=lng * u.deg, height=elev * u.m)

    astronomical_sunset, astronomical_sunrise = calculate_astronomical_night(lat, lng, elev, timezone)

    output_file = "visibility_results.txt"

    mins = IntervalType.THREE_HOUR
    datetime_list = generate_time_interval(astronomical_sunset, astronomical_sunrise, mins.value)

    master_score = 0

    with open(output_file, "w") as file:
        for dt in datetime_list:
            attributes = compute_ra_dec_attributes(dt, bucharest_location, ra, dec)
            score = compute_target_potential_score(attributes)
            master_score += score
            s_date_time = f"At datetime: {str(dt.astimezone(tz=timezone))[:-6]}"
            s_target = " - Andromeda Galaxy has the following attributes: \n"

            file.write(s_date_time)
            file.write(s_target)
            for key, value in attributes.items():
                file.write(f"{key}: {value}\n")
            file.write(f"================ SCORE:{score}\n")
            file.write("\n")


def compute_constellations_scores():
    timezone = get_location_timezone(lat, lng)
    bucharest_location = EarthLocation(lat=lat * u.deg, lon=lng * u.deg, height=elev * u.m)

    astronomical_sunset, astronomical_sunrise = calculate_astronomical_night(lat, lng, elev, timezone)
    constellations = read_constellations("constellations.json")
    output_file = "visibility_results.txt"

    mins = IntervalType.THREE_HOUR
    datetime_list = generate_time_interval(astronomical_sunset, astronomical_sunrise, mins.value)

    constellations_scores = {}

    with open(output_file, "w") as file:
        file.write(f"Master scores for tonight:\n\n")
        start_timestamp = current_milli_time()

        for constellation, brightest_star in constellations:
            star_coord = SkyCoord.from_name(brightest_star)
            master_score = 0
            for dt in datetime_list:
                attributes = compute_ra_dec_attributes(dt, bucharest_location, star_coord.ra, star_coord.dec)
                score = compute_target_potential_score(attributes)
                master_score += score
                # s_date_time = f"At datetime: {str(dt.astimezone(tz=timezone))[:-6]}"
                # s_target = " - Andromeda Galaxy has the following attributes: \n"
                #
                # file.write(s_date_time)
                # file.write(s_target)
                # for key, value in attributes.items():
                #     file.write(f"{key}: {value}\n")
                # file.write(f"================ SCORE:{score}\n")
                # file.write("\n")
            constellations_scores[brightest_star] = master_score

        sorted_constellations_scores = dict(sorted(constellations_scores.items(), key=lambda item: item[1]))

        end_timestamp = current_milli_time()

        for star, master_sc in constellations_scores.items():
            file.write(f"Star: {star}: Score {master_sc}\n")

        file.write(f"\nFinished computation in {(end_timestamp - start_timestamp) / 1000} seconds.")


def compute_constellations_validity():
    timezone = get_location_timezone(lat, lng)
    astronomical_sunset, astronomical_sunrise = calculate_astronomical_night(lat, lng, elev, timezone)
    constellations = read_constellations("constellations.json")
    output_file = "visibility_results.txt"

    mins = IntervalType.QUARTER

    datetime_list = generate_time_interval(astronomical_sunset, astronomical_sunrise, mins.value)
    bucharest_location = EarthLocation(lat=lat * u.deg, lon=lng * u.deg, height=elev * u.m)

    resultDictionary = {}

    with open(output_file, "w") as file:
        for constellation, brightest_star in constellations[:1]:
            star_coord = SkyCoord.from_name(brightest_star)

            # for dt in datetime_list:
            #     visibility_result = check_ra_dec_visibility(dt, bucharest_location, star_coord.ra, star_coord.dec)
            #     resultDictionary[constellation]
            #     for line in visible:
            #         file.write(line + "\n")


def compute_scores():
    brightness = 3

    timezone = get_location_timezone(lat, lng)
    astronomical_sunset, astronomical_sunrise = calculate_astronomical_night(lat, lng, elev, timezone)

    timestamp_file = "data/timestamps.txt"
    score_file = "data/scores.txt"

    mins = IntervalType.HOUR

    datetime_list = generate_time_interval(astronomical_sunset, astronomical_sunrise, mins.value)
    bucharest_location = EarthLocation(lat=lat * u.deg, lon=lng * u.deg, height=elev * u.m)

    begin_timestamp = current_milli_time()

    resultDictionary = {}
    timestamps = [f"Reading tables: {format_mills(begin_timestamp)}\n"]

    ngc_df = pd.read_csv("data/catalog_ngc.csv")
    ic_df = pd.read_csv("data/catalog_ic.csv")
    combined_df = pd.concat([ngc_df, ic_df], ignore_index=True)

    timestamps.append(f"Starting iteration...\n")

    night_segment_size = len(datetime_list)
    table_size = len(combined_df.index)

    current_row = 0
    current_dt = 0

    for index, row in combined_df.iterrows():
        current_row += 1
        current_dt = 0

        if row['bri'] != brightness:
            pass

        identifier = f"{row['cat']}{row['name']}"

        target_score = 0

        for dt in datetime_list:
            current_dt += 1
            timestamps.append(
                f"Computing attributes for {identifier} at {dt.time()}: {current_milli_time_formatted()}\n")
            attributes = compute_ra_dec_attributes(dt, bucharest_location, row['ra'], row['dec'])
            timestamps.append(
                f"Computing scores for {identifier} at {dt.time()}: {current_milli_time_formatted()}\n")
            score = compute_target_potential_score(attributes)
            timestamps.append(f"Resulting attributes for {dt.time()}:\n")
            for key, value in attributes.items():
                timestamps.append(
                    f"{key} : {value}\n")
            timestamps.append("\n")
            target_score += score
            print(f"{current_row}/{table_size} | {current_dt}/{night_segment_size}")

        timestamps.append(f"Writing master score for {identifier}: {current_milli_time_formatted()}\n\n")
        resultDictionary[f"{identifier}"] = target_score

    timestamps.append(f"Finishing...\n\n")
    finish_timestamp = current_milli_time()
    diff = (finish_timestamp - begin_timestamp) // 1000
    timestamps.append(f"Computation took ... {diff}s\n")
    sorted_dict = OrderedDict(sorted(resultDictionary.items(), key=itemgetter(1), reverse=True))

    with open(score_file, "w") as file:
        for key, value in sorted_dict.items():
            file.write(f"{key}: {value}" + "\n")

    with open(timestamp_file, "w") as file:
        for line in timestamps:
            file.write(line)


def main():
    compute_scores()


main()
