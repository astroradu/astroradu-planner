
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, get_sun
from astropy.time import Time
# from astropy.units import degree
import astropy.units as u
import numpy as np
from datetime import datetime, timedelta

# Observer's location: Bucharest, Romania
location = EarthLocation(lat=44.4268 * u.deg, lon=26.1025 * u.deg, height=85 * u.m)

# Observation span: Next three months
start_date = datetime(2024, 12, 25)
end_date = datetime(2025, 3, 25)
num_days = (end_date - start_date).days

# List of all northern hemisphere constellations and their brightest stars
constellation_stars = [
    ("Andromeda", "Alpheratz"), ("Auriga", "Capella"), ("Boötes", "Arcturus"),
    ("Camelopardalis", "Caph"), ("Canes Venatici", "Cor Caroli"), ("Cassiopeia", "Schedar"),
    ("Cepheus", "Alderamin"), ("Cygnus", "Deneb"), ("Draco", "Thuban"),
    ("Hercules", "Kornephoros"), ("Lacerta", "Alpha Lacertae"), ("Leo", "Regulus"),
    ("Leo Minor", "Praecipua"), ("Lyra", "Vega"), ("Monoceros", "Beta Monocerotis"),
    ("Orion", "Betelgeuse"), ("Pegasus", "Enif"), ("Perseus", "Mirfak"),
    ("Sagitta", "Delta Sagittae"), ("Taurus", "Aldebaran"), ("Triangulum", "Mothallah"),
    ("Ursa Major", "Dubhe"), ("Ursa Minor", "Polaris"), ("Vulpecula", "Anser"),
]

# Initialize dictionary to store total visibility duration
visibility_duration = {constellation: 0 for constellation, _ in constellation_stars}

# Calculate astronomical night and visibility for each day
for i in range(1):
    print(f"Processing day {i + 1} of {num_days}...")

    # Current date
    current_date = start_date + timedelta(days=i)
    time = Time(current_date)

    # Calculate astronomical twilight start and end
    midnight = Time(current_date + timedelta(hours=12))  # Roughly midday as a starting point
    sun_altaz = get_sun(midnight).transform_to(AltAz(obstime=midnight, location=location))
    astronomical_twilight_angle = -18 * u.deg

    # Calculate sunset and sunrise for the astronomical twilight
    sunset = time + (midnight - time) * (sun_altaz.alt - astronomical_twilight_angle) / sun_altaz.alt
    sunrise = time + (midnight - time) * (astronomical_twilight_angle - sun_altaz.alt) / -sun_altaz.alt

    # Observation times during the astronomical night
    obs_times = sunset + np.arange(0, (sunrise - sunset).sec, 15 * 60) * u.s

    print(obs_times)
    # Calculate visibility for each constellation
    for constellation, star_name in constellation_stars:

        if (constellation == "Andromeda"):
            # Get the star's coordinates
            star_coord = SkyCoord.from_name(star_name)

            # Create AltAz frame for the observer's location and time
            altaz_frame = AltAz(obstime=obs_times, location=location)

            # Transform star's coordinates to AltAz
            star_altaz = star_coord.transform_to(altaz_frame)

            # Check visibility: Altitude > 20 degrees
            visible = star_altaz.alt > 20 * u.deg
            visibility_duration[constellation] += np.sum(visible) * (15 * u.min).to(u.hour).value

# Sort constellations by total visibility duration in ascending order
sorted_visibility = sorted(visibility_duration.items(), key=lambda x: x[1])

# Print all constellations in reverse order
print("\nConstellations sorted by total visibility time above 20 degrees altitude (ascending):")
for i, (constellation, duration) in enumerate(sorted_visibility):
    print(f"{36 - i}. {constellation}: {duration:.2f} hours")