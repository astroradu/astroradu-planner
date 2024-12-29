from astropy.coordinates import SkyCoord, FK4, FK5
from astroquery.vizier import Vizier
import os
import pandas as pd
import astropy.units as u


def check_if_null_or_empty(file_path):
    try:
        root_dir = os.path.dirname(os.path.abspath(__file__))
        new_folder = os.path.join(root_dir, "data")
        os.makedirs(new_folder, exist_ok=True)
        with open(file_path, 'r') as file_obj:
            first_char = file_obj.read(1)
            return not first_char
    except (FileNotFoundError, IOError):
        return True


def convert_fk4_frame_to_icrs(ra_b1900, dec_b1900):
    frame_b1900 = SkyCoord(ra=ra_b1900, dec=dec_b1900, frame=FK4)
    frame_icrs = frame_b1900.transform_to(FK4(equinox='J2000'))
    return frame_icrs.ra, frame_icrs.dec


def convert_fk5_frame_to_icrs(ra_b2000, dec_b2000):
    frame_b2000 = SkyCoord(ra=ra_b2000, dec=dec_b2000, frame='fk5', equinox='B2000', unit=(u.hourangle, u.deg))
    frame_icrs = frame_b2000.transform_to('icrs')
    return frame_icrs.ra, frame_icrs.dec


def convert_fk5_frame_to_galactic(ra_b2000, dec_b2000):
    frame_b2000 = SkyCoord(ra=ra_b2000 * u.deg, dec=dec_b2000 * u.deg, frame=FK5(equinox='B2000'),
                           unit=(u.hourangle, u.deg), )
    frame_galactic = frame_b2000.transform_to('galactic')
    return frame_galactic.l.deg, frame_galactic.b.deg


def convert_ngc_ra_dec_to_icrs(row):
    ra_clean = row['ra'].replace(' ', ':')
    dec_clean = row['dec'].replace(' ', ':')
    coord = SkyCoord(ra=ra_clean, dec=dec_clean, frame=FK5(equinox='B2000'), unit=(u.hourangle, u.deg))
    frame_icrs = coord.transform_to('icrs')
    ra_dec_string = frame_icrs.to_string('hmsdms', precision=0).split()
    return ra_dec_string[0], ra_dec_string[1]


def convert_sh_ra_dec_to_icrs(row):
    ra_clean = row['ra'].replace(' ', ':')
    dec_clean = row['dec'].replace(' ', ':')
    coord = SkyCoord(ra=ra_clean, dec=dec_clean, frame=FK4(equinox='B1900'), unit=(u.hourangle, u.deg))
    frame_icrs = coord.transform_to('icrs')
    ra_dec_string = frame_icrs.to_string('hmsdms', precision=0).split()
    return ra_dec_string[0], ra_dec_string[1]


def convert_arp_ra_dec_to_icrs(row):
    ra_clean = row['ra'].replace(' ', ':')
    dec_clean = row['dec'].replace(' ', ':')
    coord = SkyCoord(ra=ra_clean, dec=dec_clean, unit=(u.hourangle, u.deg))
    ra_dec_string = coord.to_string('hmsdms', precision=0).split()
    return ra_dec_string[0], ra_dec_string[1]


def get_catalogs():
    output_catalogues_file = "data/vizier_catalogs.txt"
    if not check_if_null_or_empty(output_catalogues_file):
        return

    keywords = [
        "Messier",
        "NGC",
        "UGC",
        "Sharpless",
        "Abell",
        "peculiar galaxies"
    ]
    catalog_list = {}
    for keyword in keywords:
        found_catalogs = Vizier.find_catalogs(keyword, max_catalogs=1000)
        catalog_list.update(found_catalogs)

    with open(output_catalogues_file, "w") as file:
        for name, item in catalog_list.items():
            file.write(f"Name: {name}: {item.description}\n")


def get_arp_objects_from_vizier():
    output_arp_file = "data/catalog_arp.csv"
    catalogue_id_arp = "VII/192"
    column_filter = [
        'Arp',
        'RAJ2000',
        'DEJ2000',
        'VT',
        'dim1',
        'dim2'
    ]

    if check_if_null_or_empty(output_arp_file):
        catalogs_arp = Vizier(row_limit=40).get_catalogs(catalogue_id_arp)
        table_arp = (catalogs_arp[1]).to_pandas()
        filtered_table = table_arp[column_filter]
        table_size = len(filtered_table.index)
        arp_table = filtered_table.loc[0:table_size]
        arp_table = arp_table.copy()
        arp_table['size'] = arp_table[['dim1', 'dim2']].max(axis=1)
        arp_table = arp_table.drop(columns=['dim1', 'dim2'])
        arp_table.rename(columns={
            'Arp': 'name',
            'RAJ2000': 'ra',
            'DEJ2000': 'dec',
            'size': 'size',
            'VT': 'mag'
        }, inplace=True)
        arp_table[['ra', 'dec']] = arp_table.apply(
            lambda row: (
                convert_arp_ra_dec_to_icrs(row)
            ),
            axis=1,
            result_type='expand'
        )
        col_to_move = 'size'
        col_to_move_data = arp_table.pop(col_to_move)
        arp_table.insert(3, col_to_move, col_to_move_data)
        sorted_subset = arp_table.sort_values(by='name', ascending=True)
        filtered_arp_table = sorted_subset.dropna()
        filtered_arp_table.to_csv(output_arp_file, index=False)


# <9 e bright, >13 e faint
def get_ngc_and_ic_objects_from_vizier():
    output_ngc_file = "data/catalog_ngc.csv"
    output_ic_file = "data/catalog_ic.csv"
    catalogue_id_ngc = "VII/118"
    column_filter = [
        'Name',
        'RAB2000',
        'DEB2000',
        'size',
        'mag'
    ]

    if not check_if_null_or_empty(output_ngc_file) and not check_if_null_or_empty(output_ic_file):
        return None

    catalogs_ngc = Vizier(row_limit=300).get_catalogs(catalogue_id_ngc)
    table_ngc = (catalogs_ngc[0]).to_pandas()
    filtered_table = table_ngc[column_filter]
    table_size = len(filtered_table.index)
    ngc_table = filtered_table.loc[0:table_size]
    ngc_table = ngc_table.copy()
    ngc_table.rename(columns={
        'Name': 'name',
        'RAB2000': 'ra',
        'DEB2000': 'dec',
        'size': 'size',
        'mag': 'mag'
    }, inplace=True)

    ngc_table[['ra', 'dec']] = ngc_table.apply(
        lambda row: (
            convert_ngc_ra_dec_to_icrs(row)
        ),
        axis=1,
        result_type='expand'
    )

    ngc_table['bri'] = ngc_table['mag'].apply(
        lambda mag: 3 if mag < 9 else (1 if mag > 12 else 2)
    )

    ngc_table.drop(columns=['mag'], inplace=True)

    ic_table = ngc_table[ngc_table['name'].str.lower().str.startswith('i')]
    ic_table.loc[:, 'name'] = ic_table['name'].str.replace(r'[iI]', '', regex=True)

    ngc_table = ngc_table[~ngc_table['name'].str.lower().str.startswith('i')]

    sorted_ngc_table = ngc_table.sort_values(by='name', ascending=True)
    sorted_ic_table = ic_table.sort_values(by='name', ascending=True)

    filtered_ngc_table = sorted_ngc_table.dropna()
    filtered_ic_table = sorted_ic_table.dropna()

    filtered_ngc_table.insert(0, 'cat', 'NGC')
    filtered_ic_table.insert(0, 'cat', 'IC')

    filtered_ngc_table.to_csv(output_ngc_file, index=False, header=True)
    filtered_ic_table.to_csv(output_ic_file, index=False, header=True)


def get_sharpless_objects_from_vizier():
    output_sh_file = "data/catalog_sh.csv"
    catalogue_sh_ngc = "VII/20"
    column_filter = [
        'Sh2',
        'RA1900',
        'DE1900',
        'Diam',
        'Bright'
    ]

    if not check_if_null_or_empty(output_sh_file):
        return None

    catalog_sh = Vizier(row_limit=40).get_catalogs(catalogue_sh_ngc)
    table_sh = (catalog_sh[0]).to_pandas()
    filtered_table = table_sh[column_filter]
    table_size = len(filtered_table.index)
    sh_table = filtered_table.loc[0:table_size]
    sh_table = sh_table.copy()
    sh_table.rename(columns={
        'Sh2': 'name',
        'RA1900': 'ra',
        'DE1900': 'dec',
        'Diam': 'size',
        'Bright': 'bri'
    }, inplace=True)

    sh_table[['ra', 'dec']] = sh_table.apply(
        lambda row: (
            convert_sh_ra_dec_to_icrs(row)
        ),
        axis=1,
        result_type='expand'
    )

    sorted_sh_table = sh_table.sort_values(by='name', ascending=True)
    filtered_sh_table = sorted_sh_table.dropna()
    filtered_sh_table.insert(0, 'cat', 'SH2')
    filtered_sh_table.to_csv(output_sh_file, index=False, header=True)


def aladin():
    get_catalogs()
    get_ngc_and_ic_objects_from_vizier()
    get_sharpless_objects_from_vizier()
    get_arp_objects_from_vizier()


aladin()

# Convert B1950 to J2000

# from astLib import astCoords as coords
#
# #Convert B1950 coordinates (as given)
# clra=zeros(3)
# cldec=zeros(3)
# for i in range(len(clra)):
#     clra[i], cldec[i] = coords.convertCoords(
#                                 'B1950', 'J2000', clra1950[i],
#                                 cldec1950[i], 1950
#                             )
#
# #Convert input coords to Galactic coords
# lgal,bgal = radians(coords.convertCoords(
#                         'J2000', 'GALACTIC', ra,
#                         dec,2000
#                     ))
#
# In [1]: import numpy as np
#
# In [2]: from astropy import units as u
#
# In [3]: from astropy.coordinates import SkyCoord, FK4, FK5, Galactic
#
# In [4]: clra = np.zeros(3)
#
# In [5]: cldec = np.zeros(3)
#
# In [6]: c1 = SkyCoord(clra * u.deg, cldec * u.deg, frame=FK4)
#
# In [7]: c2 = c1.transform_to(FK5(equinox='J2000'))
#
# In [8]: c3 = c2.transform_to(Galactic)
#
# In [9]: print(c3.l.degree)
# [ 97.74220094  97.74220094  97.74220094]
#
# In [10]: print(c3.b.degree)
# [-60.18102359 -60.18102359 -60.18102359]
