from astroquery.vizier import Vizier
import os


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


def get_arp_catalog_identifiers():
    output_catalogues_file = "data/vizier_catalogs.txt"
    if check_if_null_or_empty(output_catalogues_file):
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


def get_arp_galaxies_from_vizier():
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
        catalogs_arp = Vizier(row_limit=-1).get_catalogs(catalogue_id_arp)
        table_arp = (catalogs_arp[1]).to_pandas()
        filtered_table = table_arp[column_filter]
        table_size = len(filtered_table.index)
        subset = filtered_table.loc[0:table_size]
        subset = subset.copy()
        subset['size'] = subset[['dim1', 'dim2']].max(axis=1)
        subset = subset.drop(columns=['dim1', 'dim2'])
        subset.rename(columns={
            'Arp': 'name',
            'RAJ2000': 'ra',
            'DEJ2000': 'dec',
            'size': 'size',
            'VT': 'mag'
        }, inplace=True)
        col_to_move = 'size'
        col_to_move_data = subset.pop(col_to_move)
        subset.insert(3, col_to_move, col_to_move_data)
        sorted_subset = subset.sort_values(by='name', ascending=True)
        filtered_arp_table = sorted_subset.dropna()
        filtered_arp_table.to_csv(output_arp_file, index=False)


def get_ngc_and_ic_galaxies_from_vizier():
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

    catalogs_ngc = Vizier(row_limit=-1).get_catalogs(catalogue_id_ngc)
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


def get_sharpless_galaxies_from_vizier():
    output_sh_file = "data/catalog_sh.csv"
    catalogue_sh_ngc = "VII/20"
    column_filter = [
        'Sh2',
        'RA1900',
        'DE1900',
        'Diam',
        'Bright'
    ]

    # if not check_if_null_or_empty(output_sh_file):
    #     return None

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

    sorted_sh_table = sh_table.sort_values(by='name', ascending=True)
    filtered_sh_table = sorted_sh_table.dropna()
    filtered_sh_table.insert(0, 'cat', 'SH2')
    filtered_sh_table.to_csv(output_sh_file, index=False, header=True)

def get_abell_galaxies_from_vizier():
    output_abell_file = "data/catalog_abell.csv"
    catalogue_abell_ngc = "VII/110A"
    column_filter = [
        'Sh2',
        'RA1900',
        'DE1900',
        'Diam',
        'Bright'
    ]

    # if not check_if_null_or_empty(output_sh_file):
    #     return None

    catalog_abell = Vizier(row_limit=40).get_catalogs(catalogue_abell_ngc)
    table_abell = (catalog_abell[0]).to_pandas()
    # filtered_table = table_sh[column_filter]
    # table_size = len(filtered_table.index)
    # sh_table = filtered_table.loc[0:table_size]
    # sh_table = sh_table.copy()
    # sh_table.rename(columns={
    #     'Sh2': 'name',
    #     'RA1900': 'ra',
    #     'DE1900': 'dec',
    #     'Diam': 'size',
    #     'Bright': 'bri'
    # }, inplace=True)
    #
    # sorted_sh_table = sh_table.sort_values(by='name', ascending=True)
    # filtered_sh_table = sorted_sh_table.dropna()
    # filtered_sh_table.insert(0, 'cat', 'SH2')
    table_abell.to_csv(output_abell_file, index=False, header=True)

def aladin():
    # get_arp_catalog_identifiers()
    # get_arp_galaxies_from_vizier()
    # get_ngc_and_ic_galaxies_from_vizier()
    get_abell_galaxies_from_vizier()


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