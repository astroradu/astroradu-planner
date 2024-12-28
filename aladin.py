from astroquery.vizier import Vizier

column_filter = ['Arp', 'Name', 'RAJ2000', 'DEJ2000', 'VT', 'dim1', 'dim2']


def check_if_file_empty(file_path):
    try:
        with open(file_path, 'r') as file_obj:
            first_char = file_obj.read(1)
            return not first_char
    except (FileNotFoundError, IOError):
        return True


def get_arp_galaxies_from_vizier():
    output_catalogues_file = "arp_catalogues.txt"
    output_arp_file = "catalogue_arp.csv"
    catalogue_id_arp = "VII/192"

    if check_if_file_empty(output_catalogues_file):
        catalog_list_arp = Vizier.find_catalogs("peculiar galaxies")
        with open(output_catalogues_file, "w") as file:
            for name, item in catalog_list_arp.items():
                file.write(f"Name: {name}: {item.description}\n")

    if check_if_file_empty(output_arp_file):
        catalogs_arp = Vizier(row_limit=-1).get_catalogs(catalogue_id_arp)
        table_arp = (catalogs_arp[1]).to_pandas()
        filtered_table = table_arp[column_filter]
        filtered_table.to_csv(output_arp_file, index=False)


def get_ngc_galaxies_from_vizier():
    output_ngc_file = "catalogue_ngc.csv"
    catalogue_id_ngc = "VII/118"

    if check_if_file_empty(output_ngc_file):
        catalogs_ngc = Vizier(row_limit=-1).get_catalogs(catalogue_id_ngc)
        table_ngc = (catalogs_ngc[0]).to_pandas()
        table_ngc.to_csv(output_ngc_file, index=False)


def aladin():
    # get_arp_galaxies_from_vizier()
    get_ngc_galaxies_from_vizier()


aladin()
