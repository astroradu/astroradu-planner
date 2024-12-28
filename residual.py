import json
import os


def read_constellations(file_name):
    file_path = os.path.join(os.path.dirname(__file__), file_name)
    with open(file_path, 'r') as file:
        json_data = json.load(file)
    return tuple(json_data.items())
