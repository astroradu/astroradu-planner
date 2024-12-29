import time
from datetime import datetime


def current_milli_time():
    return round(time.time() * 1000)


def current_milli_time_formatted():
    time_mills = current_milli_time()
    dt_object = datetime.fromtimestamp(time_mills / 1000.0)
    return dt_object.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]


def format_mills(mills):
    dt_object = datetime.fromtimestamp(mills / 1000.0)
    return dt_object.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
