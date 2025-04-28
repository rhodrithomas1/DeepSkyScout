from astral import LocationInfo
from astral.sun import sun
from datetime import datetime
import pytz

def get_sun_times(lat, lon, date=None):
    city = LocationInfo(name="Custom", region="Custom", timezone="Europe/London", latitude=lat, longitude=lon)
    if date is None:
        date = datetime.now(pytz.timezone("Europe/London")).date()

    s = sun(city.observer, date=date, tzinfo=pytz.timezone("Europe/London"))

    return s['sunrise'], s['sunset']  # full datetime objects with timezone
