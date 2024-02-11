import icepyx as ipx

import numpy as np
import xarray as xr
import pandas as pd

import h5py
import os,json
from pprint import pprint

region_a = ipx.Query('ATL08',[25,-31,35,-20],['2019-02-22','2023-02-28'], start_time='00:00:00', end_time='23:59:59')
region_a.earthdata_login('icepyx_devteam','icepyx.dev@gmail.com')
region_a.order_vars.avail()
region_a.order_vars.avail(options=True)
