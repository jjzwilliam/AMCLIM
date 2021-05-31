import numpy as np
import sys
## daily simulation or hourly simulation
nhours = 24
Days = 365
Hours = nhours * Days
time = Days
timestep = 24

## resolution; e.g. 0.5 degree
dlat = 0.5
dlon = 0.5

## lat/lon dimensions
lats = int(180.0/dlat)
lons = int(360.0/dlon)

## levels: this refers to how many types of practices
## currently two: housing and MMS (manure management system)
levels = 2

## array dimensions;
#mtrx = [levels,time,lats,lons]
mtrx = [time,lats,lons]

## define months info
Months_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
        'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
Months_days = [31,28,31,30,31,30,31,31,30,31,30,31]
Months_idx = [0,31,59,90,120,151,181,212,243,273,304,334,365]

## specify livestock type and production sytem
livestock = 'PIG'
production_system = 'industrial'
