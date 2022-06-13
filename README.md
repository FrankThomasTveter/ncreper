# ncreper

Welcome to the `ncreper` wiki!

## Running `ncreper.bin`


The `ncreper.bin` program can be used to preprocess NetCDF files and create new fields like wind speed from wind components, rain rate from accumulated rain, and return periods.

The `ncreper.bin` program is in the `nctools` module which can be loaded on the met.no ppi using

      module load nctools

The binary reads input from stdin. For instance:

      cd ncreper/work
      ncreper.bin < vindkast_meps_stasjoner.xml


## Program input
The input format is based on headers and data, and may look like this:

      # ?
      # exit
      #
       NCREPER V1.0 [0]
      #
       INPUT FILE (NETCDF) [1]
          input.nc
      #
       OUTPUT FILE (NETCDF) [1]
          output.nc
      #
       NUMBER OF DATA SLICES [1]VFL
          10
      #
       RETPER FILE (NETCDF) [*]
        meps_return_level_annual.nc
      #
       ZERO RETPER LEVEL [1]
          0.0D0
      #
      # IMPORT VARIABLES FROM RETPER FILE [*]
      #  land_area_fraction
      #
      #
       SPEED FROM VARX, VARY, MAX [1]VFLR
          wind_speed_10m      x_wind_10m       y_wind_10m        1000.0
          wind_speed_of_gust  x_wind_gust_10m  y_wind_gust_10m   1000.0
      #
       NEW VAR, OLD VAR, DIFFERENTIATE OVER HOURS, MIN, MAX [*]
          precipitation_60min     precipitation_amount_acc 1   0.0D0   1.0D15
          precipitation_120min    precipitation_amount_acc    120/60  0 99999
          precipitation_180min    precipitation_amount_acc    180/60  0 99999
          precipitation_360min    precipitation_amount_acc    360/60  0 99999
          precipitation_720min    precipitation_amount_acc    720/60  0 99999
          precipitation_1440min   precipitation_amount_acc   1440/60  0 99999
      #
       RETPER DIM[1]
          return_period
      #
       NEW VAR, OLD VAR, RETPER VAR [*]
          wind_speed_10m_return_period          wind_speed_10m          wind_speed_10m_return_level
          wind_speed_of_gust_return_period      wind_speed_of_gust      wind_speed_of_gust_return_level
          precipitation_60min_return_period     precipitation_60min     precipitation_amount_acc_1h_return_level
          precipitation_180min_return_period    precipitation_180min    precipitation_amount_acc_3h_return_level
          precipitation_360min_return_period    precipitation_360min    precipitation_amount_acc_6h_return_level
          precipitation_720min_return_period    precipitation_720min    precipitation_amount_acc_12h_return_level
          precipitation_1440min_return_period   precipitation_1440min   precipitation_amount_acc_24h_return_level
      #
      # PERCENTILE DIM, PERCENTILE [1]
      #   ensemble_member 95
      #
      # NEW VAR, OLD DIM VAR
      #    precipitation_60min_return_period_p95    precipitation_60min_return_period
      #    precipitation_60min_p95                  precipitation_60min

      #
       KEEP VARS [*]VFLR
          surface
          height0
          height7
          time
          forecast_reference_time
          x
          y
          latitude
          longitude
          projection_lambert
          projection_regular_ll
          land_area_fraction
      #
          precipitation_60min
          precipitation_60min_return_period    
          precipitation_180min
          precipitation_180min_return_period    
          precipitation_360min
          precipitation_360min_return_period    
          precipitation_720min
          precipitation_720min_return_period    
          precipitation_1440min
          precipitation_1440min_return_period    
      #
      #    precipitation_60min_p95
      #    precipitation_60min_return_period_p95
      #
          wind_speed_10m  wind_speed_10m_mb0
          wind_speed_10m_return_period wind_speed_10m_return_period_mb0
          wind_speed_of_gust  wind_speed_of_gust_mb0
          wind_speed_of_gust_return_period wind_speed_of_gust_return_period_mb0

Comments are preceded by `#`. The square brackets indicate how many lines the data body should consist off.

This program is very memory intensive. If you run out of memory, the program will stop, typically when reading a variable. In this case you have to request more memory from the `PPI`, or you have to slice the data more (`NUMBER OF DATA SLICES [1]VFL`).

# Local installation

To install `ncreper` locally, you first have to log onto the `PPI`,

      ssh ppi-clogin-b1

The you have to download the source code from github and compile

      git clone https://github.com/FrankThomasTveter/ncreper.git
      cd ncreper
      make
      ls ncreper/*.bin

# Updating the `nctools` module on `PPI`

The `ncreper.bin` software is part of the `nctools` module. You may install a new module by creating a new module file and directory in,

      cd /modules/centos7/user-apps/nctools
      cp -rf 0.34 XXX
      # cp   ~/ncreper.bin   XXX/bin/
      cd /modules/MET/centos7/user-modules/nctools
      cp 0.34 XXX
      emacs XXX

where XXX is the new version you have created.

