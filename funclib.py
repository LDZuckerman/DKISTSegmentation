import scipy.io as sio
from scipy.io import readsav
import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.coordinates import frames
import matplotlib.pyplot as plt
import numpy as np

# Function for out own use in converted .sav data to SunPy maps
def sav_to_maps(filename):

    fake_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime='2013-10-28 08:24',
                 observer='earth', frame=frames.Helioprojective)
    fake_header = sunpy.map.make_fitswcs_header(data=np.empty((512, 512)), coordinate=fake_coord,
                                       reference_pixel=[0, 0]*u.pixel,
                                       scale=[2, 2]*u.arcsec/u.pixel,
                                       telescope='Fake Telescope', instrument='Fake Instrument',
                                       wavelength=1000*u.angstrom)

    data = readsav(filename)
    flux_4170 = sunpy.map.Map(data['rosa_wl'], fake_header)
    flux_4310 = sunpy.map.Map(data['rosa_gband'], fake_header)
    flux_7200 = sunpy.map.Map(data['ibis_wl'], fake_header)
    velocity =  sunpy.map.Map(data['velocity'], fake_header) # velocities in km/sec in 7090 Å spectral line
    velocity_clean =  sunpy.map.Map(data['velocity_cln'], fake_header)#  velocities after removal of the large scale structures (which are related more to oscillations than the granulation/convective overshoot)
    
    return flux_4170, flux_4310, flux_4170, flux_7200, velocity, velocity_clean

def sav_to_numpy(filename, instrument, band):
    """ Read .sav file data into a numpy array.
        Parameters:
        ----------

        filename (string): Path to input data file (.sav format)
        instrument (string): The telescope/instrument (ie. IBIS, DKSIST).
        band (string): Band or data field to read.

        Returns:
        -------

        data (numpy array): Array of data values.
    """


    if instrument != 'IBIS':
        print('This functionality has so far only been implemented ' + 
              'for IBIS data.')
        sys.ext(1)
    
    bands = ['rosa_wl', 'rosa_gband', 'ibis_wl', 'velocity', 'velocity_cln']
    if band not in bands:
        print('Band must be one of: ')
        print('rosa_wl      --> 4170 A')
        print('rosa_gband   --> 4310 A')
        print('ibis_wl      --> 7200 A')
        print('velocity     --> velocity in km/sec at 7090 Å')
        print('velocity_cln --> velocity with large-scale structure removed')
        sys.exit(1)

    # TO DO: catch error if file is not in sav format
    file = sio.readsav(filename)
    data = file[band]
    
    return data

    
    
        
  

        

    
