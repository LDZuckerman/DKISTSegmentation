import scipy.io as sio
from scipy.io import readsav
import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.coordinates import frames
import matplotlib.pyplot as plt

# Function for out own use in converted .sav data to SunPy maps
def sav_to_maps(filename):

    fake_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime='2013-10-28 08:24',
                 observer='earth', frame=frames.Helioprojective)
    fake_header = sunpy.map.make_fitswcs_header(data, fake_coord,
                                       reference_pixel=[0, 0]*u.pixel,
                                       scale=[2, 2]*u.arcsec/u.pixel,
                                       telescope='Fake Telescope', instrument='UV detector',
                                       wavelength=1000*u.angstrom)

    data = readsav(filename)
    flux_4170 = sunpy.map.Map(data['rosa_wl'], fake_header)
    flux_4310 = sunpy.map.Map(data['rosa_wgband'], fake_header)
    flux_7200 = sunpy.map.Map(data['ibis_wl'], fake_header)
    velocity =  sunpy.map.Map(data['velocity'], fake_header) # velocities in km/sec in 7090 Å spectral line
    velocity_clean =  sunpy.map.Map(data['velocity_clean'], fake_header)#  velocities after removal of the large scale structures (which are related more to oscillations than the granulation/convective overshoot)
    
    return flux_4170, flux_4310, flux_4170, flux_7200, velocity, velocity_clean
  

        

    
