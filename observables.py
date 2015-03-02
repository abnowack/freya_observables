"""
TODO:
[ ] Extract terminated neutrons, record position, velocity vectors and time
[ ] Do same for photons
[ ] Create method for storing neutron and photon lists for each event
"""

import numpy as np
import matplotlib.pyplot as plt
from mcnpy.mcnp_wrapper import *
from mcnpy.ptrac import reader as preader
import h5py as h5

def parse_ptrac_to_hdf5(ptrac_filename):
    ''' save neutron and photon position, directory, energy, time '''
    with open(ptrac_filename, 'r') as ptrac:
        # parse headers and formats
        header = preader.ptrac_header(ptrac)
        input_format = preader.ptrac_input_format(ptrac)
        event_format = preader.ptrac_event_format(ptrac)
        
        with h5.File('ptrac.h5', 'w') as f:
            f_events_g = f.create_group('events')
            
            dt = np.dtype([('x', np.float), ('y', np.float), ('z', np.float),
                           ('u', np.float), ('v', np.float), ('w', np.float),
                           ('energy', np.float), ('time', np.float)])
        
            for i, history in enumerate(preader.parse_ptrac_events(ptrac, event_format)):
                f_history_g = f_events_g.create_group(str(i))
                neutrons = []
                photons = []
                for ev in history.events:
                    if ev.ipt == 1.0:
                        neutrons.append(ev)
                    elif ev.ipt == 2.0:
                        photons.append(ev)

                f_neutron_d = f_history_g.create_dataset('neutrons', (len(neutrons),), dtype=dt)
                f_photon_d = f_history_g.create_dataset('photons', (len(photons),), dtype=dt)
                
                if len(photons) > 0:
                    print i
                
                for j, neutron in enumerate(neutrons):
                    f_neutron_d[j, 'x'] = neutron.xxx
                    f_neutron_d[j, 'y'] = neutron.yyy
                    f_neutron_d[j, 'z'] = neutron.zzz
                    f_neutron_d[j, 'u'] = neutron.uuu
                    f_neutron_d[j, 'v'] = neutron.vvv
                    f_neutron_d[j, 'w'] = neutron.www
                    f_neutron_d[j, 'energy'] = neutron.erg
                    f_neutron_d[j, 'time'] = neutron.tme

                for j, photon in enumerate(photons):
                    f_photon_d[j, 'x'] = photon.xxx
                    f_photon_d[j, 'y'] = photon.yyy
                    f_photon_d[j, 'z'] = photon.zzz
                    f_photon_d[j, 'u'] = photon.uuu
                    f_photon_d[j, 'v'] = photon.vvv
                    f_photon_d[j, 'w'] = photon.www
                    f_photon_d[j, 'energy'] = photon.erg
                    f_photon_d[j, 'time'] = photon.tme
  
with open('input_cards/geometry.i', 'r') as cardfile:
    card = cardfile.read()

with run_mcnp(card, cores=4) as (status, mcnp_dir):
    parse_ptrac_to_hdf5(mcnp_dir + '\\ptrac')