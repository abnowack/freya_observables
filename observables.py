"""
TODO:
Might consider not saving in event folder, but with event field sequentially
[ ] Finalize HDF5 Format
[ ] Automated ways to plot
    [ ] Avg Neutrons
    [ ] Avg Photons
    [ ] Time Histogram of Neutrons
    [ ] Time Histogram of Photons
    [ ] Angle between neutrons
    [ ] Angle between photons
    [ ] Feynamn Distribution
"""

import numpy as np
import matplotlib.pyplot as plt
from mcnpy.mcnp_wrapper import *
from mcnpy.ptrac import reader as preader
import h5py as h5

def parse_ptrac_to_hdf5(ptrac_filename, buffer_len=1000):
    ''' save neutron and photon position, directory, energy, time '''
    with open(ptrac_filename, 'r') as ptrac:
        # parse headers and formats
        header = preader.ptrac_header(ptrac)
        input_format = preader.ptrac_input_format(ptrac)
        event_format = preader.ptrac_event_format(ptrac)
        
        with h5.File('ptrac_.h5', 'w') as f:
            f_events_g = f.create_group('events')
            
            dt = np.dtype([('evt', np.int), ('ipt', np.int), ('x', np.float), 
                           ('y', np.float), ('z', np.float), ('u', np.float), 
                           ('v', np.float), ('w', np.float),
                           ('energy', np.float), ('time', np.float)])
            
            f_events_d = f_events_g.create_dataset('events', (buffer_len,), dtype=dt, maxshape=(None,), chunks=True)
            nevt = 0
            buffer_count = 1
            for i, history in enumerate(preader.parse_ptrac_events(ptrac, event_format)):
                if (buffer_len*buffer_count - nevt - len(history.events)) < 0:
                    buffer_count += 1
                    f_events_d.resize((buffer_len * buffer_count,))                 
                
                for ev in history.events:
                    f_events_d[nevt, 'evt'] = i
                    f_events_d[nevt, 'ipt'] = int(ev.ipt)
                    f_events_d[nevt, 'x'] = ev.xxx
                    f_events_d[nevt, 'y'] = ev.yyy
                    f_events_d[nevt, 'z'] = ev.zzz
                    f_events_d[nevt, 'u'] = ev.uuu
                    f_events_d[nevt, 'v'] = ev.vvv
                    f_events_d[nevt, 'w'] = ev.www
                    f_events_d[nevt, 'energy'] = ev.erg
                    f_events_d[nevt, 'time'] = ev.tme
                    nevt += 1
            
            f_events_d.resize((nevt,))
                    
if __name__ == '__main__':
    
#    with open('input_cards/geometry.i', 'r') as cardfile:
#        card = cardfile.read()
#
#    with run_mcnp(card, cores=4) as (status, mcnp_dir):
#        parse_ptrac_to_hdf5(mcnp_dir + '\\ptrac')
    
    # calculate nu
    print 'processing'
    with h5.File('ptrac_.h5', 'r') as f:
        evts = f['events/events']
        history_index = evts[:-1, 'evt'] != evts[1:, 'evt']
        split_evts = np.array(np.split(evts, np.where(history_index)[0]))
        neutrons = split_evts * split_evts[:, 'ipt' == 1]
        photons = split_evts * split_evts[:, 'ipt' == 2]
        n_nu, n_bins = np.histogram((np.size(n) for n in neutrons), bins=50, range=[0,50])
        p_nu, p_bins = np.histogram((np.size(n) for n in neutrons), bins=50, range=[0,50])

    plt.bar(n_bins, n_nu, color='r', alpha=0.5)
    plt.bar(p_bins, p_nu, color='b', alpha=0.5)
    plt.yscale('log')
    plt.show()
    
    print np.average(np.arange(len(n_nu)), weights=n_nu)
    print np.average(np.arange(len(p_nu)), weights=p_nu)