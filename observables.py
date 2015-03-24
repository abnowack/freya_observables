"""
TODO:
Might consider not saving in event folder, but with event field sequentially
[x] Finalize HDF5 Format
[ ] Automated ways to plot
    [x] Avg Neutrons
    [x] Avg Photons
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
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('figure', fc='white', dpi=150)
rc('lines', lw=2, color='black')
rc('axes', color_cycle=['black'])

def parse_ptrac_to_hdf5(ptrac_filename, hdf5_filename='ptrac.h5', 
                        ev_buffer_len=1000, h_buffer_len=1000):
    ''' save neutron and photon position, directory, energy, time '''
    with open(ptrac_filename, 'r') as ptrac:
        # parse headers and formats
        header = preader.ptrac_header(ptrac)
        input_format = preader.ptrac_input_format(ptrac)
        event_format = preader.ptrac_event_format(ptrac)

        with h5.File(hdf5_filename, 'w') as f:
            dt = np.dtype([('evt', np.int), ('ipt', np.int), ('x', np.float),
                           ('y', np.float), ('z', np.float), ('u', np.float),
                           ('v', np.float), ('w', np.float),
                           ('energy', np.float), ('time', np.float)])
            ref_dt = h5.special_dtype(ref=h5.RegionReference)

            f_history_d = f.create_dataset('history', (h_buffer_len,), dtype=ref_dt, maxshape=(None,), chunks=True)
            f_events_d = f.create_dataset('events', (ev_buffer_len,), dtype=dt, maxshape=(None,), chunks=True)
            nevt = 0
            ev_buffer_count = 1
            h_buffer_count = 1

            for i, history in enumerate(preader.parse_ptrac_events(ptrac, event_format)):
                if (ev_buffer_len * ev_buffer_count - nevt) < len(history.events):
                    ev_buffer_count += 1
                    f_events_d.resize((ev_buffer_len * ev_buffer_count,))

                if (h_buffer_len * h_buffer_count) <= i:
                    h_buffer_len += 1
                    f_history_d.resize((h_buffer_len * h_buffer_count,))
                
                f_history_d[i] = f_events_d.regionref[nevt:nevt+len(history.events)]

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
            f_history_d.resize((i,))

if __name__ == '__main__':

    with open('input_cards/geometry.i', 'r') as cardfile:
        card = cardfile.read()

    with run_mcnp(card, cores=4) as (status, mcnp_dir):
        parse_ptrac_to_hdf5(mcnp_dir + '\\ptrac')

    # calculate nu
    print 'processing'
    with h5.File('ptrac.h5', 'r') as f:
        hist = f['history']
        evts = f['events']
        
        # avg neutron, photon
        n_n = np.zeros(np.size(hist))
        n_p = np.zeros(np.size(hist))
        
        t_n, t_p = [], []
        
        dt_n, dt_p = [], []

        for i, h in enumerate(hist):
            evt = evts[h]
            neutrons = evt[evt['ipt'] == 1]
            photons = evt[evt['ipt'] == 2]
            
            n_n[i] = neutrons.size
            n_p[i] = photons.size
            
            t_n.extend(neutrons['time'])
            t_p.extend(photons['time'])
            
            if n_n[i] > 1:
                for n1 in neutrons:
                    for n2 in neutrons[1:]:
                        dt_n.append(n1['time'] - n2['time'])
            
            if n_p[i] > 1:
                for p1 in photons:
                    for p2 in photons[1:]:
                        dt_p.append(p1['time'] - p2['time'])

    fig = plt.figure()
    plt.hist(n_n, bins=30, range=[0, 30], label='Neutron', histtype='step', linestyle='dashed')
    plt.hist(n_p, bins=30, range=[0, 30], label='Photon', histtype='step', linestyle='solid')
    plt.yscale('log')
    plt.xlabel('Frequency')
    plt.ylabel('Count')
    plt.title('# Escaped Per Event')
    plt.legend([plt.Line2D([0], [0], color='black', linestyle='dashed'),
                plt.Line2D([0], [0], color='black', linestyle='solid')], 
               [r'Neutron $\nu$ = ' + '{0:.2f}'.format(np.average(n_n)), 
                r'Photon  $\nu$ = ' + '{0:.2f}'.format(np.average(n_p))], frameon=False)
    plt.show()
    
    fig = plt.figure()
    plt.hist(t_n, bins=50, range=[0, 50], histtype='step', label='Neutron', linestyle='dashed')
    plt.hist(t_p, bins=50, range=[0, 50], histtype='step', label='Photon', linestyle='solid')
    plt.legend([plt.Line2D([0], [0], color='black', linestyle='dashed'),
                plt.Line2D([0], [0], color='black', linestyle='solid')], 
               [r'Neutron', r'Photon'], frameon=False) 
    plt.xlabel('Time Since Fission (ns)')
    plt.ylabel('Count')
    plt.title('Detection Time Distribution')
    plt.show()
    
    fig = plt.figure()
    plt.hist(dt_n, bins=100, range=[-50, 50], histtype='step', label='Neutron', linestyle='dashed')
    plt.hist(dt_p, bins=100, range=[-50, 50], histtype='step', label='Photon', linestyle='solid')
    plt.legend([plt.Line2D([0], [0], color='black', linestyle='dashed'),
                plt.Line2D([0], [0], color='black', linestyle='solid')], 
               [r'Neutron', r'Photon'], frameon=False) 
    plt.xlabel('Time Between Events (ns)')
    plt.ylabel('Count')
    plt.title('Timing Correlation')
    plt.show()