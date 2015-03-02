"""
TODO:
[ ] FIX mcnpy/ptrac/reader.py input_format indexing issue
[ ] Extract terminated neutrons, record position, velocity vectors and time
[ ] Do same for photons
[ ] Create method for storing neutron and photon lists for each event
"""

import numpy as np
import matplotlib.pyplot as plt
from mcnpy.mcnp_wrapper import *
from mcnpy.ptrac import reader as preader

def parse_ptrac_events(ptrac_filename):
    positions = []
    
    with open(ptrac_filename, 'r') as ptrac:
        # parse headers and formats
        header = preader.ptrac_header(ptrac)
        input_format = preader.ptrac_input_format(ptrac)
        event_format = preader.ptrac_event_format(ptrac)
        
        for history in preader.parse_ptrac_events(ptrac, event_format):
            for ev in history.events:
                print ev
  
with open('input_cards/geometry.i', 'r') as cardfile:
    card = cardfile.read()

with run_mcnp(card, cores=4) as (status, mcnp_dir):
    parse_ptrac_events(mcnp_dir + '\\ptrac')