import h5py
import numpy as np
from yaff.analysis.utils import get_slice

f = h5py.File('restart_1.h5',mode='a')
tgrp = f['trajectory']

old_counter = np.array(tgrp['counter'])
shape = old_counter.shape
tgrp.create_dataset('econs_counter', shape, data=old_counter+1)
f.close()
