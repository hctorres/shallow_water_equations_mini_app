import yt
#from yt.frontends import boxlib
#from yt.frontends.boxlib.data_structures import AMReXDataset
from yt.frontends import amrex
from yt.frontends.amrex.data_structures import AMReXDataset

ds = AMReXDataset("plt00000")
#ds = AMReXDataset("plt01000")

ds.field_list

#sl = yt.SlicePlot(ds, 2, ('boxlib', 'phi'))
sl = yt.SlicePlot(ds, 'z', ('boxlib', 'phi'))

sl.save()
