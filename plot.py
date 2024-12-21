import yt
#from yt.frontends import boxlib
#from yt.frontends.boxlib.data_structures import AMReXDataset
from yt.frontends import amrex
from yt.frontends.amrex.data_structures import AMReXDataset

#yt.toggle_interactivity()

ds = AMReXDataset("plt00000")
#ds = AMReXDataset("plt01000")

ds.field_list

#sl = yt.SlicePlot(ds, 'z', ('boxlib', 'phi'), origin="native", axes_unit="unitary")

#sl = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'phi'), origin="native", axes_unit="unitary")
sl = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'p'), origin="native", axes_unit="unitary")

sl.save()
