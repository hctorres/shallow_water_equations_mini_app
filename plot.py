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
sl_phi_old = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'phi_old'), origin="native", axes_unit="unitary")
sl_p = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'p'), origin="native", axes_unit="unitary")
sl_u = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'u'), origin="native", axes_unit="unitary")
sl_v = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'v'), origin="native", axes_unit="unitary")

#sl.plots[("boxlib", "output_array")].axes.grid(True)

sl_phi_old.save()
sl_p.save()
sl_u.save()
sl_v.save()