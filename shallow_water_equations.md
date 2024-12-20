# Shallow Water Equations Mini-App



$$ \psi = a \sin(x) \cos(y) $$
$$ p = C_{pcf}( \cos(2 x) + \cos(2 y) ) + 5000 $$

From the steam funcion we compute the velociy field.

$$ u = -\frac{\partial \psi}{\partial y} $$
$$ v = \frac{\partial \psi}{\partial x} $$

Note that this is different from the usual definiton of stream function that I have seen before. Where the signs are switched.

$$ u = \frac{\partial \psi}{\partial y} $$
$$ v = -\frac{\partial \psi}{\partial x} $$

Continuity is still automatically satisified either way but maybe look into this.