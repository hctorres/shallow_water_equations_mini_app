# Shallow Water Equations Mini-App



$$ \psi = a \sin(x) \cos(y) $$

From the steam funcion we compute the velociy field. Is it a stream function? The value $a = 1000000$ is used.

$$ u = -\frac{\partial \psi}{\partial y} $$
$$ v = \frac{\partial \psi}{\partial x} $$

Note that this is different from the usual definiton of stream function that I have seen before. Where the signs are switched.

$$ u = \frac{\partial \psi}{\partial y} $$
$$ v = -\frac{\partial \psi}{\partial x} $$

Continuity is still automatically satisified either way but maybe look into this.

The pressure? is initialized to:
$$ p = C_{\mathrm{pcf}}( \cos(2 x) + \cos(2 y) ) + 5000 $$

where $C_{\mathrm{pcf}} = \frac{\pi^2 a^2}{e_l}$ where $e_l = n * dx$ but dx is hard set to a value of 100000 for some reason. Maybe the length scale is not 0 to $2\pi$ as associated with the initialization of the stream funcion but is instead $(0, N*dx)$ and $(0, N*dy)$. Where N and M are set according to the inputfile and dx = dy = 100000 is hard set.



