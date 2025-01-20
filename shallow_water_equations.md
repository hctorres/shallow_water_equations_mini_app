# Shallow Water Equations Mini-App


## Initializtion

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

where $C_{\mathrm{pcf}} = \frac{\pi^2 a^2}{e_l^2}$ where $e_l = n * dx$ but dx is hard set to a value of 100000 for some reason. Maybe the length scale is not 0 to $2\pi$ as associated with the initialization of the stream funcion but is instead $(0, N*dx)$ and $(0, N*dy)$. Where N and M are set according to the inputfile and dx = dy = 100000 is hard set.


## Governing Equation
The shallow water equations are:
$$ \frac{\partial{\vec{V}}}{\partial t} + \eta \hat{N} \times P \vec{V} + \nabla \left( P + \frac{1}{2}\vec{V} \cdot \vec{V} \right) = 0$$

$$ \frac{\partial{P}}{\partial t} + \nabla \cdot \left( P \vec{V} \right) = 0$$

where:

$$ \eta = \frac{\nabla \times \vec{V}}{P} $$

Where in component form we get:
$$ \frac{\partial u}{\partial t} = \eta P v - \frac{\partial}{\partial x} \left( P + \frac{1}{2} (u^2 + v^2) \right)$$

$$ \frac{\partial v}{\partial t} = -\eta P u - \frac{\partial}{\partial y} \left( P + \frac{1}{2} (u^2 + v^2) \right)$$

$$ \frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}\left(Pu\right) - \frac{\partial}{\partial y} (Pv) $$ 

$$ \eta = \frac{1}{P} \left( \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y} \right) = \frac{ \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y} }{P} $$

Note that the three equaitons are evaluated at different locaitons. The equation to time step u is evaluated at the y faces. Likewise the equation to time step v is evaluated at the x faces. The equation for eta is evaluated at the cell centers.

We need to evaluete the same term both the x and y velocity equations:
$$ h = P + \frac{1}{2} (u^2 + v^2) $$
Since we need $\frac{\partial h}{\partial x}$ at the y faces (location of $u$) and $\frac{\partial h}{\partial y}$ at the y faces (location of $v$). We will compute h at the nodal locations (same as P).

Taking keeping the indexing scheme in mind ()

$$ h_{i,j} = P_{i,j} + \frac{1}{2} \left( u^2|_{\mathrm{node}} + v^2|_{\mathrm{node}} \right) $$

$$ h_{i,j} = P_{i,j} + \frac{1}{2} \left( \frac{u_{i-1,j}^2 + u_{i,j}^2}{2} + \frac{v_{i,j-1}^2 + v_{i,j}^2}{2} \right) $$

$$ h_{i,j} = P_{i,j} + \frac{1}{4} \left( u_{i-1,j}^2 + u_{i,j}^2 + v_{i,j-1}^2 + v_{i,j}^2 \right) $$

In the code there is also a term for the pressure velcity product that is evaluated at the same location as the velocity. Called cu and cv in the code.

$$ [Pu] |_{\mathrm{y\_face}}  = [cu]_{i,j} = \frac{P_{i,j} + P_{i+1,j}}{2} u_{i,j}$$

$$ [Pv] |_{\mathrm{x\_face}} = \frac{P_{i,j} + P_{i,j+1}}{2} v_{i,j}$$

$\eta$ is called z in the code and evaluated at the cell centers. 

$$ \eta = \frac{ \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y} }{P} $$
$$ \eta_{i,j} = \frac{ \frac{v_{i+1,j} - v_{i,j}}{dx} + \frac{u_{i,j+1} - u_{i,j}}{dy} }{\frac{P_{i,j} + P_{i+1,j} + P_{i,j+1} + P_{i+1,j+1}}{4}} $$

$$ \eta_{i,j} = \frac{ \frac{4}{dx} (v_{i+1,j} - v_{i,j}) + \frac{4}{dy} (u_{i,j+1} - u_{i,j}) }{P_{i,j} + P_{i+1,j} + P_{i,j+1} + P_{i+1,j+1}} $$

Now the semi-descritized equation for the x velocity is:

$$ \frac{\partial u}{\partial t} = \eta P v - \frac{\partial}{\partial x} \left( P + \frac{1}{2} (u^2 + v^2) \right)$$

$$ \frac{\partial u}{\partial t} = \eta [P v] - \frac{\partial h}{\partial x} $$

$$ \frac{\partial u}{\partial t} = \left( \frac{\eta_{i,j-1} + \eta_{i,j}}{2} \right) \left( \frac{ [P v]_{i,j-1} + [P v]_{i,j} + [P v]_{i+1,j-1} + [P v]_{i+1,j} }{4} \right) - \left(\frac{h_{i+1,j} - h_{i,j}}{dx} \right)  $$

For the y velocity: 

$$ \frac{\partial v}{\partial t} = -\eta P u - \frac{\partial}{\partial y} \left( P + \frac{1}{2} (u^2 + v^2) \right)$$


$$ \frac{\partial v}{\partial t} = -\eta [P u] - \frac{\partial h}{\partial y}$$

$$ \frac{\partial v}{\partial t} = \left( \frac{\eta_{i-1,j} + \eta_{i,j}}{2} \right) \left( \frac{ [P u]_{i-1,j} + [P u]_{i-1,j+1} + [P u]_{i,j} + [P u]_{i,j+1} }{4} \right) - \left(\frac{h_{i,j+1} - h_{i,j}}{dy} \right)  $$

For the pressure:

$$ \frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}\left(Pu\right) - \frac{\partial}{\partial y} (Pv) $$ 

$$ \frac{\partial P_{i,j}}{\partial t} = - \left( \frac{ [Pu]_{i,j} - [Pu]_{i-1,j} }{dx} \right) - \left(  \frac{[Pv]_{i,j} - [Pu]_{i,j-1}}{dy} \right) $$