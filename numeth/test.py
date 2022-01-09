from numeth.Euler import euler 

# a, b, N, w = 1, 2, 5, 1
# test = euler(a, b, N, w).Method()

from numeth.Taylor import taylor

# a, b, N, w = 1, 2, 10, 0
# test = taylor(a, b, N, w).Method(deg=4)

from numeth.Midpoint import midpoint

# a, b, N, w = 0, 1, 2, 0
# test = midpoint(a, b, N, w).Method()
from numeth.Mod_Euler import mod_euler

a, b, N, w = 0, 1, 1, 10
test = mod_euler(a, b, w, N).Method()