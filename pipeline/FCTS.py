import numpy as np

# Forward Time Central Space

# ht <= 1/(2*alfa) * (1/hx**2 + 1/hy**2)**(-1)

# u_{n+1} = u_n + alfa * ( ht/ hx**2 * (u{i+1,j} - 2*u{i,j} + u{i-1,j}) + ht/ hy**2 * (u{i,j+1} - 2*u{i,j} + u{i,j-1})) +
# + ht * f{i,j,n}

#if
#print(ht <= 1/(2*alfa) * (1/hx**2 + 1/hy**2)**(-1))


lx = alfa * ht/hx**2
ly = alfa * ht/hy**2


for i in range(t):




