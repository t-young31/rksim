from rksim import *

data = Data()
# Add data from a csv file containing a time series of [R]
# where the first column is the time (s) and the second
# the concentration of R (mol dm^-3) at that time
data += extract_data('example_data.csv', name='R')

# Form the system of reactions. Here only R -> P
reaction = Irreversible(Reactant('R'), Product('P'))
system = System(reaction)

# Fit the system of reactions to the data
data.fit(system)

# Plot the simulated components and the rate constant
system.plot(name='example')
print('k_R->P', system.rate_constant('R', 'P'))
