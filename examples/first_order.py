from rksim import *

data = Data()
data += extract_data('first_order.csv', names=['P', 'R'])

reaction = Irreversible(Reactant('R'), Product('P'))
system = System(reaction)

data.fit(system)
system.plot(name='first_order')
print('k_R->P', system.rate_constant('R', 'P'))
