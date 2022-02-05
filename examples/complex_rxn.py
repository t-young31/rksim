from rksim import *


if __name__ == '__main__':

    system = System(Reaction('a+A->b'),
                    Reaction('a+B->d'),
                    Reaction('b->C+a'),
                    Reaction('d+C->E+a'),
                    Reaction('d+A->D+a'),
                    Reaction('a+D->c'),
                    Reaction('c->F+a'))

    system.set_initial_concentration('A', 0.05)
    system.set_initial_concentration('B', 0.05)
    system.set_initial_concentration('a', 0.05*2/100)  # 2 mol %

    system.set_rate_constants(ks=[100, 100, 0.027, 0.0, 0.1, 100, 0.027])
    system.fit(data=None, max_time=30*60*60)  # 30h with no true data

    # Save all the time series for further analysis
    for species in system.species:
        species.simulated_series.save()
