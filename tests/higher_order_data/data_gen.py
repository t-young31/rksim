"""
Integrated rate equations from
https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118468586.epoc1012
which references

R. Livingstone, Evaluation and Interpretation of Rate Data, inInvestigation of
Rates and Mechanisms of Reactions(S. L. Friess, E. S. Lewis, and A.Weissberger,
Eds.), 2nd ed., Part I, Interscience, New York, 1961, p. 114
"""

import numpy as np


def generate_2nd_order_data():
    """Generate data for a 2nd order reaction: 2R -> P

    d[R]/dt = -2k [R]^2             :   1/[R] = 1/[R}0 + 2kt
    d[P]/dt = k [R]^2              :   [P] = [R]0 - 2x[R]

    k    = 1 s^-1
    [R]0 = 1 mol dm^-3
    [P]0 = 0

    Data printed in a .csv file in the format time, [P], [R]
    """
    times = np.linspace(0, 10, num=100)            # s
    r0 = 0.6                                       # mol dm^-3
    k = 2                                          # s^-1

    with open('second_order.csv', 'w') as data_file:
        print('Data for R -> P where v = k[R]^2', file=data_file)

        for t in times:
            r = (1.0 / r0 + 2.0 * k * t)**-1
            p = (r0 - r) / 2

            #       Time,     [R],    [P]
            print(f'{t:.6f},{r:.6f},{p:.6f}', file=data_file)

    return None


def generate_ab_data():
    """
    Generate data for a second order reaction A + B -> P

    d[A]/dt = -k[A][B]
    d[B]/dt = -k[A][B]
    d[P]/dt = k[A][B]

    [P] = ([B]0 - [A]0 h(t)) / (1 - h(t))  where
    h(t) = ([B]0 / [A]0) e^(kt ([B]0 - [A]0))

    Data printed in a .csv file
    """
    times = np.linspace(0, 10, num=100)            # s
    a0 = 0.6                                       # mol dm^-3
    b0 = 0.5                                       # mol dm^-3

    k = 1.7                                        # mol^-1 dm^3 s^-1

    with open('ab.csv', 'w') as data_file:
        print('Data for A + B -> P where v = k[A][B]', file=data_file)

        for i, t in enumerate(times):

            h = (b0 / a0) * np.exp(k * t * (b0 - a0))
            p = (b0 - a0 * h) / (1.0 - h)

            a = a0 - p
            b = b0 - p

            #       Time,     [A],    [B],   [P]
            print(f'{t:.6f},{a:.6f},{b:.6f},{p:.6f}', file=data_file)

    return None


if __name__ == '__main__':

    generate_2nd_order_data()
    generate_ab_data()
