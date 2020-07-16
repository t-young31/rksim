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


if __name__ == '__main__':

    generate_2nd_order_data()
