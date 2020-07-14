import numpy as np


def generate_1st_order_data():
    """Generate data for a first order reaction: R -> P

    d[P]/dt = k [R]            :   [P] = [R]0(1 - e^(-kt))
    d[R]/dt = -k [R]           :   [R] = [R]0 - [P]

    k    = 1 s^-1
    [R]0 = 1 mol dm^-3
    [P]0 = 0

    Data printed in a .csv file in the format time, [P], [R]
    """
    times = np.linspace(0, 10, num=100)            # s
    r0 = 1                                          # mol dm^-3
    k = 1                                           # s^-1

    with open('first_order.csv', 'w') as data_file:
        print('Data for R -> P where v = k[R]', file=data_file)

        for t in times:
            p = r0 * (1.0 - np.exp(-k*t))
            r = r0 - p

            #       Time,     [P],    [R]
            print(f'{t:.6f},{p:.6f},{r:.6f}', file=data_file)

    return None


if __name__ == '__main__':

    generate_1st_order_data()
