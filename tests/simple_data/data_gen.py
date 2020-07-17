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


def generate_reversible_1st_order_data():
    """Generate data for a first order reaction: R <-> P as
    generate_1st_order_data()

    [R] = ([R]0 / (kf + kb)) (kb + kfe^(-t(kf + kb)))
    [P] = [R]0 - R

    where [P]0 = 0

    Data printed in a .csv file in the format time, [P], [R]
    """
    times = np.linspace(0, 10, num=50)                # s
    r0 = 1                                            # mol dm^-3
    kf = 0.5                                          # s^-1
    kb = 0.1                                          # s^-1

    with open('first_order_reversible.csv', 'w') as data_file:
        print('Data for R <-> P where v = k[R]', file=data_file)

        for t in times:
            r = (r0 / (kf + kb)) * (kb + kf * np.exp(-t * (kf + kb)))
            p = r0 - r

            #       Time,     [P],    [R]
            print(f'{t:.6f},{p:.6f},{r:.6f}', file=data_file)

    return None


def generate_two_1st_order_data():
    """
    Generate data for two first order reactions. A -> B and A -> C i.e.
    the network: B <-- A --> C
                    k1   k2

    From https://en.wikipedia.org/wiki/Rate_equation the integrated rate
    equations are

    [A] = [A]0 e^-t(k1 + K2)
    [B] = (k1 / (k1 + k2)) [A]0 (1 - e^(-t(k1+k2)))
    [C] = ((k2 / (k1 + k2)) [A]0 (1 - e^(-t(k1+k2)))

    where [B]0 = [C]0 = 0
    """
    times = np.linspace(0, 10, num=50)                # s

    a0 = 2.0                                          # mol dm^-3
    k1 = 0.4                                          # s^-1
    k2 = 0.8                                          # s^-1

    with open('double_first_order.csv', 'w') as data_file:
        print('Data for B <-- A --> C', file=data_file)

        for t in times:
            a = a0 * np.exp(-t * (k1 + k2))
            b = a0 * (k1 / (k2 + k1)) * (1.0 - np.exp(-t * (k1 + k2)))
            c = a0 * (k2 / (k2 + k1)) * (1.0 - np.exp(-t * (k1 + k2)))

            #       Time,   [A],   [B],   [C]
            print(f'{t:.6f},{a:.6f},{b:.6f},{c:.6f}', file=data_file)

    return None


def generate_consecutive_data():
    """
    Generate data for A --> B --> C
                         k1   k2

    d[A]/dt = -k1 [A]          :   [A] = [A]0 e^(-t k1)

    d[B]/dt = k1[A] - k2[B]    :
    [B] = [A]0 (k1 / (k2 - k1)) (e^(-t k1) - e^(-t k2)) + [B]0e^(-t k2)

    d[C]/dt = k2[B]            :
    [C] = [A]0 (1 + ((k1e^(-t k2) - k2e^(-t k1)) / (k2 - k1)) +
          [B]0(1 - e^(-t k2)) + [C]0
    """
    times = np.linspace(0, 10, num=50)  # s

    a0 = 2.0  # mol dm^-3
    b0 = 1.7  # mol dm^-3
    c0 = 0.1  # mol dm^-3

    k1 = 0.5  # s^-1
    k2 = 0.8  # s^-1

    with open('consecutive_first_order.csv', 'w') as data_file:
        print('Data for A --> B --> C', file=data_file)

        for t in times:
            exp1 = np.exp(-t * k1)
            exp2 = np.exp(-t * k2)

            a = a0 * exp1
            b = a0 * (k1 / (k2 - k1)) * (exp1 - exp2) + b0 * exp2

            frac = (k1 * exp2 - k2 * exp1) / (k2 - k1)
            c = a0 * (1.0 + frac) + b0 * (1 - exp2) + c0

            #       Time,   [A],   [B],   [C]
            print(f'{t:.6f},{a:.6f},{b:.6f},{c:.6f}', file=data_file)

    return None


if __name__ == '__main__':

    generate_1st_order_data()
    generate_reversible_1st_order_data()
    generate_two_1st_order_data()
    generate_consecutive_data()
