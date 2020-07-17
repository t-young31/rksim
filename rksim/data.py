import os
import numpy as np
import rksim.exceptions as ex
from rksim.fit import fit
from rksim.plotting import plot


class TimeSeries:

    def _check(self):
        """Check the time and concentration arrays are the correct format"""
        # Arrays should be np.ndarrays and populated with floats
        for array in [self.times, self.concentrations]:
            assert type(array) is np.ndarray

            # Should be able to convert all values in the array to a float
            for value in array:
                float(value)

        # Times and concentrations need to be the same length to be plotted
        # against one another
        assert self.times.shape == self.concentrations.shape

        # Times should be monotonically increasing so the difference in
        # consecutive elements should all be positive
        assert np.min(np.diff(self.times)) > 0

    def __init__(self, name, times, concentrations, simulated=False):
        """
        Concentration of a reactant/product/intermediate as a function of time.

        :param name: (str) Name of this component in the reaction

        :param times: (np.ndarray) Monotonically increasing array of times in
                      seconds (s)

        :param concentrations: (np.ndarray) Concentration of this component
                               as a function of time in mol dm^-3. Shape should
                               be identical to times

        :param simulated: (bool) Is this time series simulated or experimental
        """
        self.name = str(name)
        self.is_simulated = simulated

        self.times = times
        self.concentrations = concentrations

        self._check()


class Data:

    def __add__(self, other):
        """Add another time-series or list of time series onto this dataset"""

        # Add a list or set of time series
        if type(other) is list or type(other) is set:
            assert all(isinstance(item, TimeSeries) for item in other)
            self._list += other

        # Add a time series
        if isinstance(other, TimeSeries):
            self._list.append(other)

        # Add another set of data
        if isinstance(other, Data):
            self._list += other._list

        return self

    def __getitem__(self, item):
        return self._list[item]

    def __iter__(self):
        return iter(self._list)

    def __len__(self):
        return len(self._list)

    def max_time(self):
        """Get the maximum time found in these data"""
        return max(series.times[-1] for series in self._list)

    def assign(self, system):
        """Assign time series to components in a system of reactions"""

        for time_series in self._list:
            for i in system.network.nodes:
                node = system.network.nodes[i]
                species = node['species']

                # If the name of the series is the same as the species
                # then assign it a time series
                if time_series.name == species.name:
                    species.series = time_series

                    # Initial concentration for this component
                    node['c0'] = time_series.concentrations[0]

        return None

    def fit(self, system, optimise=True):
        """Fit a system of reactions/equations to these data"""
        self.assign(system)
        return fit(self, system, optimise)

    def plot(self, name=None, dpi=400):
        """Plot the data with matplotlib"""
        return plot(self._list, name=name, dpi=dpi)

    def __init__(self, *args):
        """
        General data structure used for fitting. May be initialised from
        file formats including .csv

        :param args: (str) Name of the files to extract data from
        """
        self._list = []
        add_data_from_files(self, filenames=args)

        # Systems of equations (rksim.systems.System) fit for these data
        self.fits = None


def add_data_from_files(data, filenames):
    """From a list of files extract the data and add it"""

    for filename in filenames:
        data += extract_data(filename)

    return None


def extract_data(filename, **kwargs):
    """
    Extract data from a file. Expecting a format similar to the following .csv
    file:

    A title line
    t0, r0, p0
    t1, r2, p1
    .   .   .
    .   .   .

    where t, r and p are the time in seconds, concentration of reactant and
    product at that time in mol dm^-3. Supports any number of columns for
    different reactant components, but must have at least time and one
    concentration (i.e. no fewer than two columns).


    :param filename: (str) Name of the file to extract data from

    :optional names: (list(str)) Names of the time series in the file
              for example the above example would have names=['R', 'P']

    :optional name: (str) Name of the time series in the file

    :return: (list(rksim.data.TimeSeries)) List of time-series
    """
    # Ensure the file exists
    if not os.path.exists(filename):
        raise ex.RKSimCritical(f'Data file {filename} does not exist')

    array = get_raw_data_array(filename)
    n_rows, n_columns = array.shape

    # Data must have at least 2 columns
    if n_columns < 2:
        raise ex.DataMalformatted('Data must have at least time and one'
                                  'concentration column')

    # Times are the first column of the array
    times = array[:, 0]

    if n_columns == 2:
        name = filename if 'name' not in kwargs else kwargs['name']
        return TimeSeries(name, times, concentrations=array[:, 1])

    # There are more than 1 concentrations given as a function of time
    series_list = []

    # Iterate through the columns adding a time series for that component
    for i in range(1, n_columns):

        # Name can be defined in the keyword arguments
        if 'names' in kwargs:
            name = kwargs['names'][i-1]

        else:
            name = f'{filename}_{i}'

        series = TimeSeries(name, times, concentrations=array[:, i])
        series_list.append(series)

    return series_list


def get_raw_data_array(filename):
    """Get a data array (matrix) from a file"""

    data_file = open(filename, 'r')

    if filename.endswith('.csv'):
        array = get_raw_data_from_csv(csv_file=data_file)

    else:
        raise NotImplementedError

    data_file.close()

    # Ensure the array is regular..
    first_row_length = len(array[0])

    # All the rows must be the same length
    if not all(len(row) == first_row_length for row in array):
        raise ex.DataMalformatted

    return array


def get_raw_data_from_csv(csv_file):
    """
    Extract the data from a .csv file into a numpy array

    :param csv_file: (file)
    :return: (np.ndarray)
    """
    array = []

    for line in csv_file:
        items = line.split(',')

        # Try and convert all the items in this row of the csv to a float,
        # continue if this is not possible
        try:
            # Remove any whitespace from the data
            row = [float(item.strip()) for item in items]
            array.append(row)

        except ValueError:
            continue

    if len(array) == 0:
        raise ex.DataMalformatted

    return np.array(array)
