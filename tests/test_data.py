from rksim.data import Data, TimeSeries
from rksim.systems import *
from rksim.data import extract_data
from rksim.exceptions import RKSimCritical, DataMalformatted
from rksim.fit import get_initial_concentrations
import pytest
import numpy as np
import os

here = os.path.abspath(os.path.dirname(__file__))


def test_time_series():
    times = np.linspace(0, 1, 3)
    concs = np.array([0.0, 1.0, 0.0])

    with pytest.raises(AssertionError):
        _ = TimeSeries(name='test', times=None, concentrations=None)
        _ = TimeSeries(name='test', times=times, concentrations=None)
        _ = TimeSeries(name='test', times=None, concentrations=concs)

        _ = TimeSeries(name='test', times=[0, 1, 2], concentrations=concs)
        _ = TimeSeries(name='test', times=times,
                       concentrations=[0.0, 1.0, 0.0])

        # Cannot have decreasing times
        _ = TimeSeries(name='test',
                       times=np.linspace(0, -1, 3),
                       concentrations=concs)

        # Cannot set arrays of strings for times or concentrations
        _ = TimeSeries(name='test', times=np.array(['a', 'b', 'c']),
                       concentrations=concs)

    ts = TimeSeries(name='test', times=times, concentrations=concs)
    assert ts.times.shape == (3,)
    assert ts.concentrations.shape == (3,)


def test_data():

    with pytest.raises(RKSimCritical):
        _ = Data('a_filename_that_doesnt_exist')

    with pytest.raises(DataMalformatted):
        data_path = os.path.join(here, 'data', 'first_order_broken.csv')
        _ = Data(data_path)

    # First order kinetic data
    data_path = os.path.join(here, 'data', 'first_order.csv')
    data = Data(data_path)

    # Test plotting the data
    data.plot(name='test.png', dpi=100)
    assert os.path.exists('test.png')
    os.remove('test.png')


def test_data_with_names():

    # Two column data of R and P
    data_path = os.path.join(here, 'data', 'first_order.csv')

    data = Data()
    data += extract_data(data_path, names=['R', 'P'])

    assert data[0].name == 'R'
    assert data[1].name == 'P'

    # Single column data
    data_path = os.path.join(here, 'data', 'first_order_only_product.csv')

    data = Data()
    series = extract_data(data_path, name='P')

    data += series

    assert data[0].name == 'P'


def test_simple_fit():
    data_path = os.path.join(here, 'data', 'first_order.csv')

    data = Data()
    # Extract the data in the correct order...!
    data += extract_data(data_path, names=['P', 'R'])

    assert 9.9 < data.get_max_time() < 10.1

    system = System(Irreversible(Reactant(name='R'),
                                 Product(name='P')))
    data.assign_series(system)

    r0, p0 = get_initial_concentrations(system)
    assert 0.999 < r0 < 1.001
    assert -0.001 < p0 < 0.001

    data.fit(system)
    system.plot()
    assert os.path.exists('system.png')
    os.remove('system.png')