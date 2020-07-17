import matplotlib as mpl
from rksim.data import TimeSeries
from rksim.plotting import plot
from rksim.exceptions import PlottingFailed
import pytest
import numpy as np
import os

mpl.use('Template')


def test_plot_series():

    series = TimeSeries(name='test',
                        concentrations=np.linspace(0, 1, 10),
                        times=np.linspace(0, 1, 10))

    # Should be able to plot a time series with no name (and use plt.show())
    plot([series])

    # If a name is specified then it should be saved to the current dir
    plot([series], name='test')
    assert os.path.exists('test.png')
    os.remove('test.png')

    with pytest.raises(PlottingFailed):
        plot([None])
