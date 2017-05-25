__author__ = 'mjohnpayne'

# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

import plotly.plotly as py
from plotly.graph_objs import *
py.sign_in('mjohnpayne', 'u9h7zt53k3')

x = ['day 1', 'day 1', 'day 1', 'day 1', 'day 1', 'day 1',
     'day 2', 'day 2', 'day 2', 'day 2', 'day 2', 'day 2']

trace0 = Box(
    y=[0.2, 0.2, 0.6, 1.0, 0.5, 0.4, 0.2, 0.7, 0.9, 0.1, 0.5, 0.3],
    x=x,
    name='kale',
    marker=Marker(
        color='#3D9970'
    )
)
trace1 = Box(
    y=[0.6, 0.7, 0.3, 0.6, 0.0, 0.5, 0.7, 0.9, 0.5, 0.8, 0.7, 0.2],
    x=x,
    name='radishes',
    marker=Marker(
        color='#FF4136'
    )
)
trace2 = Box(
    y=[0.1, 0.3, 0.1, 0.9, 0.6, 0.6, 0.9, 1.0, 0.3, 0.6, 0.8, 0.5],
    x=x,
    name='carrots',
    marker=Marker(
        color='#FF851B'
    )
)
data = Data([trace0, trace1, trace2])
layout = Layout(
    yaxis=YAxis(
        title='normalized moisture',
        zeroline=False
    ),
    boxmode='group'
)
fig = Figure(data=data, layout=layout)
plot_url = py.plot(fig, filename='box-grouped')