# -*- coding: utf-8 -*-
"""
Created on Wed May 17 04:28:47 2023

@author: SSI-User
"""
# Import packages
from dash import Dash, html, dash_table
import pandas as pd

# Incorporate data
df = pd.read_csv('https://raw.githubusercontent.com/DrucillaLiu/potential-distribution/main/potentials.csv')

# Initialize the app
app = Dash(__name__)

# App layout
app.layout = html.Div([
    html.Div(children='My First App with Data'),
    dash_table.DataTable(data=df.to_dict('records'),page_size=13)
])

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=False)