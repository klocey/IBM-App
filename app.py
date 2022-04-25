import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
from dash import dash_table
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go

import math
import os
import numpy as np
import pandas as pd
import scipy as sc
from scipy import stats
from scipy.stats import gaussian_kde

from random import randint, choice, shuffle, sample, seed
import warnings
import sys
import re
import time
import base64
import io
import json
import ast
import time
from functools import reduce
import operator

#########################################################################################
################################# CONFIG APP ############################################
#########################################################################################

FONT_AWESOME = "https://use.fontawesome.com/releases/v5.10.2/css/all.css"

warnings.filterwarnings('ignore')
pd.set_option('display.max_columns', None)

external_stylesheets=[dbc.themes.BOOTSTRAP, FONT_AWESOME]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.config.suppress_callback_exceptions = True
server = app.server


#########################################################################################
########################### CUSTOM FUNCTIONS ############################################
#########################################################################################

def get_kdens_choose_kernel(_list, kernel=0.5):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    #xs = np.linspace(0, 1, n)
    xs = np.linspace(min(_list), max(_list), n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D
    
#########################################################################################
#################### DASH APP CONTROL CARDS  ############################################
#########################################################################################


def description_card1():
    """
    :return: A Div containing dashboard title & descriptions.
    """
    return html.Div(
        id="description-card1",
        children=[
            html.H5("IBM Workbench: Explore relationships of metabolic activity and biodiversity using individual-based models (IBMs)",
                    style={
            'textAlign': 'center',
            }),
        ],
    )


def control_card1():
    return html.Div(
        id="control-card1",
        children=[
            html.H5("Customizable parameters", style={'display': 'inline-block',
            'width': '80%'},),
            html.I(className="fas fa-question-circle fa-lg", id="target1",
                style={'display': 'inline-block', 'width': '10%', 'color':'#99ccff'},
                ),
            dbc.Tooltip("Change these while an IBM is running, and the IBM will automatically respond.", target="target1",
                style = {'font-size': 12},
                ),
            html.Hr(),
            #html.Br(),
            html.Div(
                id="Number of species",
                children=[
                        html.P('Species', style={'display': 'inline-block',
                                                 'font-size': 17,
                                                 'width': '51%'},
                                                 ),
                        html.I(className="fas fa-question-circle fa-lg", id="target2",
                            style={'display': 'inline-block', 'width': '20%', 'color':'#cccccc'},
                            ),
                        dbc.Tooltip("The number of species the model starts with. If immigration is allowed, this is the number of species in the metacommunity. Species traits are randomly parameterized. Starting with a greater number of species means starting with greater diversity of trait combinations.", target="target2",
                            style = {'font-size': 12},
                            ),
                        dcc.Input(id='S',
                            type='number',
                            value=100,
                            min=1, max=1000, step=1),
                        ],
                    style={ 'width': '40%',
                            'display': 'inline-block',
                            'margin-right': '20px',
                    },
                ),
            html.Div(
                id="Resource inflow",
                children=[
                        html.P('Resource inflow', style={'display': 'inline-block',
                                                       'font-size': 17,
                                                       'width': '80%'},
                                                      ),
                        html.I(className="fas fa-question-circle fa-lg", id="target4",
                            style={'display': 'inline-block', 'width': '20%', 'color':'#cccccc'},
                            ),
                        dbc.Tooltip("Volume of resources per unit influent.", target="target4",
                            style = {'font-size': 12},
                            ),
                            
                        dcc.Input(id='R',
                            type='number',
                            value=100,
                            min=0, max=1000, step=1),
                        ],
                    style={'width': '50%',
                            'display': 'inline-block',
                            #'margin-right': '13px',
                    },
                ),
            
            html.Hr(),
            html.Br(),
            
            html.Div(
            id="Flow rate",
            children=[
                    html.P('Flow rate', style={'display': 'inline-block',
                                             'font-size': 17,
                                             'width': '52%'},
                                             ),
                    html.I(className="fas fa-question-circle fa-lg", id="target3",
                        style={'display': 'inline-block', 'width': '20%', 'color':'#cccccc'},
                        ),
                    dbc.Tooltip("Units of volume flowing in (influent) and out per time step. Each unit is 1/100 of the system volume.", target="target3",
                        style = {'font-size': 12},
                        ),
                    dcc.Input(id='Q',
                        type='number',
                        value=5,
                        min=1, max=20, step=1),
                    ],
                style={ 'width': '45%',
                        'display': 'inline-block',
                        #'margin-right': '20px',
                },
            ),
            
            html.Div(
            id="Immigration",
            children=[
                    html.P('Immigration rate', style={'display': 'inline-block',
                                                   'font-size': 17,
                                                   'width': '80%'},
                                                  ),
                    html.I(className="fas fa-question-circle fa-lg", id="target5",
                        style={'display': 'inline-block', 'width': '20%', 'color':'#cccccc'},
                        ),
                    dbc.Tooltip("Individuals per unit influent. Immigrants are chosen at random according to species-specific rates.", target="target5",
                        style = {'font-size': 12},
                        ),
                        
                    dcc.Input(id='immigration',
                        type='number',
                        value=1,
                        min=0, max=10, step=1),
                    ],
                style={'width': '50%',
                        'display': 'inline-block',
                        #'margin-right': '13px',
                },
            ),
            
            html.Hr(),
            #html.Br(),
            
            html.Div(
                id="Immigration on/off",
                    children=[
                            html.P('Immigration', style={'display': 'inline-block',
                                                           'font-size': 17,
                                                           'width': '80%'},
                                                          ),
                            dcc.RadioItems(id = 'immigration_on_off',
                                    options=[{"label": i, "value": i} for i in [' on', ' off']],
                                    value=' on',
                                    ),
                            ],
                    style={'width': '45%',
                            'display': 'inline-block',
                            #'margin-right': '13px',
                    },
                ),
            html.Div(
                id="Reproduction on/off",
                    children=[
                            html.P('Reproduction', style={'display': 'inline-block',
                                                           'font-size': 17,
                                                           'width': '80%'},
                                                          ),
                            dcc.RadioItems(id = 'reproduction_on_off',
                                    options=[{"label": i, "value": i} for i in [' on', ' off']],
                                    value=' on',
                                    ),
                            ],
                    style={'width': '45%',
                            'display': 'inline-block',
                            #'margin-right': '13px',
                    },
                ),
            
            html.Hr(),
            html.Div(
                id="Active dispersal on/off",
                    children=[
                            html.P('Active disperal', style={'display': 'inline-block',
                                                           'font-size': 17,
                                                           'width': '80%'},
                                                          ),
                            dcc.RadioItems(id = 'active_dispersal_on_off',
                                    options=[{"label": i, "value": i} for i in [' on', ' off']],
                                    value=' on',
                                    ),
                            ],
                    style={'width': '45%',
                            'display': 'inline-block',
                            #'margin-right': '13px',
                    },
                ),
            html.Div(
            id="Death on/off",
                children=[
                        html.P('Death', style={'display': 'inline-block',
                                                       'font-size': 17,
                                                       'width': '80%'},
                                                      ),
                        dcc.RadioItems(id = 'death_on_off',
                                options=[{"label": i, "value": i} for i in [' on', ' off']],
                                value=' on',
                                ),
                        ],
                style={'width': '45%',
                        'display': 'inline-block',
                        #'margin-right': '13px',
                },
            ),
            
            html.Hr(),
            html.Br(),
            
            html.H5("IBM controls", style={'width': '42%',
                    'display': 'inline-block',
                    #'margin-right': '13px',
                    },),
            html.I(className="fas fa-question-circle fa-lg", id="target_ibm_controls",
                style={'display': 'inline-block', 'width': '20%', 'color':'#99ccff'},
                ),
            dbc.Tooltip("IBMs will run more slowly with several thousands of individuals. You can rarefy to 1,000 randomly chosen individuals to stop run-away growth. You should parameterize a system that fluctuates below 5K individuals or at no more than 10K.", target="target_ibm_controls",
                style = {'font-size': 12},
                ),
                
            
            html.Button('Run new IBM', id='btn1', n_clicks=0,
                style={#'width': '100%',
                        'display': 'inline-block',
                        'margin-right': '2%',
                        'margin-left': '3%',
                },
                ),
            html.Button('Pause/Run', id='btn2', n_clicks=0,
            style={'width': '45%',
                    'display': 'inline-block',
            },
            ),
            html.Hr(),
            html.Button('Clear/Reset', id='btn3', n_clicks=0,
            style={#'width': '90%',
                    'display': 'inline-block',
                    'margin-right': '2%',
                    'margin-left': '3%',
                    
            },
            ),
            html.Button('Rarefy to 1K', id='btn-rarefy', n_clicks=0,
            style={'width': '45%',
                    'display': 'inline-block',
                    'margin-right': '2%',
                    #'margin-left': '3%',
            },
            ),
            ],
        )

def time_series1():
    return html.Div(id="right-column1", className="one columns",
            children=[
                        html.Div(
                            id="time_series_options",
                            children=[
                            html.H5("Plot a time series"),
                            dcc.Dropdown(
                                id='plot_by2',
                                options=[{"label": i, "value": i} for i in ['Total abundance (N)', 'Species richness (S)', 'Total resources']],
                                value=None,
                                style={'display': 'inline-block',
                                    'width': '80%',
                                    'margin-right': '0px',
                                    'font-size': "100%"},
                                ),
                            html.Button('Plot', id='btn4', #n_clicks=0,
                            style={#'width': '10%',
                                    'display': 'inline-block',
                                    'margin-right': '0px',
                                    'margin-left': '0px',
                                    #'horizontal-align': 'left',
                                    'vertical-align': 'top',
                                },
                            ),
                            ],style={'display': 'inline-block', 'width':'75%'},
                            ),
                        html.Hr(),
                        dcc.Graph(id='time_series_fig', style={'width': '100%',
                            #'display': 'inline-block',
                            'background-color': '#f0f0f0','padding': '0px', 'margin-bottom': '0px',
                            'margin-right': '0px','margin-left': '0px','height': '400px',
                            },),
                        ],
                    style={'width': '100%',
                        #'display': 'inline-block',
                        'background-color': '#f0f0f0','padding': '0px', 'margin-bottom': '0px',
                        'margin-right': '0px','margin-left': '0px','height': '500px',
                        },
                  )
                  

def distribution_1():
    return html.Div(id="right-column3", className="one columns",
    children=[
                html.Div(
                    id="distribution_options",
                    children=[
                    html.H5("Plot a distribution"),
                    dcc.Dropdown(
                        id='plot_by3',
                        options=[{"label": i, "value": i} for i in ['growth rate', 'active dispersal rate', 'resuscitation rate', 'basal metabolic rate', 'bmr reduction in dormancy', 'immigration rate', 'resource quota', 'body size']],
                        value=None,
                        style={'display': 'inline-block',
                            'width': '80%',
                            'font-size': "100%"},
                        ),
                    html.Button('Plot', id='btn5', n_clicks=0,
                    style={#'width': '10%',
                            'display': 'inline-block',
                            #'margin-right': '2%',
                            #'margin-left': '-140px',
                            'horizontal-align': 'left',
                            'vertical-align': 'top',
                        },
                    ),
                    ],style={'display': 'inline-block', 'width': '75%'},
                    ),
                html.Hr(),
                dcc.Graph(id='distribution_fig', style={'width': '100%',
                #'display': 'inline-block',
                'background-color': '#f0f0f0','padding': '0px', 'margin-bottom': '0px',
                'margin-right': '0px','margin-left': '0px','height': '400px',
                },),],
            style={'width': '100%',
                #'display': 'inline-block',
                'background-color': '#f0f0f0','padding': '0px', 'margin-bottom': '0px',
                'margin-right': '0px','margin-left': '0px','height': '500px',
                },
          )
          
def xy_1():
    return html.Div(id="right-column4", className="one columns",
    children=[html.Div(id="xy_options",
                    children=[
                    html.H5("Explore x-y relationships"),
                    html.Div(id="x_options",
                        children=[
                        html.B("X-variable", style={'font-size': "100%", 'margin-right': '10px'}),
                        dcc.Dropdown(
                            id='plot_by4',
                            options=[{"label": i, "value": i} for i in ['growth rate', 'active dispersal rate', 'resuscitation rate', 'basal metabolic rate', 'bmr reduction in dormancy', 'immigration rate', 'resource quota', 'body size']],
                            value=None,
                            style={'width': '180px','font-size': "100%",},
                            ),
                            ],
                            style={'display': 'inline-block',
                                     'width': '170px','font-size': "100%",'margin-right': '40px'
                                      },
                        ),
                    html.Div(id="y_options",
                        children=[
                        html.B("Y-variable", style={'font-size': "100%", 'margin-right': '10px'}),
                        dcc.Dropdown(
                            id='plot_by5',
                            options=[{"label": i, "value": i} for i in ['growth rate', 'active dispersal rate', 'resuscitation rate', 'basal metabolic rate', 'bmr reduction in dormancy', 'immigration rate', 'resource quota', 'body size']],
                            value=None,
                            style={'width': '180px','font-size': "100%",},
                            ),
                            ],
                            style={'display': 'inline-block',
                                     'width': '170px','font-size': "100%",'margin-right': '40px'
                                      },
                    ),
                    html.Div(id="model_options",
                        children=[
                        html.B("Model", style={'font-size': "100%", 'margin-right': '10px'}),
                        dcc.Dropdown(
                            id='model',
                            options=[{"label": i, "value": i} for i in ['Linear', 'Quadratic', 'Cubic', 'Power law']],
                            value='Linear',
                            style={'width': '140px','font-size': "100%",},
                            ),
                            ],
                            style={'display': 'inline-block',
                                     'width': '170px','font-size': "100%",'margin-right': '40px'
                                      },
                    ),
                    html.Hr(),
                    html.Button('Plot', id='btn6', n_clicks=0,
                    style={#'width': '10%',
                            'display': 'inline-block',
                            #'margin-right': '2%',
                            #'margin-left': '-10px',
                            #'horizontal-align': 'left',
                            #'vertical-align': 'top',
                        },
                    ),
                html.Hr(),
                dcc.Graph(id='xy_fig', style={'width': '100%', 'background-color': '#f0f0f0','padding': '0px', 'margin-bottom': '0px', 'margin-right': '0px','margin-left': '0px','height': '440px', #'display': 'inline-block',
                },),],
            #style={'width': '100%', 'background-color': '#f0f0f0','padding': '0px', 'margin-bottom': '0px',
            #    'margin-right': '0px','margin-left': '0px','height': '500px', #'display': 'inline-block',
            #    },
          ),],
          style={'width': '100%', 'background-color': '#f0f0f0','padding': '0px', 'margin-bottom': '0px',
          'margin-right': '0px','margin-left': '0px','height': '600px', #'display': 'inline-block',
          },
          )
#########################################################################################
############################### IBM ANIMATION CARD ######################################
#########################################################################################

    
#########################################################################################
######################### DASH APP TABLE FUNCTIONS ######################################
#########################################################################################


    
#########################################################################################
################################# DASH APP LAYOUT #######################################
#########################################################################################


app.layout = html.Div([
    
    dcc.Store(id='main_df1', storage_type='memory'),
    dcc.Store(id='main_df2', storage_type='memory'),
    dcc.Store(id='species_1', storage_type='memory'),
    dcc.Store(id='species_2', storage_type='memory'),
    dcc.Store(id='resources_1', storage_type='memory'),
    dcc.Store(id='resources_2', storage_type='memory'),
    
    html.Div(id='placeholder1', style={'display': 'none'}),
    
    html.Div(id='N1', style={'display': 'none'}),
    html.Div(id='N2', style={'display': 'none'}),
    html.Div(id='S1', style={'display': 'none'}),
    html.Div(id='S2', style={'display': 'none'}),
    html.Div(id='R1', style={'display': 'none'}),
    html.Div(id='R2', style={'display': 'none'}),
    
    
    html.Div(
            style={'background-color': '#f9f9f9'},
            id="banner1",
            className="banner",
            children=[html.Img(src=app.get_asset_url("plotly_logo.png"),
                               style={'textAlign': 'right'})],
        ),
    
    html.Div(
            id="top-column1",
            className="ten columns",
            children=[description_card1()],
            style={'width': '95.3%',
                    'display': 'inline-block',
                    'border-radius': '15px',
                    'box-shadow': '1px 1px 1px grey',
                    'background-color': '#f0f0f0',
                    'padding': '10px',
                    'margin-bottom': '10px',
            },
        ),
    
    html.Div(id="left-column1", className="one columns",
            children=[control_card1()],
            style={'width': '24%',
                    'display': 'inline-block',
                    'border-radius': '15px',
                    'box-shadow': '1px 1px 1px grey',
                    'background-color': '#f0f0f0',
                    'padding': '10px',
                    'margin-bottom': '10px'},
        ),
    
    html.Div(id="right-column2", className="one columns",
            children=[
                html.Div(
                id="IBM_animation",
                children=[html.Div(
                            id="model_animation",
                            children=[
                                html.Div(
                                    id="plot_by_box",
                                    children=[
                                    html.B("Plot individuals by",
                                        style={'display': 'inline-block',
                                       'margin-right': '20px',
                                       'width': '100%',
                                    },),
                                    dcc.Dropdown(
                                        id='plot_by',
                                        options=[{"label": i, "value": i} for i in ['body size', 'resource quota']],
                                        value='body size',
                                        style={#'display': 'inline-block',
                                           #'vertical-align': 'top',
                                           #'margin-right': '40px',
                                           'width': '70%',
                                        },
                                        ),
                                    ],
                                    style={'display': 'inline-block',
                                            'vertical-align': 'top',
                                            'margin-right': '40px',
                                            'width': '30%',
                                         },
                                    ),
                                
                                html.Div(
                                    id="update_text_box1",
                                    children=[
                                        html.H6(id='Nc_S_R',
                                            style={
                                            #'display': 'inline-block',
                                            #'width': '50%',
                                            #'font-size': '24',
                                            }),
                                        
                                        ],
                                    style={'display': 'inline-block'},
                                    ),
                                    
                                html.Hr(),
                                dcc.Graph(id='model_animation_fig')],
                                
                            style={#'width': '100%',
                                #'display': 'inline-block',
                                'background-color': '#f0f0f0',
                                'padding': '0px',
                                'margin-bottom': '0px',
                                'margin-right': '0px',
                                'margin-left': '0px',
                                'height': '570px',
                                },
                          ),
                          dcc.Interval(
                              id='interval',
                              interval = 1000,
                              n_intervals = 0,
                              max_intervals = -1,
                              disabled = True,
                          ),
                          ],
                        ),
                ],
                style={'width': '69.3%',
                        'height': '587px',
                        'display': 'inline-block',
                        'border-radius': '15px',
                        'box-shadow': '1px 1px 1px grey',
                        'background-color': '#f0f0f0',
                        'padding': '10px',
                        'margin-bottom': '10px',
                    },
            ),
            
            
    html.Div(id="Time-series1", className="one columns",
        children=[time_series1()],
        style={'width': '47%',
                'display': 'inline-block',
                'border-radius': '15px',
                'box-shadow': '1px 1px 1px grey',
                'background-color': '#f0f0f0',
                'padding': '10px',
                'margin-bottom': '10px'},
    ),
    
    html.Div(id="distribution_1", className="one columns",
        children=[distribution_1()],
        style={'width': '47%',
                'display': 'inline-block',
                'border-radius': '15px',
                'box-shadow': '1px 1px 1px grey',
                'background-color': '#f0f0f0',
                'padding': '10px',
                'margin-bottom': '10px'},
    ),
    
    html.Div(id="xy_1", className="one columns",
        children=[xy_1()],
        style={'width': '47%',
                'display': 'inline-block',
                'border-radius': '15px',
                'box-shadow': '1px 1px 1px grey',
                'background-color': '#f0f0f0',
                'padding': '10px',
                'margin-bottom': '10px'},
    ),
    
])


#########################################################################################
############################    Callbacks   #############################################
#########################################################################################

@app.callback([Output('placeholder1', 'children'),
               Output('interval', 'disabled'),
               Output('main_df1', 'clear_data'),
               Output('species_1', 'clear_data'),
               Output('resources_1', 'clear_data'),
               Output('main_df2', 'clear_data'),
               Output('species_2', 'clear_data'),
               Output('resources_2', 'clear_data'),
               Output('btn2', 'n_clicks'),
               Output('btn3', 'n_clicks'),
               Output('btn4', 'n_clicks'),
               Output('interval', 'n_intervals'),
               ],
              [Input('btn1', 'n_clicks')],
              prevent_initial_call=True,
    )
def update_df(n_clicks1):
    return 1, False, True, True, True, True, True, True, 0, 0, 0, 0
    
    
    
@app.callback([Output('main_df1', 'data'),
               Output('species_1', 'data'),
               Output('resources_1', 'data'),
               Output('N1', 'children'),
               Output('S1', 'children'),
               Output('R1', 'children'),
               Output('btn-rarefy', 'n_clicks'),],
              [Input('main_df2', 'data'),
               Input('species_2', 'data'),
               Input('resources_2', 'data'),
               Input('N2', 'children'),
               Input('S2', 'children'),
               Input('R2', 'children')],
    )
def update_df(df, species, resources, N, S, R):
    return df, species, resources, N, S, R, 0
    

@app.callback([Output('model_animation_fig', 'figure'),
               Output('main_df2', 'data'),
               Output('species_2', 'data'),
               Output('resources_2', 'data'),
               Output('interval', 'interval'),
               Output('Nc_S_R', 'children'),
               Output('N2', 'children'),
               Output('S2', 'children'),
               Output('R2', 'children'),
               ],
              [Input('interval', 'n_intervals')],
              [State('interval', 'max_intervals'),
               State('interval', 'disabled'),
               State('main_df1', 'data'),
               State('species_1', 'data'),
               State('resources_1', 'data'),
               State('S', 'value'),
               State('Q', 'value'),
               State('R', 'value'),
               State('btn2', 'n_clicks'),
               State('btn3', 'n_clicks'),
               State('plot_by', 'value'),
               State('immigration', 'value'),
               State('immigration_on_off', 'value'),
               State('reproduction_on_off', 'value'),
               State('death_on_off', 'value'),
               State('active_dispersal_on_off', 'value'),
               State('N1', 'children'),
               State('S1', 'children'),
               State('R1', 'children'),
               State('btn-rarefy', 'n_clicks'),
              ],
            )
def run_model(n_intervals, max_intervals, disabled, individuals, species, resources, S, Q, R0, n_clicks2, n_clicks3, plot_by, immigration_rate, imm_toggle, repr_toggle, death_toggle, act_disp_toggle, N1, S1, R1, n_clicks4):
    
    print('rarefy clicks:', n_clicks4)
    w = 100
    h = 50
    N = 1000
    
    if Q is None or math.isnan(Q) == True:
        Q = 1
    if R0 is None or math.isnan(R0) == True:
        R0 = 0
    
    if resources is None:
        resources = pd.DataFrame(columns=['Resource ID', 'x_coord', 'y_coord', 'size'])

    else:
        resources = pd.read_json(resources)
        
    R = 0
    if resources.shape[0] > 0:
        R = np.sum(resources['size'])
        
    fig_data ={}
    figure = go.Figure(
            data = fig_data,
            layout = go.Layout(
                
                xaxis = dict(
                    title = dict(
                        text = None,
                        font = dict(
                            family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                            " Helvetica, Arial, sans-serif",
                            size = 18,
                        ),
                    ),
                    #rangemode="tozero",
                    #zeroline=True,
                    visible=False,
                    showticklabels = False,
                ),
                                
                yaxis = dict(
                    title = dict(
                        text = None,
                        font = dict(
                            family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                            " Helvetica, Arial, sans-serif",
                            size = 18,
                        ),
                    ),
                    #rangemode="tozero",
                    #zeroline=True,
                    visible=False,
                    showticklabels = False,
                ),
                                
                margin = dict(l=0, r=0, b=0, t=0),
                showlegend = False,
                height = 500,
                paper_bgcolor = "rgb(245, 247, 249)",
                plot_bgcolor = "rgb(245, 247, 249)",
            ),
        )
    
    figure.update_xaxes(range=[0, w])
    figure.update_yaxes(range=[0, h])
    
    if disabled is True or n_clicks3 & 1 == True:
        Nc_S_R = 'N = 0' + ' | ' + 'S = 0' + ' | ' + 'Total resources = 0'
        return figure, None, None, None, 1000, Nc_S_R, [0], [0], [0]
    
    if n_clicks2 & 1 == True:
        raise PreventUpdate

    if species is None:
        # declare initial dataframe
        species = pd.DataFrame(columns=['Species ID'])
        
        # assign species IDs and traits
        species['Species ID'] = np.random.randint(0, 0xFFFFFF, size=S)
        species['growth rate'] = np.random.uniform(0.001, 1, size=S)
        species['active dispersal rate'] = np.random.uniform(0, 20, size=S)
        species['resuscitation rate'] = np.random.uniform(0.001, 1, size=S)
        species['basal metabolic rate'] = np.random.uniform(0.001, 1, size=S)
        species['bmr reduction in dormancy'] = np.random.uniform(0.001, 1, size=S)
        species['immigration rate'] = np.random.uniform(0.001, 1, size=S) #1 - species['active dispersal rate']/np.sum(species['active dispersal rate'])
        
        for i in list(range(1, 2)):
            species['resource efficiency ' + str(i)] = np.random.uniform(0.001, 1, size=S)
        
        sp_ids = species['Species ID'].tolist()
        clrs = []
        for id in sp_ids:
            clr = "#" + "%06x" % id
            clrs.append(clr)
        species['color'] = clrs
        
    else:
        species = pd.read_json(species)
            
    if individuals is None and n_intervals < 2:
        individuals = species.copy(deep=True)
        individuals['Ind ID'] = list(range(S))
        individuals['age'] = [0] * individuals.shape[0]
        individuals['x_coord'] = 0
        individuals['y_coord'] = np.random.uniform(0, h, size=S)
        individuals['resource quota'] = 10 #np.random.uniform(0, 100, size=S)
        individuals['body size'] = [10]*S
        individuals['metabolic state'] = [1] * S # 0 = dormant, 1 = active
            
    elif individuals is None and n_intervals >= 1:
        individuals = species.copy(deep=True)
        individuals = individuals[individuals['growth rate'] < 0]
        
    else:
        individuals = pd.read_json(individuals)
        if individuals.shape[0] == 0:
            raise PreventUpdate
        elif individuals.shape[0] > 1000 and n_clicks4 > 0:
            individuals = individuals.sample(n=1000, replace=False)
        
        
            
    ####################################################
    ########### SIMULATE RESOURCE INFLOW ###############
    ####################################################
    
    r2 = pd.DataFrame(columns=['Resource ID', 'x_coord', 'y_coord', 'size'])
    r2['Resource ID'] = [1]
    r2['x_coord'] = 0
    r2['y_coord'] = np.random.uniform(0, h, size=1)
    r2['size'] = [R0*Q]
    resources = pd.concat([resources, r2], ignore_index=True)
    
    ####################################################
    ############## SIMULATE IMMIGRATION ################
    ####################################################
    
    if immigration_rate > 0 and imm_toggle == ' on':
        im = int(immigration_rate*Q)
        
        maxID = 0
        if individuals.shape[0] > 0:
            maxID = 1 + np.max(individuals['Ind ID'])
            
        if im > 0:
            i2 = species.sample(n=im, replace=True, weights=species['immigration rate'])
            i2['Ind ID'] = list(range(im))
            i2['Ind ID'] = i2['Ind ID'] + maxID
            i2['age'] = [0] * i2.shape[0]
            i2['x_coord'] = 0
            i2['y_coord'] = np.random.uniform(0, h, size=im)
            i2['resource quota'] = 10
            i2['body size'] = [10]*im
            i2['metabolic state'] = [1] * im
            individuals = pd.concat([individuals, i2], ignore_index=True)
    
    ####################################################
    ### SIMULATE PROCESSES AMONG ACTIVE INDIVIDUALS ####
    ####################################################
    
    df = []
    if individuals.shape[0] > 0:
        individuals = individuals[individuals['resource quota'] >= 0]
        individuals['x_coord'] = individuals['x_coord'] + (Q*0.01)*w
        
        df_a = individuals[individuals['metabolic state'] == 1]
        
        n = df_a.shape[0]
        if n > 0:
            
            # resource consumption
            if resources.shape[0] > 0:
                res = resources['size']
                res_weights = res/np.sum(res)
                
                R = np.sum(res)
                r_per_capita = R/n
                D = R#/(w)
                p = np.random.binomial(1, D/(1 + D), size=n)
                total_consumed = np.minimum(r_per_capita, df_a['resource efficiency 1'] * df_a['body size']) * p
                df_a['resource quota'] = df_a['resource quota'] + total_consumed
                
                R -= np.sum(total_consumed)
                if R < 0:
                    R = 0
                    
                resources['size'] = R*res_weights
            
            # growth
            g = np.minimum(df_a['body size'] * df_a['growth rate'], df_a['resource quota'])
            df_a['body size'] = df_a['body size'] + g #* df_a['metabolic state']
            df_a['resource quota'] = df_a['resource quota'] - g #* df_a['metabolic state']
            
            # active dispersal inside the system
            if act_disp_toggle == ' on':
                d = np.minimum(df_a['x_coord'], df_a['active dispersal rate'])
                d = np.minimum(d, df_a['resource quota'])
                df_a['x_coord'] = df_a['x_coord'] - d
                df_a['resource quota'] = df_a['resource quota'] - d/w
            
            # active maintenance
            df_a['resource quota'] = df_a['resource quota'] - df_a['basal metabolic rate']
            
            # increase age
            df_a['age'] = df_a['age'] + 1
            
            # death
            if death_toggle == ' on':
                df_a = df_a[df_a['resource quota'] >= 0]
            
            # outflow of active individuals
            if df_a.shape[0] > 0:
                df_a = df_a[df_a['x_coord'] <= w]
            

            n = df_a.shape[0]
            if n > 0:
                df_a['resource quota'] = df_a['resource quota'].apply(lambda x : x if x > 0 else 0)
                df_a['resource quota'].fillna(0, inplace=True)
                
                # reproduction
                if repr_toggle == ' on':
                    ri = df_a['resource quota']/df_a['basal metabolic rate']
                    p = ri/(1 + ri) * df_a['body size']/(20 + df_a['body size']) * df_a['age']/(20 + df_a['age'])
                    df_a['reproduce'] = p
                    df_a['reproduce'].replace(np.inf, 0, inplace=True)
                    df_a['reproduce'].replace(-np.Inf, 0, inplace=True)
                    df_a['reproduce'].fillna(0, inplace=True)
                    
                    x = df_a['reproduce'].tolist()
                    df_a['reproduce'] = np.random.binomial(1, df_a['reproduce'], size = n)
                    
                    reproduce_no = df_a[df_a['reproduce'] == 0]
                    reproduce_no.drop(labels=['reproduce'], axis=1, inplace=True)
                    
                    reproduce_yes = df_a[df_a['reproduce'] == 1]
                    reproduce_yes.drop(labels=['reproduce'], axis=1, inplace=True)
                    
                    if reproduce_yes.shape[0] > 0:
                        reproduce_yes['body size'] = reproduce_yes['body size']/2
                        reproduce_yes['resource quota'] = reproduce_yes['resource quota']/2
                        progeny = reproduce_yes.copy(deep=True)
                        progeny['age'] = 0
                        progeny['Ind ID'] = progeny['Ind ID'] + np.max(individuals['Ind ID'])
                        
                        progeny['y_coord'] = progeny['y_coord'] + np.random.uniform(-1, 1, size=progeny.shape[0])
                        progeny['y_coord'].values[progeny['y_coord'].values < 0] = 0
                        progeny['y_coord'].values[progeny['y_coord'].values > h] = h
                        
                        # merge dataframes of active individuals
                        df_a = pd.concat([reproduce_no, reproduce_yes, progeny], ignore_index=True)
                    else:
                        df_a = pd.concat([reproduce_no, reproduce_yes], ignore_index=True)
                
                # transition to dormancy
                lambda_ = df_a['resource quota']/df_a['basal metabolic rate']
                p = 1/(1+lambda_) * df_a['age']/(10+df_a['age'])
                df_a['metabolic state'] = 1 - np.random.binomial(1, p, size = df_a.shape[0])
                df_a['symbol'] = ['circle']*df_a.shape[0]
            
        ####################################################
        #### SIMULATE PROCESS AMONG DORMANT INDIVIDUALS ####
        ####################################################
        
        df_d = individuals[individuals['metabolic state'] == 0]
        n = df_d.shape[0]
        if n > 0:
        
            # dormant maintenance
            df_d['resource quota'] = df_d['resource quota'] - df_d['basal metabolic rate'] * df_d['bmr reduction in dormancy']
            
            df_d['resource quota'] = df_d['resource quota'].apply(lambda x : x if x > 0 else 0)
            df_d['resource quota'].fillna(0, inplace=True)
            
            # death
            if death_toggle == ' on':
                df_d = df_d[df_d['resource quota'] >= 0]
                df_d['symbol'] = ['circle-open']*df_d.shape[0]
            
            if df_d.shape[0] > 0:
                # transition to activity
                df_d['metabolic state'] = np.random.binomial(1, df_d['resuscitation rate'], size = df_d.shape[0])
                
                # increase age
                df_d['age'] = df_d['age'] + 1
                
                # outflow of individuals
                df_d = df_d[df_d['x_coord'] <= w]
        
    
        if df_a.shape[0] > 0 and df_d.shape[0] > 0:
            df = pd.concat([df_a, df_d], ignore_index=True)
        elif df_a.shape[0] > 0:
            df = df_a.copy(deep=True)
        elif df_d.shape[0] > 0:
            df = df_d.copy(deep=True)
        else:
            df = None
          
    if isinstance(df, list) == True:
        df = None

    ########### SIMULATE RESOURCE FLOW ##############
    if resources.shape[0] > 0:
        resources['x_coord'] = resources['x_coord'] + (Q*0.01)*w
        resources = resources[resources['x_coord'] <= w]
        #print(n_intervals, ' resources:', np.sum(resources['size']), '\n')
    
    ####################################################
    ############### CHECK DATAFRAMES ###################
    ####################################################
    
    if resources.shape[0] == 0:
        resources = None
        
    if df is None and resources is None:
        Nc_S_R = 'N = 0' + ' | S = 0' + ' | Total resources = 0'
        N1.append(0)
        S1.append(0)
        R1.append(0)
        return figure, df, species.to_json(), resources, 1000, Nc_S_R, N1, S1, R1
        
    elif df is None:
        R = np.sum(resources['size'])
        Nc_S_R = 'N = 0' + ' | ' + 'S = 0' + ' | ' + 'Total resources = ' + str(np.round(R,3))
        N1.append(0)
        S1.append(0)
        R1.append(R)
        return figure, df, species.to_json(), resources.to_json(), 1000, Nc_S_R, N1, S1, R1
    
    ####################################################
    ################ GENERATE FIGURE ###################
    ####################################################
    
    if df.shape[0] > 0:
        if df.shape[0] > 1000 and n_clicks4 > 0:
            df = df.sample(n=1000, replace=False)
            
        fig_data = []
        
        fig_data.append(go.Scatter(
                            x = df['x_coord'],
                            y = df['y_coord'],
                            text = 'Body size: ' + np.round(df['body size'], 3).astype(str) + '<br>' + 'Resource quota: ' + np.round(df['resource quota'], 3).astype(str) + '<br>' + 'BMR: ' + np.round(df['basal metabolic rate'], 3).astype(str) + '<br>' + 'BMR reduction in dormancy: ' + np.round(df['bmr reduction in dormancy'], 3).astype(str) + '<br>' + 'Resource use efficiency: ' + np.round(df['resource efficiency 1'], 3).astype(str) + '<br>' + 'Resuscitation rate: ' + np.round(df['resuscitation rate'], 3).astype(str) + '<br>' + 'Active dispersal rate: ' + np.round(df['active dispersal rate'], 3).astype(str) + '<br>' + 'Growth rate: ' + np.round(df['growth rate'], 3).astype(str),
                            mode = "markers",
                            marker_size= 4 + df[plot_by]**0.75,
                            marker_color=df['color'],
                            marker_symbol=df['symbol'],
                            
                        )
                    )
                        
        figure = go.Figure(
                data = fig_data,
                layout = go.Layout(
                    xaxis = dict(
                        title = dict(
                            text = None,
                            font = dict(
                                family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                                " Helvetica, Arial, sans-serif",
                                size = 18,
                            ),
                        ),
                        rangemode="tozero",
                        zeroline=True,
                        #visible=False,
                        showticklabels = False,
                    ),
                                    
                    yaxis = dict(
                        title = dict(
                            text = None,
                            font = dict(
                                family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                                " Helvetica, Arial, sans-serif",
                                size = 18,
                            ),
                        ),
                        rangemode="tozero",
                        zeroline=True,
                        #visible=False,
                        showticklabels = False,
                    ),
                                    
                    margin = dict(l=0, r=0, b=0, t=0),
                    showlegend = False,
                    height = 500,
                    paper_bgcolor = "rgb(245, 247, 249)",
                    plot_bgcolor = "rgb(245, 247, 249)",
                ),
            )
        
        figure.update_xaxes(range=[0, w])
        #figure.update_xaxes(visible=False)
        figure.update_yaxes(range=[0, h])
        #figure.update_yaxes(visible=False)
        
        Nc = str(df.shape[0])
        S = str(len(list(set(df['Species ID'].tolist()))))
        R = str(np.round(np.sum(resources['size']), 3))
        
    interval = max([500, df.shape[0]**0.95])
    Nc_S_R = 'N = ' + Nc + ' | ' + 'S = ' + S + ' | ' + 'Total resources = ' + R #+ ' | Time step = ' + str(n_intervals)
    
    N1.append(float(Nc))
    S1.append(float(S))
    R1.append(float(R))
    
    if resources is None:
        return figure, df.to_json(), species.to_json(), resources, interval, Nc_S_R, N1, S1, R1
        
    return figure, df.to_json(), species.to_json(), resources.to_json(), interval, Nc_S_R, N1, S1, R1
    
    
@app.callback(Output('time_series_fig', 'figure'),
             [Input('btn4', 'n_clicks')],
             [State('plot_by2', 'value'),
              State('N1', 'children'),
              State('S1', 'children'),
              State('R1', 'children')],
              )
def time_series_plot(n_clicks, var_lab, N, S, R):
    x = []
    
    if var_lab == 'Total abundance (N)':
        if N is None:
            x = []
        else:
            x = list(N)
    elif var_lab == 'Species richness (S)':
        if S is None:
            x = []
        else:
            x = list(S)
    elif var_lab == 'Total resources':
        if R is None:
            x = []
        else:
            x = list(R)

    fig_data = []
    fig_data.append(go.Scatter(
                        x = list(range(len(x))),
                        y = x,
                        mode='lines+markers',
                        marker_size= 10,
                        marker_color='#99ccff',
                        marker_symbol='circle-open',
                        
                    )
                )
                    
    figure = go.Figure(
            data = fig_data,
            layout = go.Layout(
                xaxis = dict(
                    title = dict(
                        text = 'Time',
                        font = dict(
                            family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                            " Helvetica, Arial, sans-serif",
                            size = 18,
                        ),
                    ),
                    rangemode="tozero",
                    zeroline=True,
                    visible=True,
                    showticklabels = True,
                ),
                                
                yaxis = dict(
                    title = dict(
                        text = var_lab,
                        font = dict(
                            family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                            " Helvetica, Arial, sans-serif",
                            size = 18,
                        ),
                    ),
                    rangemode="tozero",
                    zeroline=True,
                    visible=True,
                    showticklabels = True,
                ),
                                
                margin = dict(l=0, r=0, b=0, t=0),
                showlegend = False,
                height = 400,
                paper_bgcolor = "rgb(245, 247, 249)",
                plot_bgcolor = "rgb(245, 247, 249)",
            ),
        )
    
    return figure
    


@app.callback(Output('distribution_fig', 'figure'),
            [Input('btn5', 'n_clicks')],
            [State('plot_by3', 'value'),
             State('main_df1', 'data'),
            ],
            )
def distribution_plot(n_clicks, var_lab, main_df):
        x = []
        
        if main_df is None:
            x, y = [0]*100, [0]*100
        else:
            main_df = pd.read_json(main_df)
            x = main_df[var_lab].dropna()
            x, y = get_kdens_choose_kernel(x, 0.5)
            
        fig_data = []
        fig_data.append(go.Scatter(
                                x = x,
                                y = y,
                                mode='lines',
                                #marker_size= 10,
                                #marker_color='#99ccff',
                                #marker_symbol='circle-open',
                            )
                        )
                            
        figure = go.Figure(
                    data = fig_data,
                    layout = go.Layout(
                        xaxis = dict(
                            title = dict(
                                text = var_lab,
                                font = dict(
                                    family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                                    " Helvetica, Arial, sans-serif",
                                    size = 18,
                                ),
                            ),
                            #rangemode="tozero",
                            #zeroline=True,
                            visible=True,
                            showticklabels = True,
                        ),
                                        
                        yaxis = dict(
                            title = dict(
                                text = 'kernel density',
                                font = dict(
                                    family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                                    " Helvetica, Arial, sans-serif",
                                    size = 18,
                                ),
                            ),
                            #rangemode="tozero",
                            #zeroline=True,
                            visible=True,
                            showticklabels = True,
                        ),
                                        
                        margin = dict(l=0, r=0, b=0, t=0),
                        showlegend = False,
                        height = 400,
                        paper_bgcolor = "rgb(245, 247, 249)",
                        plot_bgcolor = "rgb(245, 247, 249)",
                    ),
                )
            
        return figure



@app.callback(Output('xy_fig', 'figure'),
            [Input('btn6', 'n_clicks')],
            [State('plot_by4', 'value'),
             State('plot_by5', 'value'),
             State('main_df1', 'data'),
            ],
            )
def distribution_plot(n_clicks, x_var, y_var, main_df):
        x = []
        
        if main_df is None:
            x, y = [0]*100, [0]*100
        else:
            main_df = pd.read_json(main_df)
            tdf = main_df.filter(items=[x_var, y_var], axis=1)
            del main_df
            
            tdf.dropna(inplace=True)
            x = tdf[x_var]
            y = tdf[y_var]
            
            
        fig_data = []
        fig_data.append(go.Scatter(
                                x = x,
                                y = y,
                                mode='markers',
                                marker_size= 10,
                                marker_color='#99ccff',
                                marker_symbol='circle',
                            )
                        )
                            
        figure = go.Figure(
                    data = fig_data,
                    layout = go.Layout(
                        xaxis = dict(
                            title = dict(
                                text = x_var,
                                font = dict(
                                    family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                                    " Helvetica, Arial, sans-serif",
                                    size = 18,
                                ),
                            ),
                            #rangemode="tozero",
                            #zeroline=True,
                            visible=True,
                            showticklabels = True,
                        ),
                                        
                        yaxis = dict(
                            title = dict(
                                text = y_var,
                                font = dict(
                                    family = '"Open Sans", "HelveticaNeue", "Helvetica Neue",'
                                    " Helvetica, Arial, sans-serif",
                                    size = 18,
                                ),
                            ),
                            #rangemode="tozero",
                            #zeroline=True,
                            visible=True,
                            showticklabels = True,
                        ),
                                        
                        margin = dict(l=0, r=0, b=0, t=0),
                        showlegend = False,
                        height = 440,
                        paper_bgcolor = "rgb(245, 247, 249)",
                        plot_bgcolor = "rgb(245, 247, 249)",
                    ),
                )
            
        return figure
        
#########################################################################################
############################# Run the server ############################################
#########################################################################################

if __name__ == "__main__":
    app.run_server(host='0.0.0.0', debug = False) # modified to run on linux server
