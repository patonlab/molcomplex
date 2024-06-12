import pandas as pd
import sys
import dash
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
import dash_dangerously_set_inner_html as dhtml
from dash.dependencies import Input, Output, State

from settings import config
from complex_funcs import *

# Code adapted from https://towardsdatascience.com/web-development-with-python-dash-complete-tutorial-6716186e09b3 
# and https://iwatobipen.wordpress.com/2019/02/16/make-interactive-dashboard-with-dash2-chemoinformatcs-rdkit/
#by Guilian Luchini

# mcdir = '../molcomplex/'
# sys.path.insert(1, mcdir)

import plotly.graph_objs as go

from molcomplex.complex_object import mol_complex
a = mol_complex(['C'], twc=False, linked=True)
deslist = a.columns.values[2:]

#initialize app
app = dash.Dash(
    name=config.app_name,
    assets_folder="static",
     external_stylesheets=[
        dbc.themes.LUX,
        config.fontawesome
    ]
)
app.title = config.app_name


########################## Navbar ##########################
"""
Displays About and Links pages
-Possibly add wiki pages later for various complexity scores
"""

# Output
navbar = dbc.Nav(className="nav nav-pills", children=[
    
    ## about
    dbc.NavItem(html.Div([
        dbc.NavLink("About", href="/", id="about-popover", active=False),
        dbc.Popover(id="about", is_open=False, target="about-popover", children=[
            dbc.PopoverHeader("How it works"), dbc.PopoverBody(config.about)
        ])
    ])),
    ## links
    dbc.DropdownMenu(label="Links", nav=True, children=[
        dbc.DropdownMenuItem([html.I(className="fa fa-window-maximize"), "  C-CAS"], href=config.contacts, target="_blank"), 
        dbc.DropdownMenuItem([html.I(className="fa fa-github"), "  Code"], href=config.code, target="_blank")
    ])
])


# Callbacks
@app.callback(output=[Output(component_id="about", component_property="is_open"), 
                      Output(component_id="about-popover", component_property="active")], 
              inputs=[Input(component_id="about-popover", component_property="n_clicks")], 
              state=[State("about","is_open"), State("about-popover","active")])
def about_popover(n, is_open, active):
    if n:
        return not is_open, active
    return is_open, active



########################## Body ##########################
# Input
inputs = dbc.Col([
    dbc.Label('Input SMILES',html_for='input-SMILES'),
    dcc.Input(
            id="input-SMILES",
            type='text',
            placeholder="Enter SMILES ( For Example: CC(=O)OC1=CC=CC=C1C(=O)O ) ",
            style={
                    'width': '100%',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'textAlign': 'left',
            },
            debounce = True
        ),
    dbc.Row([
        dbc.Spinner([
            html.Div([html.Div(id="inputmolimg")],
                style={
                        'width': '310px',
                        'height': '260px',
                        'display': 'flex',
                        'justify-content': 'center',
                        'align-items': 'center',
                        'lineHeight': '60px',
                        'borderWidth': '5px',
                        'borderStyle': 'inset',
                        'borderRadius': '1px',
                        #'borderColor': '#4ccce6',
                        'textAlign': 'center',
                        'margin-top': '5px',
                        'margin-bottom': '5px',
                        'margin-left': 'auto',
                        'margin-right': 'auto'
                })],
            color='primary',type='grow'
            )
        ]),

    html.Br(),

    dbc.Label("Complexity Score", html_for="score"), 
    dcc.Dropdown(id='score',
            value='BALABAN',
            options=[{'label': key, 'value': key} for key in deslist],
            style={'width':'100%', 'textAlign':'center'},
            ),

    dbc.Label("Sort By", html_for="sortby"),    
    dcc.Dropdown(id='sortby',
            value='Increasing complexity',
            options=[{'label': key, 'value': key} for key in ['Increasing complexity','Decreasing complexity']],
            style={'width':'100%', 'textAlign':'center'},
            ),

    dbc.Label("Number of Bonds to Break", html_for="num-bonds"), 
    dbc.Input(id="num-bonds", placeholder="n bonds broken", 
              type="number", value="1",min=1, max=10,step=1),
    
    dbc.Label("Number of Molecules Per Row", html_for="mols-per-row"), 
    dbc.Input(id="mols-per-row", placeholder="n mols per row", 
              type="number", value="3",min=1,step=1),

    html.Br(),html.Br(),
    dbc.Col(dbc.Button("run", id="run", color="primary",style={ 'width': '100%',}))
])

# Output
body = dbc.Row([
        ## input
        dbc.Col(md=3, children=[
            inputs, 
            html.Br(),html.Br(),html.Br(),],
            style={'overflow':'auto'},
        ),
        ## output
        dbc.Col(md=9, children=[
            dbc.Spinner([
                ### title
                # html.H2('Enter a smiles and press "run" to begin',id="title"),
                ### download

                dbc.Badge(html.A('Download Data', id='download-excel', download="tables.xlsx", 
                    href="", 
                    target="_blank"), color="primary", pill=True,),

        
                ### plot
                #Mol Grid Image
                html.Div([html.Div(id="molgridimg")],
                    style={
                            'display': 'flex',
                            'justify-content': 'center',
                            'width': '100%',
                            'height': '600px',
                            'textAlign': 'center',
                            'margin-top': '5px',
                            'margin-bottom': '5px',
                            'margin-left': '5px',
                            'margin-right': '5px'
                    },
                    ),
                        

            ], color="primary", type="grow"), 
        ],style={'overflow':'auto'},)
])

#Callbacks - what happens upon interaction with any widgets

#single mol img
@app.callback(
    Output('inputmolimg', 'children'),
    [Input('input-SMILES', 'value'),
    Input('score', 'value'),
    ]
)
def update_img(smiles,feature):
    if smiles == None or smiles == '':
        return None
    else:
        try:
            smiles = canonicalize_smiles(smiles)
        except:
            #print error message - invalid smiles
            return None
        
    df = mol_complex([smiles],  twc=False, linked=True)

    feature_value = df[feature].values[0]
    #try:
    svg = smi2svg(smiles,feature_value,feature)
    # except:
    #     svg = 'Input a SMILES string to display molecule. (Example: CCCC)'
    return dhtml.DangerouslySetInnerHTML(svg)

@app.callback(
    output=[
        # Output(component_id='title',component_property='children'),
        Output(component_id='molgridimg', component_property='children'),
        Output(component_id='download-excel',component_property='href')
    ],
    inputs=[
        Input(component_id='run', component_property='n_clicks')
    ],
    state=[
        State('input-SMILES','value'),
        State('score', 'value'),
        State('num-bonds','value'),
        State('sortby', 'value'),
        State('mols-per-row','value')
     ])
def results(n_clicks, smiles, score, nbonds, sortby, mols_per_row):
    if smiles == None or smiles == '':
        return None,None
    else:
        try:
            smiles = canonicalize_smiles(smiles)
        except:
            #print error message - invalid smiles
            return None,None
    
    

    mols, smi_list, highlight_bonds  = parse_contents(smiles,int(nbonds))
   
    df = mol_complex(smi_list, twc=False, linked=True)
    print(df)
    df['highlight_bonds'] = highlight_bonds

    if sortby == 'Decreasing complexity':
        sorted = df.sort_values(by=score,ascending=False)
    else:
        sorted = df.sort_values(by=score,ascending=True)

    feature_list = sorted[score].values
    bonds = [list(val) for val in sorted['highlight_bonds'].values]

    svg = smis2svg(smiles,bonds,feature_list,score,sortby,int(mols_per_row))
    #transparent background
    svg = svg.replace("rect style='opacity:1.0","rect style='opacity:0.0")

    sorted = sorted.drop('highlight_bonds',axis=1)
    input_df = mol_complex([smiles], twc=False, linked=True)
    output = pd.concat([input_df,sorted])
    

    return dhtml.DangerouslySetInnerHTML(svg), download_file(output)
    
colors = {
    'background': '#111111',
    'text': '#7FDBFF'
}

app.layout = dbc.Container(
    fluid=True,
    style = {
        'background-image':'url(/static/bg-wave-white.png)',
        'background-size': 'cover',
        'background-repeat': 'no-repeat',
        },
    children=[
    dbc.NavItem(html.Img(src=app.get_asset_url("logo-removebg-preview.png"), height="200px")),
    navbar,
    html.Br(),
    body,
    
])

# link for white wave https://unsplash.com/photos/GJKx5lhwU3M?utm_source=unsplash&utm_medium=referral&utm_content=creditShareLink

if __name__=='__main__':
    app.run_server(debug=True)
