import streamlit as st
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import os
import sys
baseName = os.path.basename(__file__)
dirName = os.path.dirname(__file__)
sys.path.append(dirName + r'./')
import scalepyBack

# default empty figure
defaultFig = go.Figure()
defaultFig.update_xaxes(
    title_text = 'Period (sec)',
    range=[0,4],
    tickvals=np.arange(0,4.5,0.5),
    dtick = 1,
    showgrid = True,
    zeroline=True,
    zerolinewidth=1
)

defaultFig.update_yaxes(
    title_text = 'pSa (g)',
    range=[0,3],
    showgrid = True,
    zeroline=True,
    zerolinewidth=1
)

defaultFig.update_layout(showlegend=False, template=None, plot_bgcolor = "#F0F2F6", width=1100,height=600, title_text='No Data', title_x=0.5, legend=dict(
    yanchor="top",
    x = 1,
    xanchor="right"
    ))

with st.sidebar:
    st.markdown("# ScalePy")
    st.markdown("💻 Click [here](https://github.com/dogukankaratas/scalepy) for source code.")
    st.markdown("🧑‍💻 Use [here](https://www.linkedin.com/in/dogukankaratas/) to reach me.")


inputCol, graphCol = st.columns([1, 2.7])

with inputCol:
    st.title("AFAD Values")
    with st.form("locationForm"):
        lat = st.number_input("Latitude", 34.25, 42.95, 36.0, 0.5)
        lon = st.number_input("Longitude", 24.55, 45.95, 42.0, 0.5)
        intensity = st.selectbox("Intensity Level", ["DD1", "DD2", "DD3", "DD4"], 1)
        soil = st.selectbox('Soil Type', ('ZA', 'ZB', 'ZC', 'ZD', 'ZE'), 2)
        responseCreate = st.form_submit_button("Create Response Spectrum")

if responseCreate:
    values = scalepyBack.parameterCreator(lat, lon, intensity)
    ss = values['Ss']
    s1 = values['S1']
    defaultTarget = scalepyBack.targetSpectrum(ss, s1, soil)
    defaultTarget = defaultTarget.rename(columns={'T': 'Period (sec)', 'Sa': 'pSa (g)'})
    defaultFig = go.Figure()
    defaultFig.add_trace(go.Scatter(x = defaultTarget['Period (sec)'],
                                        y=defaultTarget['pSa (g)'],
                                        name='Response Spectrum', line=dict(color='red')))

    defaultFig.update_xaxes(
            title_text = 'Period (sec)',
            range=[0,4],
            tickvals=np.arange(0,4.5,0.5),
            dtick = 1,
            showgrid = True,
            zeroline=True,
            zerolinewidth=1
        )

    defaultFig.update_yaxes(
            title_text = 'pSa (g)',
            range=[0,3],
            showgrid = True,
            zeroline=True,
            zerolinewidth=1
        )

    defaultFig.update_layout(showlegend=True, template=None, plot_bgcolor = "#F0F2F6", width=1100,height=600, 
                                title_text='Response Spectrum', title_x=0.5, legend=dict(
                                                                yanchor="top",
                                                                x = 1,
                                                                xanchor="right")
                                )

    writeCol, tableCol = st.columns([1,3])
    with writeCol:
        st.write(values)

    with tableCol:
        st.table(defaultTarget)
    
with graphCol:
    st.plotly_chart(defaultFig)


