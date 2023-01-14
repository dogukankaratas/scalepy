import streamlit as st
import scalepyBack
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import base64

# default empty figure
defaultFig = go.Figure()
defaultFig.update_xaxes(
    title_text = 'Period (sec)',
    range=[0,4],
    tickvals=np.arange(0,4.5,0.5),
    dtick = 1
)

defaultFig.update_yaxes(
    title_text = 'pSa (g)',
    range=[0,3]
)

defaultFig.update_layout(showlegend=False, width=1000,height=600, title="No Data", legend=dict(
    yanchor="top",
    x = 1,
    xanchor="right",
    ))

st.set_page_config(page_title="Scalepy GUI",layout='wide')

st.markdown("""
    <style>
      section[data-testid="stSidebar"][aria-expanded="true"]{
        width: 60% !important;
      }
      section[data-testid="stSidebar"][aria-expanded="false"]{
        width: 60% !important;
      }
    </style>""", unsafe_allow_html=True)

with st.sidebar:
    st.markdown("## ScalePy - Ground Motion Selection and Scaling Framework")
    responseTab, selectionTab, scalingTab = st.tabs(["Define Response Spectrum", "Filter Ground Motions", "Scale Ground Motions"])
    with responseTab:
        tbdyTab, asceTab, userDefinedTab = st.tabs(["acc. to TBDY-2018", "acc. to ASCE7-22", "User Defined"])
        with tbdyTab:
            with st.form("tbdyForm"):
                Ss = st.number_input('Spectral Acceleration at Short Periods (Ss)', value=0.8)
                S1 = st.number_input('Spectral Acceleration at 1 sec (S1)', value=0.4)
                soil = st.selectbox('Soil Type', ('ZA', 'ZB', 'ZC', 'ZD', 'ZE'), 3)
                responseButton = st.form_submit_button('Create Response Spectrum')
            st.markdown("#### Click [here](https://github.com/dogukankaratas/scalepy) for the source code :rocket:")
        with asceTab:
            st.markdown("## Still under work. Maybe you can contribute this section?")
        with userDefinedTab:
            st.markdown("## Still under work. Maybe you can contribute this section?")
        
    with selectionTab:
        with st.form("selectionForm"):
            period = st.number_input("Structure Period", 0.0, 10.0, 1.0, 0.1)
            magnitudeRange = st.slider('Magnitude Range', 0.0, 12.0, (4.0, 9.0), step=0.2)
            vs30Range = st.slider('Vs30 Range', 0, 1500, (180, 360), step=10)
            rjbRange = st.slider('RJB Range', 0, 3000, (0, 500), step=10)
            faultMechanism = st.selectbox('Fault Mechanism', ["Strike - Slip", "Normal", "Reverse", "Oblique", "Reverse - Oblique", "Normal - Oblique"])
            duration575Range = st.slider('%5-%75 Duration Range', 0, 100, (0, 50), step=5)
            duration595Range = st.slider('%5-%95 Duration Range', 0, 100, (0, 50), step=5)
            ariasIntensity = st.slider('Arias Intensity Range', 0, 10, (0, 5), step=1)
            filterButton = st.form_submit_button("Filter Ground Motions")
            numberRecords = st.number_input("Number of Ground Motions to be Scaled", 1, value=11, step=1)
            selectButton = st.form_submit_button("Find Optimum Selected Ground Motions")

    with scalingTab:
        with st.form("scalingForm"):
            spectralOrdinate = st.selectbox("Spectral Ordinate", ["SRSS", "RotD50", "RotD100"])
            targetShift = st.number_input("Target Spectrum Shift", 0.1, value=1.3, step=0.1)
            periodRange = st.slider("Period Range of Interest Coefficients", 0.1, 3.0, (0.2, 1.5), 0.1)
            scaleButton = st.form_submit_button("Perform Amplitude Scaling")

st.markdown("#### ScalePy")

if responseButton:
    defaultTarget = scalepyBack.targetSpectrum(Ss, S1, soil)
    defaultTarget = defaultTarget.rename(columns={'T': 'Period (sec)', 'Sa': 'pSa (g)'})
    defaultFig = go.Figure()
    defaultFig.add_trace(go.Scatter(x = defaultTarget['Period (sec)'],
                                    y=defaultTarget['pSa (g)'],
                                    name='Response Spectrum', line=dict(color='red')))

    defaultFig.update_xaxes(
        title_text = 'Period (sec)',
        range=[0,4],
        tickvals=np.arange(0,4.5,0.5),
        dtick = 1
    )

    defaultFig.update_yaxes(
        title_text = 'pSa (g)',
        range=[0,3]
    )

    defaultFig.update_layout(showlegend=True, width=1000,height=600, title='Response Spectrum', legend=dict(
    yanchor="top",
    xanchor="right",
    ))

if filterButton:

    def tupleToStr(tup):
        return'{value1} {value2}'.format(value1= tup[0], value2 = tup[1])
    selectedTarget = scalepyBack.targetSpectrum(Ss, S1, soil)
    selected_keys, eqe_selected_x, eqe_selected_y, rsn_selected, t, eqe_s = scalepyBack.recordSelection(tupleToStr(magnitudeRange), 
                                                                                             tupleToStr(vs30Range), 
                                                                                             tupleToStr(rjbRange), 
                                                                                             faultMechanism, 
                                                                                             tupleToStr(duration575Range), 
                                                                                             tupleToStr(duration595Range), 
                                                                                             tupleToStr(ariasIntensity), 
                                                                                             selectedTarget, 
                                                                                             "Any", 
                                                                                             period)
    defaultFig = go.Figure()
    for name in rsn_selected:
        defaultFig.add_trace(go.Scatter(x = t,
                                    y=eqe_selected_x[name], line=dict(color='gray'), showlegend=False))
    defaultFig.add_trace(go.Scatter(x = t,
                                    y=eqe_selected_y[name], line=dict(color='gray'), showlegend=False))

    defaultFig.add_trace(go.Scatter(x = selectedTarget['T'],
                                    y=selectedTarget['Sa'],
                                    name='Response Spectrum', line=dict(color='red')))

    defaultFig.update_layout(showlegend=True, width=1000,height=600, title="Filtered Records", legend=dict(
            yanchor="top",
            x = 1,
            xanchor="right"))

    defaultFig.update_xaxes(
            title_text = 'Period (sec)',
            range=[0,4],
            tickvals=np.arange(0,4.5,0.5),
            dtick = 1)

    defaultFig.update_yaxes(
        title_text = 'pSa (g)',
        range=[0,3])

    defaultFig.add_trace(go.Scatter(
        x = [None],
        y = [None],
        mode = 'lines',
        name = "Filtered Records",
        line=dict(color='gray')
    ))

if selectButton:

    def tupleToStr(tup):
        return'{value1} {value2}'.format(value1= tup[0], value2 = tup[1])
    selectedTarget = scalepyBack.targetSpectrum(Ss, S1, soil)
    selected_keys, eqe_selected_x, eqe_selected_y, rsn_selected, t, eqe_s = scalepyBack.recordSelection(tupleToStr(magnitudeRange), 
                                                                                             tupleToStr(vs30Range), 
                                                                                             tupleToStr(rjbRange), 
                                                                                             faultMechanism, 
                                                                                             tupleToStr(duration575Range), 
                                                                                             tupleToStr(duration595Range), 
                                                                                             tupleToStr(ariasIntensity), 
                                                                                             selectedTarget, 
                                                                                             "Any", 
                                                                                             period,
                                                                                             numberRecords)
    defaultFig = go.Figure()
    for name in selected_keys:
        defaultFig.add_trace(go.Scatter(x = t,
                                    y=eqe_selected_x[name], line=dict(color='gray'), showlegend=False))
    defaultFig.add_trace(go.Scatter(x = t,
                                    y=eqe_selected_y[name], line=dict(color='gray'), showlegend=False))

    defaultFig.add_trace(go.Scatter(x = selectedTarget['T'],
                                    y=selectedTarget['Sa'],
                                    name='Response Spectrum', line=dict(color='red')))

    defaultFig.update_layout(showlegend=True, width=1000,height=600, title="Optimum Selected Records", legend=dict(
            yanchor="top",
            x = 1,
            xanchor="right"))

    defaultFig.update_xaxes(
            title_text = 'Period (sec)',
            range=[0,4],
            tickvals=np.arange(0,4.5,0.5),
            dtick = 1)

    defaultFig.update_yaxes(
        title_text = 'pSa (g)',
        range=[0,3]
        )

    defaultFig.add_trace(go.Scatter(
        x = [None],
        y = [None],
        mode = 'lines',
        name = "Selected Records",
        line=dict(color='gray')
    ))

if scaleButton:

    def tupleToStr(tup):
        return'{value1} {value2}'.format(value1= tup[0], value2 = tup[1])
    selectedTarget = scalepyBack.targetSpectrum(Ss, S1, soil)
    selected_keys, eqe_selected_x, eqe_selected_y, rsn_selected, t, eqe_s = scalepyBack.recordSelection(tupleToStr(magnitudeRange), 
                                                                                             tupleToStr(vs30Range), 
                                                                                             tupleToStr(rjbRange), 
                                                                                             faultMechanism, 
                                                                                             tupleToStr(duration575Range), 
                                                                                             tupleToStr(duration595Range), 
                                                                                             tupleToStr(ariasIntensity), 
                                                                                             selectedTarget, 
                                                                                             "Any", 
                                                                                             period,
                                                                                             numberRecords)
    
    if spectralOrdinate == "SRSS":

        sf_dict, multiplied_selected_x, multiplied_selected_y, geo_mean_1st_scaled_df, srss_mean_df, srss_mean_scaled_df = scalepyBack.amplitudeScaling(selected_keys, selectedTarget, period, targetShift, periodRange[0], periodRange[1], 'srss')

        defaultFig = go.Figure()
        for name in selected_keys:
            defaultFig.add_trace(go.Scatter(x = t,
                                        y= eqe_selected_x[name], line=dict(color='gray'), showlegend=False))
        defaultFig.add_trace(go.Scatter(x = t,
                                        y = eqe_selected_y[name], line=dict(color='gray'), showlegend=False))

        defaultFig.add_trace(go.Scatter(x = selectedTarget['T'],
                                        y = selectedTarget['Sa'],
                                        name='Response Spectrum', line=dict(color='red', dash='dash')))

        defaultFig.add_trace(go.Scatter(x = selectedTarget['T'],
                                        y = selectedTarget['Sa']*targetShift,
                                        name='Shifted Response Spectrum', line=dict(color='red')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = geo_mean_1st_scaled_df[ "Mean"],
                                        name='Geometric Mean Scaled', line=dict(color='blue')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = srss_mean_df[ "Mean"],
                                        name='SRSS Scaled', line=dict(color='black', dash='dash')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = srss_mean_scaled_df[ "Mean"],
                                        name='SRSS Mean Scaled', line=dict(color='black')))

        defaultFig.update_layout(showlegend=True, width=1000,height=600, title="Scaled Ground Motions", legend=dict(
                yanchor="top",
                x = 1,
                xanchor="right"))

        defaultFig.update_xaxes(
                title_text = 'Period (sec)',
                range=[0,4],
                tickvals=np.arange(0,4.5,0.5),
                dtick = 1)

        defaultFig.update_yaxes(
            title_text = 'pSa (g)',
            range=[0,3])

        defaultFig.add_trace(go.Scatter(
            x = [None],
            y = [None],
            mode = 'lines',
            name = "Selected Records",
            line=dict(color='gray')
        ))

        defaultFig.add_vrect(x0=periodRange[0]*period, x1=periodRange[1]*period, 
              fillcolor="yellow", opacity=0.1, line_width=0)

    if spectralOrdinate == "RotD50":
        sf_dict, multiplied_selected_x, multiplied_selected_y, geo_mean_1st_scaled_df, rotd50_mean_df, rotd50_mean_scaled_df = scalepyBack.amplitudeScaling(selected_keys, selectedTarget, period, targetShift, periodRange[0], periodRange[1], 'rotd50')

        defaultFig = go.Figure()
        for name in selected_keys:
            defaultFig.add_trace(go.Scatter(x = t,
                                        y= eqe_selected_x[name], line=dict(color='gray'), showlegend=False))
        defaultFig.add_trace(go.Scatter(x = t,
                                        y = eqe_selected_y[name], line=dict(color='gray'), showlegend=False))

        defaultFig.add_trace(go.Scatter(x = selectedTarget['T'],
                                        y = selectedTarget['Sa'],
                                        name='Response Spectrum', line=dict(color='red', dash='dash')))

        defaultFig.add_trace(go.Scatter(x = selectedTarget['T'],
                                        y = selectedTarget['Sa']*targetShift,
                                        name='Shifted Response Spectrum', line=dict(color='red')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = geo_mean_1st_scaled_df[ "Mean"],
                                        name='Geometric Mean Scaled', line=dict(color='blue')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = rotd50_mean_df[ "Mean"],
                                        name='RotD50 Scaled', line=dict(color='black', dash='dash')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = rotd50_mean_scaled_df[ "Mean"],
                                        name='RotD50 Mean Scaled', line=dict(color='black')))

        defaultFig.update_layout(showlegend=True, width=1000,height=600, title="Scaled Ground Motions", legend=dict(
                yanchor="top",
                x = 1,
                xanchor="right"))

        defaultFig.update_xaxes(
                title_text = 'Period (sec)',
                range=[0,4],
                tickvals=np.arange(0,4.5,0.5),
                dtick = 1)

        defaultFig.update_yaxes(
            title_text = 'pSa (g)',
            range=[0,3])

        defaultFig.add_trace(go.Scatter(
            x = [None],
            y = [None],
            mode = 'lines',
            name = "Selected Records",
            line=dict(color='gray')
        ))

        defaultFig.add_vrect(x0=periodRange[0]*period, x1=periodRange[1]*period, 
              fillcolor="yellow", opacity=0.1, line_width=0)

    if spectralOrdinate == "RotD100":
        sf_dict, multiplied_selected_x, multiplied_selected_y, geo_mean_1st_scaled_df, rotd100_mean_df, rotd100_mean_scaled_df = scalepyBack.amplitudeScaling(selected_keys, selectedTarget, period, targetShift, periodRange[0], periodRange[1], 'rotd100')

        defaultFig = go.Figure()
        for name in selected_keys:
            defaultFig.add_trace(go.Scatter(x = t,
                                        y= eqe_selected_x[name], line=dict(color='gray'), showlegend=False))
        defaultFig.add_trace(go.Scatter(x = t,
                                        y = eqe_selected_y[name], line=dict(color='gray'), showlegend=False))

        defaultFig.add_trace(go.Scatter(x = selectedTarget['T'],
                                        y = selectedTarget['Sa'],
                                        name='Response Spectrum', line=dict(color='red', dash='dash')))

        defaultFig.add_trace(go.Scatter(x = selectedTarget['T'],
                                        y = selectedTarget['Sa']*targetShift,
                                        name='Shifted Response Spectrum', line=dict(color='red')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = geo_mean_1st_scaled_df[ "Mean"],
                                        name='Geometric Mean Scaled', line=dict(color='blue')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = rotd100_mean_df[ "Mean"],
                                        name='RotD100 Scaled', line=dict(color='black', dash='dash')))

        defaultFig.add_trace(go.Scatter(x = t,
                                        y = rotd100_mean_scaled_df[ "Mean"],
                                        name='RotD100 Mean Scaled', line=dict(color='black')))

        defaultFig.update_layout(showlegend=True, width=1000,height=600, title="Scaled Ground Motions", legend=dict(
                yanchor="top",
                x = 1,
                xanchor="right"))

        defaultFig.update_xaxes(
                title_text = 'Period (sec)',
                range=[0,4],
                tickvals=np.arange(0,4.5,0.5),
                dtick = 1)

        defaultFig.update_yaxes(
            title_text = 'pSa (g)',
            range=[0,3])

        defaultFig.add_trace(go.Scatter(
            x = [None],
            y = [None],
            mode = 'lines',
            name = "Selected Records",
            line=dict(color='gray')
        ))
        
        defaultFig.add_vrect(x0=periodRange[0]*period, x1=periodRange[1]*period, 
              fillcolor="yellow", opacity=0.1, line_width=0)

st.plotly_chart(defaultFig)

if selectButton:
    defaultFrame = pd.DataFrame(columns=["Record Sequence Number", "Earthquake Name", "Station Name", "Scale Factor"], index=pd.RangeIndex(start=1, stop = numberRecords + 1, name='index'))
    eqe_s_filtered = eqe_s[ eqe_s["RecordSequenceNumber"].isin(selected_keys)]
    defaultFrame["Earthquake Name"] = eqe_s_filtered["EarthquakeName"].to_list()
    defaultFrame["Station Name"] = eqe_s_filtered["StationName"].to_list()
    defaultFrame["Record Sequence Number"] = selected_keys
    st.table(defaultFrame)

if scaleButton:

    defaultFrame = pd.DataFrame(columns=["Record Sequence Number", "Earthquake Name", "Station Name", "Scale Factor"], index=pd.RangeIndex(start=1, stop = numberRecords + 1, name='index'))
    eqe_s_filtered = eqe_s[ eqe_s["RecordSequenceNumber"].isin(selected_keys)]
    defaultFrame["Earthquake Name"] = eqe_s_filtered["EarthquakeName"].to_list()
    defaultFrame["Station Name"] = eqe_s_filtered["StationName"].to_list()
    defaultFrame["Record Sequence Number"] = selected_keys

    if spectralOrdinate == "SRSS":

        sf_dict, multiplied_selected_x, multiplied_selected_y, geo_mean_1st_scaled_df, srss_mean_df, srss_mean_scaled_df = scalepyBack.amplitudeScaling(selected_keys, selectedTarget, period, targetShift, periodRange[0], periodRange[1], 'srss')
        defaultFrame["Scale Factor"] = list(sf_dict.values())
        st.table(defaultFrame)

    if spectralOrdinate == "RotD50":

        sf_dict, multiplied_selected_x, multiplied_selected_y, geo_mean_1st_scaled_df, rotd50_mean_df, rotd50_mean_scaled_df = scalepyBack.amplitudeScaling(selected_keys, selectedTarget, period, targetShift, periodRange[0], periodRange[1], 'rotd50')
        defaultFrame["Scale Factor"] = list(sf_dict.values())
        st.table(defaultFrame)

    if spectralOrdinate == "RotD100":
        sf_dict, multiplied_selected_x, multiplied_selected_y, geo_mean_1st_scaled_df, rotd100_mean_df, rotd100_mean_scaled_df = scalepyBack.amplitudeScaling(selected_keys, selectedTarget, period, targetShift, periodRange[0], periodRange[1], 'rotd100')
        defaultFrame["Scale Factor"] = list(sf_dict.values())
        st.table(defaultFrame)



