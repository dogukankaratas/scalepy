import streamlit as st

st.set_page_config(page_title="Scalepy GUI", layout="wide")

with st.sidebar:
    st.markdown("# ScalePy")
    st.markdown("Click [here](https://github.com/dogukankaratas/scalepy) for source code.")
    st.markdown("Use [here](https://www.linkedin.com/in/dogukankaratas/) to reach me.")

st.markdown("# Welcome to ScalePy 👋")
st.write("ScalePy is an open-source ground motion selection and scaling framework developed in Python.")
st.write("Package can be use as an API, which is located in GitHub page, and also a GUI has developed for the end-user.")

st.markdown('## Guide')
st.markdown('Follow the steps below to get started.')

st.markdown('### Step 1: Generate Response Spectrum')
st.markdown("First generate a response spectrum with using built-in functions or upload a file.")
col1, col2 = st.columns(2)
with col1:
    st.image('assets/responseSpectrum.png')
with col2:
    pass

st.markdown('### Step 2: Filter Record Database and Select Record Set')
st.markdown('Filter your records acc. to the properties of your case and select an optimum set using ScalePy algorithm.')
col3, col4 = st.columns(2)
with col3:
    st.image('assets/filterRecords.png')
with col4:
    st.image('assets/selectRecords.png')
    
st.markdown('### Step 3: Amplitude Scaling')
st.markdown('Perform amplitude scaling with selected spectral ordinate.')
col5, col6 = st.columns(2)
with col5:
    st.image('assets/scaleRecords.png')
with col6:
    st.image('assets/rotScale.png')
    
st.markdown('### Results')
st.markdown('For every step, your data will be visualized in the chart.')
col7, col8 = st.columns(2)
with col7:
    st.image('assets/optimumRecords.png')
with col8:
    st.image('assets/scaledRecordsChart.png')

st.markdown("Also you can read your scale factors and selected records from the table.")
st.image('assets/resultTable.png')