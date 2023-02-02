import streamlit as st
import base64

st.set_page_config(page_title="Scalepy GUI", layout="wide")

with st.sidebar:
    st.markdown("# ScalePy")
    st.markdown("💻 Click [here](https://github.com/dogukankaratas/scalepy) for source code.")
    st.markdown("🧑‍💻 Use [here](https://www.linkedin.com/in/dogukankaratas/) to reach me.")

st.markdown("# Welcome to ScalePy 👋")
st.write("ScalePy is an open-source ground motion selection and scaling framework developed in Python.")
st.write("Package can be use as an API, which is located in GitHub page, and also a GUI has developed for the end-user.")

st.markdown('## Response Spectrum Definition')
st.markdown('Define your target spectrum via using one of the functions of ScalePy below.')

st.image("https://media.giphy.com/media/o0FcopcmojZRYw0lrK/giphy.gif")

st.markdown('## Filter Record Database')
st.markdown('Use the given parameters to filter ground motion database records to find similiar ground motions with your case.')

st.image("https://media.giphy.com/media/eX8k3hxYu5pDv2Ubjg/giphy.gif")

st.markdown('## Find Optimum Set')
st.markdown('Use ScalePy similarity algorithm to find optimum ground motion data set.')

st.image("https://media.giphy.com/media/31KMao2AiZO5TpQAso/giphy.gif")

st.markdown('## Perform Amplitude Scaling')
st.markdown('Use the user inputs to perform amplitude scaling and find optimum scale factors for your set.')

st.image("https://media.giphy.com/media/O8Af6Yivtoi9g2u46W/giphy.gif")