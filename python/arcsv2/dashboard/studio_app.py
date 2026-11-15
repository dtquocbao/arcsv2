"""ARC-SV Studio — guided app (U6). Launch via ``arcsv2 dashboard``."""

from __future__ import annotations

import os

import streamlit as st

from arcsv2 import __version__
from arcsv2.dashboard import studio_ui as ui

st.set_page_config(page_title="ARC-SV Studio", page_icon="🧬", layout="wide")

nav = st.navigation(
    [
        st.Page(ui.page_example, title="Try example", icon="🚀", default=True),
        st.Page(ui.page_setup, title="Setup check", icon="🩺"),
        st.Page(ui.page_new_analysis, title="New analysis", icon="▶️"),
        st.Page(ui.page_analyses, title="My analyses", icon="📋"),
        st.Page(ui.page_results, title="Results", icon="📊"),
        st.Page(ui.page_help, title="Help", icon="❓"),
    ]
)

with st.sidebar:
    st.title("ARC-SV Studio")
    st.caption(f"Version {__version__}")
    root = os.environ.get("ARCSV2_RESULTS", "").strip() or str(ui.results_path())
    st.markdown(f"**Data folder**  \n`{root}`")
    st.caption("Developer view: `arcsv2 dashboard --legacy`")

nav.run()
