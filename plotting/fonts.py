import matplotlib.font_manager as fm
from pyfonts import load_google_font
import streamlit as st

PUBLICATION_FONTS = ["Arial", "Helvetica", "Times New Roman"]

@st.cache_resource
def get_available_fonts():
    fm._load_fontmanager(try_read_cache=False)
    all_fonts = sorted(set([f.name for f in fm.fontManager.ttflist]))
    publication_fonts, missing_fonts = [], []
    google_fonts_map = {
        "Arial": "Arimo",
        "Helvetica": "Roboto",
        "Times New Roman": "Tinos"
    }
    for pub_font in PUBLICATION_FONTS:
        found = [f for f in all_fonts if pub_font.lower() in f.lower()]
        if found:
            publication_fonts.append(found[0])
        else:
            gf = google_fonts_map.get(pub_font)
            if gf:
                try:
                    path = load_google_font(gf)
                    fm.fontManager.addfont(path)
                    publication_fonts.append(gf)
                except Exception:
                    missing_fonts.append(pub_font)
    fm._load_fontmanager(try_read_cache=False)
    if not publication_fonts:
        publication_fonts = ["DejaVu Sans"]
    return sorted(set(publication_fonts)), missing_fonts