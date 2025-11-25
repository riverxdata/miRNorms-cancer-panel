import streamlit as st

st.set_page_config(page_title="Contact & Citation", page_icon="📧", layout="wide")

st.title("📧 Contact & Citation")

st.markdown("---")

# Citation Section
st.header("📰 How to Cite")

st.markdown(
    """
If you use **RIVER-VIS** in your research, please cite our preprint:

### BioRxiv Preprint
"""
)

# Citation box
st.info(
    """
**Citation Text:**

```
Rapid Integration and Visualization for Enhanced Research (RIVER) Platform: 
A Data Platform for bioinformatics applications. bioRxiv (2024).
```

**BibTeX:**

```bibtex
@article{RIVER2024,
  title={Rapid Integration and Visualization for Enhanced Research (RIVER) Platform: 
         A Data Platform for bioinformatics applications},
  author={Loi, Luu-Phuc and Nguyen, Thanh-Giang Tan},
  journal={bioRxiv},
  year={2024},
  note={Preprint}
}
```

**Full Citation (APA):**

```
Loi, L. P., & Nguyen, T. G. T. (2024). Rapid Integration and Visualization for Enhanced Research 
(RIVER) Platform: A Data Platform for bioinformatics applications. bioRxiv.
```

**Full Citation (Chicago):**

```
Loi, Luu-Phuc, and Thanh-Giang Tan Nguyen. 2024. "Rapid Integration and Visualization for Enhanced 
Research (RIVER) Platform: A Data Platform for bioinformatics applications." bioRxiv.
```
    """
)

st.markdown("---")

# Contact Section
st.header("📬 Contact Information")

st.markdown(
    """
For questions, feedback, or collaboration inquiries, please reach out to us:
"""
)

col1, col2 = st.columns(2)

with col1:
    st.markdown(
        """
    ### 👨‍⚕️ Primary Contact
    
    **Dr. Luu-Phuc Loi**
    - Position: Head of the Research Office
    - Institution: Institute for Applied Research in Health Sciences and Aging
    - Organization: Thong Nhat Hospital
    - Email: [luu.p.loi@googlemail.com](mailto:luu.p.loi@googlemail.com)
    - Role: Project Lead & Principal Investigator
    """
    )

with col2:
    st.markdown(
        """
    ### 👨‍💻 Secondary Contact
    
    **Mr. Thanh-Giang Tan Nguyen**
    - Position: Founder & Lead Developer
    - Platform: RIVERXDATA
    - Email: [nttg8100@gmail.com](mailto:nttg8100@gmail.com)
    - Role: Owner & Developer
    """
    )

st.markdown("---")

# Support Section
st.header("💬 Support & Feedback")

st.markdown(
    """
### Report Issues & Get Support
For bugs, issues, or technical support:

1. **Email Dr. Loi:** [luu.p.loi@googlemail.com](mailto:luu.p.loi@googlemail.com)
2. **Email Mr. Nguyen:** [nttg8100@gmail.com](mailto:nttg8100@gmail.com)
3. **GitHub Issues:** [Submit an issue on GitHub](#)

### Feature Requests
Have a feature request or suggestion?

- **Email:** Send your ideas to either contact
- **GitHub Discussions:** [Post your idea](#)

### Collaboration & Research Inquiries
Interested in collaboration or using RIVER-VIS for your research?

**Contact:** Dr. Luu-Phuc Loi
- Email: [luu.p.loi@googlemail.com](mailto:luu.p.loi@googlemail.com)
- Institution: Institute for Applied Research in Health Sciences and Aging, Thong Nhat Hospital

**Technical Implementation:** Mr. Thanh-Giang Tan Nguyen
- Email: [nttg8100@gmail.com](mailto:nttg8100@gmail.com)
- Platform: RIVERXDATA
"""
)

st.markdown("---")

# Acknowledgments Section
st.header("🙏 Acknowledgments")

st.markdown(
    """
### Project Leadership
- **Dr. Luu-Phuc Loi** - Principal Investigator & Research Director
- **Mr. Thanh-Giang Tan Nguyen** - Founder & Development Lead

### Institutions
- Institute for Applied Research in Health Sciences and Aging
- Thong Nhat Hospital
- RIVERXDATA Platform

### Technologies & Libraries
RIVER-VIS is built using:
- **Streamlit** - Interactive web app framework
- **Pandas** - Data manipulation
- **Plotly** - Interactive visualizations
- **SciPy** - Statistical analysis
- **Matplotlib & Seaborn** - Publication-quality plots
- **NumPy** - Numerical computing
- **Scikit-learn** - Machine learning utilities

### Data Sources
- **GEO Database (GSE211692)** - Serum microRNA profiles of 9,921 patients with 13 types of human solid cancers
- NCBI Gene Expression Omnibus (GEO)

### Related Resources
- [Streamlit Documentation](https://docs.streamlit.io/)
- [SciPy Statistics](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [GEO Database](https://www.ncbi.nlm.nih.gov/geo/)
- [BioRxiv](https://www.biorxiv.org/)
"""
)

st.markdown("---")

# Footer
st.markdown(
    """
### 📄 Additional Information

- **Current Version:** 0.1.0 (Beta)
- **Last Updated:** November 2024
- **Status:** Pre-publication (bioRxiv Preprint)
- **Project:** Rapid Integration and Visualization for Enhanced Research (RIVER) Platform

### Citation Status
This is a **preprint** available on **bioRxiv**. Please cite accordingly and note that the work 
is currently under peer review for a full journal publication.

---

*RIVER-VIS is a project of the Institute for Applied Research in Health Sciences and Aging at 
Thong Nhat Hospital, developed under the RIVERXDATA Platform. All rights reserved.*
"""
)

st.markdown("---")
