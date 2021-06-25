# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 21:28:07 2021

@author: thomas
"""
#imports
import pandas as pd
import streamlit as st
import base64
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from PIL import Image

# Page title
st.title('Sequence data checker :microscope:')
st.markdown('''
            Last edited on June 25 2021 by Thomas Hansen. Contact 21th14 at queensu dot ca
            
            For now, this is a personal learning experience, not a productivity tool. Use at your own risk.
            ***
            ''')

# Functions: Trim sequence
def trimstart(seq):
    i=0
    j=0
    while i < (len(seq)-1) and j < 100:
        if seq[i]=='N':
            j=0
        i+=1
        j+=1
    return int(i-j+1)

def trimend(seq,start):
    i=seq.find('N',start)
    if i==-1:
        i=len(seq)
    return i

# Sidebar for uploading
with st.sidebar.header('Sequencing data goes here'):
    uploaded_files = st.sidebar.file_uploader("Only .ab1 for now", accept_multiple_files=True)
    # st.sidebar.markdown("""Query (optional) <- future plan""")
    query=st.sidebar.text_input("Query (optional) DOES NOTHING FOR NOW")

# Main panel
if st.sidebar.button('Submit'):
    for file in uploaded_files:
        # Determine filetype in the future
        filetype='abi'
        for record in SeqIO.parse(file,filetype):
            st.write("Filename:", file.name)
            first_nt = trimstart(record.seq)
            last_nt = trimend(record.seq,first_nt)
            trimmed = record.seq[first_nt:last_nt]
            st.write("Trimmed length:",str(len(trimmed)),", Start:",str(first_nt),", End:",str(last_nt))
            st.text_area("Trimmed forward (5'→3') sequence:",trimmed, height=50)
            st.text_area("Trimmed reverse (3'→5') complement:",trimmed.reverse_complement(), height=50)
            for frame_start in range(3):
                translated=trimmed[frame_start:].translate()
                if query:
                    try:
                        querytranslated=Seq(query).translate()
                        score=pairwise2.align.localxx(translated,querytranslated,score_only=True)
                    except:
                        score='?'
                    st.write("Translation in Frame",str(frame_start+1),", Stops:",str(translated.count('*')),", Identities",str(score))
                else:
                    st.write("Translation in Frame",str(frame_start+1),", Stops:",str(translated.count('*')))
                st.write(translated)
        st.write('---')
else:
    st.info('Waiting for sequence data')