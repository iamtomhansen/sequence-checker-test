# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 21:28:07 2021
@author: thomas
"""
#imports
# import pandas as pd # for future use of dataframes
import streamlit as st
# import base64 # for future downloads to csv/xlsx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
# from PIL import Image # for future images

# Page config
st.set_page_config(page_title='Sequence data checker', page_icon='🔬', layout='centered', initial_sidebar_state='auto')

# Page title
st.title('Sequence data checker :microscope:')
st.markdown('''
            Last edited on June 29 2021 by Thomas Hansen. Contact 21th14 at queensu dot ca
            
            For now, this is a personal learning experience, not a productivity tool. Use at your own risk.
            ***
            ''')

# Functions: Trim sequence
def trimstart(seq):
    i=0
    j=0
    while i < (len(seq)-1) and j < 100: # Scan through the chain until the end or until 100 non-N residues are passed
        if seq[i]=='N':
            j=0
        i+=1
        j+=1
        if i==len(seq): # In case there are no N's in the sequence
            i=0
            j=0
    return int(i-j+1)

def trimend(seq,start):
    i=seq.find('N',start)
    if i==-1:
        i=len(seq)
    return i

def seq_filetype(file):
    filetypes={
        '.abi':'abi',
        '.ab1':'abi',
        '.aln':'clustal',
        '.fa':'fasta',
        '.fasta':'fasta',
        '.fna':'fasta',
        '.fastq':'fastq',
        '.fq':'fastq',
        '.embl':'genbank',
        '.gb':'genbank',
        '.gbk':'genbank',
        '.genbank':'genbank',
        '.scf':'snapgene'
    }
    dotposition = file.name.rfind('.')
    filetype = filetypes.get(file.name[dotposition:len(file.name)],'unknown')
    return filetype

def query_aligner(query,seq_alignto):
    try: # Future: Could support protein queries (and make it force uppercase if I do)
        if query_type=='Nucleic':
            query=Seq(query).translate() # translate if it's nucleic
        score=seq_alignto.count_overlap(query)
    except:
        score='?'
    return score

# Sidebar for uploading
with st.sidebar.header('Sequencing data goes here'):
    uploaded_files = st.sidebar.file_uploader("Only abi, fasta, genbank, and snapgene files work for now.", accept_multiple_files=True)
    query=st.sidebar.text_input("Query (optional)")
if query:
    query_type = st.sidebar.radio('',['Nucleic','Amino'])
    filter_query = st.sidebar.checkbox('Exclude mismatches',value=False)
 
# Main panel

if st.sidebar.button('Submit'): # On clicking the button
    for file in uploaded_files: # For each file
        # Determine filetype in the future
        filetype=seq_filetype(file)
        if filetype=='unknown': # skip the file if it can't handle it
            st.header("Sorry, couldn't read this file format: "+file.name)
            st.write('---')
            continue
        for record in SeqIO.parse(file,filetype): # For each record in the file
            st.header("Filename: "+file.name)
            first_nt = trimstart(record.seq) # Get through the N's at the start
            last_nt = trimend(record.seq,first_nt) # First N after the resolved region
            trimmed = record.seq[first_nt:last_nt]
            st.write("Trimmed length:",str(len(trimmed)),", Start:",str(first_nt),", End:",str(last_nt))
            st.text_area("Trimmed forward (5'→3') sequence:",trimmed, height=50)
            st.text_area("Trimmed reverse (3'→5') complement:",trimmed.reverse_complement(), height=50)
            for frame_start in range(3):
                # translate the forward strand
                translated=trimmed[frame_start:].translate()
                if query: # is query is not empty, will return True
                    score=query_aligner(query,translated)
                    if score != 0:
                        st.write("Forward (5'→3') translation in Frame",str(frame_start+1),", Stops:",str(translated.count('*')),", Matches:",str(score))
                        st.write(translated)
                else:
                    st.write("Forward (5'→3') translation in Frame",str(frame_start+1),", Stops:",str(translated.count('*')))
                    st.write(translated)
                # Now the reverse complement
                translated=trimmed.reverse_complement()[frame_start:].translate()
                if query:
                    score=query_aligner(query,translated)
                    if score != 0:
                        st.write("Reverse (3'→5') translation in Frame",str(frame_start+1),", Stops:",str(translated.count('*')),", Matches:",str(score))
                        st.write(translated)
                else:
                    st.write("Reverse (3'→5') translation in Frame",str(frame_start+1),", Stops:",str(translated.count('*')))
                    st.write(translated)
        st.write('---') # Separate files
else:
    st.info('Waiting for sequence data') # Before the button is clicked