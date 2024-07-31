# ProStruct
ProStruct: A Python-based Tool for Homology Modeling and 3D Structure Prediction
## Abstract:
Homology modeling is a widely used computational technique for predicting the three-dimensional (3D) structures of proteins based on known templates, leveraging evolutionary relationships to provide structural insights critical for understanding protein function, interactions, and potential therapeutic targets. However, existing tools often require significant expertise and computational resources, presenting a barrier for many researchers. This project aims to develop a user-friendly, automated homology modeling tool that streamlines the process of protein structure prediction by integrating sequence alignment, template selection, and preliminary structure validation into a cohesive workflow. Utilizing Biopython for sequence alignment and BLAST searches to identify suitable template structures from the Protein Data Bank (PDB), the tool features automated template selection based on alignment scores, customizable alignment parameters, and a modular design that allows integration with additional structure prediction and validation tools. The pipeline is designed to efficiently handle both single and multiple sequence inputs. This tool simplifies the homology modeling process, making it accessible to researchers with varying levels of computational expertise, enabling high-throughput structural biology research, and facilitating the exploration of protein functions and interactions. By providing a flexible framework and user-friendly interface, the tool represents a significant advancement in the accessibility and efficiency of protein structure prediction. Its open-source nature encourages community collaboration and continuous improvement, ensuring it remains a cutting-edge solution for structural biology research.
### Keywords: Homology modeling, protein structure prediction, Biopython, BLAST, structural open-source tool.

### Objectives
Develop an Automated Tool: Create a homology modeling tool that automates the key steps of protein structure prediction.

User-Friendly Interface: Design the tool to be easily accessible and usable by researchers with varying levels of computational expertise.

Workflow: Integrate sequence alignment, template selection, and preliminary structure validation into a single, cohesive workflow.

Utilize Open-Source Libraries: Leverage open-source libraries such as Biopython to ensure flexibility, accessibility, and scalability.

Facilitate High-Throughput Research: Enable the tool to handle multiple sequences efficiently, supporting large-scale structural biology studies.

Enhance Structural Insights: Provide accurate and reliable protein structure predictions to aid in understanding protein functions, interactions, and potential therapeutic targets.

Make the tool open-source to foster collaboration, continuous improvement, and the integration of additional features by the research community.

### Workflow
![image](https://github.com/omicscodeathon/prostruct/blob/main/workflow/WhatsApp%20Image%202024-07-29%20at%2001.31.46.jpeg)
Figure1.workflow



### Methods :

### 1. Tool Design and Architecture:
The Prostruc tool combines sequence alignment, template selection, model creation, and validation to expedite homology modeling of protein structures. The tool was implemented using Biopython for sequence handling and integration with Rosetta and PyRosetta for structural modeling tasks.

### 2. Sequence Retrieval and Input:
Protein sequences were retrieved from databases such as UniProt or provided by users in FASTA format via the tool’s command-line interface (CLI). Input sequences underwent initial validation for correct format and completeness using Biopython's SeqIO module.

### 3. Sequence Alignment:
Sequence alignment was performed using BLAST for pairwise alignments and the MSA (Multiple Sequence Alignment) module for multiple sequence alignments. Alignment parameters were optimized for sequence identity thresholds, E-value and alignment quality to ensure accurate template identification.

### 4. Template Selection:
Templates for homology modeling were selected based on sequence similarity and structural relevance using PyRosetta's homology modeling protocols. Templates were identified from the Protein Data Bank (PDB) using BLAST and aligned sequences.

### 5. Building Model:
PyRosetta’s comparative modeling tools were employed for initial model generation, incorporating aligned sequences and selected templates.


##  Team Members

### 1. Shivani Pawar
### 2. Wilson Sena Kwaku Banini
### 3. Musa Muhammad Shamsuddeen
### 4. Toheeb Jumah
### 5. Nigel Dolling
### 6. Abdulwasiu Tiamiyu
