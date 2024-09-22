# Prostruc
> Prostruc: A Python-based Tool for Homology Modeling and 3D Structure Prediction

## Introduction to Prostruc Structural Prediction Tool
<p align="justify"> Homology modeling is a widely used computational technique for predicting the three-dimensional (3D) structures of proteins based on known templates, leveraging evolutionary relationships to provide structural insights critical for understanding protein function, interactions, and potential therapeutic targets. However, existing tools often require significant expertise and computational resources, presenting a barrier for many researchers. This project aims to develop a user-friendly, automated homology modeling tool that streamlines the process of protein structure prediction by integrating sequence alignment, template selection, and preliminary structure validation into a cohesive workflow. Utilizing Biopython for sequence alignment and BLAST searches to identify suitable template structures from the Protein Data Bank (PDB), the tool features automated template selection based on alignment scores, customizable alignment parameters, and a modular design that allows integration with additional structure prediction and validation tools. The pipeline is designed to efficiently handle both single and multiple sequence inputs. This tool simplifies the homology modeling process, making it accessible to researchers with varying levels of computational expertise, enabling high-throughput structural biology research, and facilitating the exploration of protein functions and interactions. By providing a flexible framework and user-friendly interface, the tool represents a significant advancement in the accessibility and efficiency of protein structure prediction. Its open-source nature encourages community collaboration and continuous improvement, ensuring it remains a cutting-edge solution for structural biology research.</p>

## Prostruc Objectives
1. **Develop an Automated Tool:** Create a homology modeling tool that automates the key steps of protein structure prediction.
2. **User-Friendly Interface:** Design the tool to be easily accessible and usable by researchers with varying levels of computational expertise.
3. **Workflow:** Integrate sequence alignment, template selection, and preliminary structure validation into a single, cohesive workflow.
4. **Utilize Open-Source Libraries:** Leverage open-source libraries such as Biopython to ensure flexibility, accessibility, and scalability.
5. **Facilitate High-Throughput Research:** Enable the tool to handle multiple sequences efficiently, supporting large-scale structural biology studies.
6. **Enhance Structural Insights:** Provide accurate and reliable protein structure predictions to aid in understanding protein functions, interactions, and potential therapeutic targets.
7. **Make the tool open-source:** To foster collaboration, continuous improvement, and the integration of additional features by the research community.

## Workflow
![General Workflow](https://github.com/user-attachments/assets/6d2eb1b0-ec3d-4ed6-9d55-82bc0f77b068)

## Pipeline Development Procedure

### 1. Tool Design and Architecture
The Prostruc tool combines sequence alignment, template selection, model creation, and validation to expedite homology modeling of protein structures. The tool was implemented using Biopython for sequence handling and integration with Rosetta and PyRosetta for structural modeling tasks.

### 2. Sequence Retrieval and Input
Protein sequences were retrieved from databases such as UniProt or provided by users in FASTA format via the tool’s command-line interface (CLI). Input sequences underwent initial validation for correct format and completeness using Biopython's SeqIO module.

### 3. Template Selection
Templates for homology modeling were selected based on sequence similarity and structural relevance using PyRosetta's homology modeling protocols. Templates were identified from the Protein Data Bank (PDB) using BLAST and aligned sequences.

### 4. Sequence Alignment
Sequence alignment was performed using BLAST for pairwise alignments and the MSA (Multiple Sequence Alignment) module for multiple sequence alignments. Alignment parameters were optimized for sequence identity thresholds, E-value and alignment quality to ensure accurate template identification.

### 5. Building Model
PyRosetta’s comparative modeling tools were employed for initial model generation, incorporating aligned sequences and selected templates.

### 6. Validation and Optimization
Validation of the Prostruc tool is performed by comparing Prostruc to popular models such as  SWISS MODEL, MODELLER and PRIMO. The metrics analyzed on the various 3D models generated by these tools against the Prostruc models are:
1. QMEANDisCo (QMEAN Distance Constraints)
2. Root Mean Square Deviations (RMSD)
3. DOPE Scores
4. MolProbity
5. TM Scores

## Results
![image](https://github.com/omicscodeathon/prostruct/blob/main/output/template.png)
![image](https://github.com/omicscodeathon/prostruct/blob/main/output/image4.png)
https://github.com/omicscodeathon/prostruct/tree/main/scripts/prostruc#readme 

## Usage
### Web tool
- The web tool can be found at using the following address: http://149.165.154.75:8501

### Python package
- For users who prefer running Prostruc locally or want to integrate it into their workflows, the Python package provides full functionality via a command-line interface.
- Install the Prostruc package via pip:
  ```bash
  pip install prostruc
  
- Example:
  ```bash
  prostruc --sequence "AAAA" --job_name "new_protein" --email "user@example.com"

## Package Requirements

To ensure that Prostruc functions correctly, the following requirements must be met:

- **Python 3.6+**:  
  Prostruc requires Python version 3.6 or above. Ensure that you have the appropriate Python version installed. You can check your Python version with the following command:
  ```bash
  python --version

- **Docker**:
  Docker is necessary for managing the computational workloads, including modeling and validation processes. Make sure Docker is installed and actively running in the background.
  Verify Docker installation and status using:
  ```bash
  docker --version 

- **Internet**:
  An active internet connection is required for Prostruc to perform BLAST searches, retrieve templates, and complete various prediction     tasks.
  
##  Team Members
- [Shivani Pawar](https://github.com/ShivMC): Department of Biotechnology and Bioinformatics, Deogiri College, Auranagabad, Maharashtra, India.
- [Wilson Sena Kwaku Banini](https://github.com/wilson743): Department of Theoretical and Applied Biology, College of Science, Kwame Nkrumah University of Science and Technology, Kumasi, Ghana.            
- [Musa Muhammad Shamsuddeen](https://github.com/Shamss99): Faculty of Health Sciences, Department of Public Health, National Open University of Nigeria.
- [Toheeb Jumah](https://github.com/Toheeb27): School of Collective Intelligence, University Mohammed VI Polytechnic, Rabat, Morocco.
- [Nigel Dolling](https://github.com/NigelDolling): Department of Parasitology, Noguchi Memorial Institute for Medical Research, University of Ghana, Legon.
- [Abdulwasiu Tiamiyu](https://github.com/Tiamiyu1): School of Collective Intelligence, University Mohammed VI Polytechnic, Rabat, Morocco.
- [Olaitan I. Awe](https://github.com/laitanawe): Training officer, ASBCB, Cape Town, South Africa.
  
## Acknowledgment
- [African Society for Bioinformatics and Computational Biology, Cape Town, South Africa](https://www.asbcb.org/)
- [National Institutes of Health (NIH) Office of Data Science Strategy (ODSS)](https://datascience.nih.gov/)
