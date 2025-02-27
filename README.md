# TranscriptionDynamicModeling

Repository for "Highly Constrained Kinetic Models for Quantitative Single-Cell Gene Expression Analysis"

## Authors: 
Hyeon Jin Cho1,2 
Christopher H. Borher2 
Pawel Trzaskoma5 
Rob Patro3 
Carson C. Chow4,* 
Daniel R. Larson2,6,* 


## Affiliations:
1. Department of Cell Biology and Molecular Genetics, University of Maryland, College Park, MD, USA
2. Laboratory of Receptor Biology and Gene Expression, National Cancer Institute, NIH, Bethesda, MD, USA
3. Department of Computer Science and Center for Bioinformatics and Computational Biology, University of Maryland, College Park, MD, USA
4. Laboratory of Biological Modeling, National Institute of Diabetes and Digestive and Kidney Diseases, NIH, Bethesda, MD, USA
5. Laboratory of Molecular Immunogenetics, Genomics and Immunity Section, National Institute of Arthritis and Musculoskeletal and Skin Diseases, NIH, Bethesda, MD, USA
6. Lead contact
* Correspondence: Hyeon Jin Cho (hyeonjin.cho@nih.gov)
* Correspondence: Carson C. Chow (carsonc@niddk.nih.gov)
* Correspondence: Daniel R. Larson (dan.larson@nih.gov)

## This repository includes:
1. Code used for preprocessing of raw scRNA-Seq data
   * HBEC scRNA-seq data (10X)
     - Seurat_workflow.R: filtering step using Seurat package (R)
   * TF Atlas Data:
     - 1_TF_Atlas_workflow.py: extract count matrix of TF of interest
     - 2_rawCounts2freq.py: make individual PDF files from count matrix
2. Codes used to generate all plots.
3. Code used to run StochasticGene.jl (using smFISH, bursting data, and TF dwell time data as inputs)
   * StochasticGene(v1.1.7)
     - makeswarm.sh: make swarm file (this will make example.swarm)
     - fitscript_smFISH.jl: adjusted StochasticGene.jl 

All the rates were inferred using StochasticGene v.1.1.7

Julia package StochasticGene can be installed directly from Julia using:
```
Pkg.add(name="StochasticGene", version="1.1.7")
```
and is also available at: https://github.com/nih-niddk-mbs/StochasticGene.jl
