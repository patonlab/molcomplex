          ___           ___                         ___           ___           ___           ___                       ___           ___      
         /\  \         /\  \                       /\__\         /\  \         /\  \         /\  \                     /\__\         /|  |     
        |::\  \       /::\  \                     /:/  /        /::\  \       |::\  \       /::\  \                   /:/ _/_       |:|  |     
        |:|:\  \     /:/\:\  \                   /:/  /        /:/\:\  \      |:|:\  \     /:/\:\__\                 /:/ /\__\      |:|  |     
      __|:|\:\  \   /:/  \:\  \   ___     ___   /:/  /  ___   /:/  \:\  \   __|:|\:\  \   /:/ /:/  /  ___     ___   /:/ /:/ _/_   __|:|__|     
     /::::|_\:\__\ /:/__/ \:\__\ /\  \   /\__\ /:/__/  /\__\ /:/__/ \:\__\ /::::|_\:\__\ /:/_/:/  /  /\  \   /\__\ /:/_/:/ /\__\ /::::\__\_____
     \:\~~\  \/__/ \:\  \ /:/  / \:\  \ /:/  / \:\  \ /:/  / \:\  \ /:/  / \:\~~\  \/__/ \:\/:/  /   \:\  \ /:/  / \:\/:/ /:/  / ~~~~\::::/___/
      \:\  \        \:\  /:/  /   \:\  /:/  /   \:\  /:/  /   \:\  /:/  /   \:\  \        \::/__/     \:\  /:/  /   \::/_/:/  /      |:|~~|    
       \:\  \        \:\/:/  /     \:\/:/  /     \:\/:/  /     \:\/:/  /     \:\  \        \:\  \      \:\/:/  /     \:\/:/  /       |:|  |    
        \:\__\        \::/  /       \::/  /       \::/  /       \::/  /       \:\__\        \:\__\      \::/  /       \::/  /        |:|__|    
         \/__/         \/__/         \/__/         \/__/         \/__/         \/__/         \/__/       \/__/         \/__/         |/__/     
                                                                                                                                                                                                                                                                                  
Implementing a variety of complementary metrics for molecular complexity and synthetic accessibility.

A collaboration with the Sarpong group to understand complexity of molecules

# Requirements
- numpy, pandas
- rdkit
- openbabel
- mordred 
- SYBA (conda install -c lich syba)

## Set up conda environment directly using the yml file:
To install the required packages through Conda, use the env.yml file as follows and the activate the environment: 
1. `conda env create -f env.yml`
2. `conda activate mc`
This will set up the environment with molcomplex installed.

## For installation by cloning the GitHub folder, perform the follwoing steps:
1. Download the zipped folder or clone using: `git clone https://github.com/patonlab/molcomplex.git`  
2. Navigate to the installed folder and run: `python setup.py install`. This will install `molcomplex` in the environment you are present in. 
3. Install necessary dependencies using the following: `conda install -c lich syba`, `conda install -c conda-forge rdkit`, and `conda install -c conda-forge openbabel`

## Recommended installation and update guide (under works)
In a nutshell, `molcomplex` and its dependencies are installed/updated as follows:  
1. Install using conda-forge: `conda install -c conda-forge molcomplex`  
2. Update to the latest version: `pip install molcomplex --upgrade` 

## Usage
To display the options type:

``python -m molcomplex -h``

The `molcomplex` package can be utilised as follows to obtain a csv with complexity scores.

``python -m molcomplex -f examples/test.txt``

To write to CSV add in the following:

``python -m molcomplex -f examples/test.txt --csv``

To perform a retro analysis by breaking down bonds to get complexity scores for precursors of the input SMILES add the following option:

``python -m molcomplex -f examples/test.txt --csv --retro``

## Usage APP
To run the web app perform the following steps:
1. Navigate to the webapp folder: `cd mcwebapp`
2. Run the app as follows: `python molcomplexapp.py`
3. copy paste the `http://127.0.0.1:8050/` or similar into web browser to utilise as an app.


# Metrics implemented
- Bertz Complexity (CT) Score (JACS 1981, 103, 3241-3243)
- Balaban J Score (Chem. Phys. Lett. 1982, 89, 399-404)
- Coley SCScore (J. Chem. Inf. Model. 2018, 58, 2, 252)
- IPC: Bonchev & Trinajstic's information content of the coefficients of the characteristic polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule (J. Chem. Phys. 1977, 67, 4517-4533)
- Ertl SA_Score (J. Cheminform. 2009, 1, 8)
- Boettcher Score (J. Chem. Inf. Model. 2016, 56, 3, 462–470)
- Rücker's total walk count (twc) index: Rücker, G.; Rücker, C. Counts of All Walks as Atomic and Molecular Descriptors. (J. Chem. Inf. Comput. Sci. 1993, 33, 683-695)
- Proudfoot's Cm index based on atom environments: Proudfoot, J. R. A path based approach to assessing molecular complexity. Bioorganic Med. Chem. Lett. 27, 2014–2017 (2017)
- Kappa Shape Indices 1, 2 & 3 (Quant. Struct. Act. Relat. 1986, 5, 1-7)
- McGowan Volume (Chromatographia, 1987, 23, 243-246)
- Labute Approximate Surface Area (Methods Mol Biol 2004, 275, 261-78)
- Van der Waals Volume Atomic and Bond Contributions (J. Org. Chem. 2003, 68, 7368-7373).
- Zagreb Index 
- MOE Type Desciptors (Labute ASA, PEOE VSA, SMR VSA, SLogP VSA)
- SYBA Score (J. Cheminformatics 2020, 12, 35)
- Multiple additional 2D metrics.

# To do list
- compare against human metric
- comment out descriptors that are highly correlated (say > 0.9)

# Currently broken
- Kier's alpha-modified shape indices

# Metrics to implement:

- Bertz’s Ns and Nt index: Bertz, S. H. & Sommer, T. J. Rigorous mathematical approaches to strategic bonds and synthetic analysis based on conceptually simple new complexity indices. Chem. Commun. 16, 2409–2410 (1997).

- Randić's zeta index: Randić, M. & Plavšić, D. Characterization of molecular complexity. Int. J. Quantum Chem. 91, 20–31 (2002).

- https://www.nature.com/articles/s41598-018-37253-8

Two noteworthy substructure-based methods are:
- Barone, R. & Chanon, M. A new and simple approach to chemical complexity. Application to the synthesis of natural products. J. Chem. Inf. Comput. Sci. 41, 269–272 (2001).

- Whitlock, H. W. On the structure of total synthesis of complex natural products. J. Org. Chem. 63, 7982–7989 (1998).
