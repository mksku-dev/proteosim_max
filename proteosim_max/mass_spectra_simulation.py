amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    Takes a peptide seuqence and calculates the molecular mass of the sequence 

    Parameters
    ----------
    peptide_seq : str
        Sequene of the peptide (one letter code)
    amino_acid_mass_dict : dict
        Is a dictonary which gives the molecular mass of every aminoacid

    Returns
    -------
    peptide_mass : dict
        Returns the peptide and the molecular mass of the petide sequence 
    """
    peptide_mass = 0
    for x in peptide_seq:
        if x in amino_acid_mass_dict:
            peptide_mass += amino_acid_mass_dict[x]
    
    mass_peptide_dic = {peptide_seq : peptide_mass}
    
    return mass_peptide_dic

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Calculates the molecular mass of multiple peptides in a list 
    
    Parameters
    ----------
    peptides : list
        List of peptide sequences 
    amino_acid_mass_dict : dict

    Returns
    -------
    peptide_mass : dict
        Returns a dictonary of the peptide sequences and the corressponding molecular mass

    """

    mass_peptide_map = {}

    for pep in peptides:
        peptide_mass = 0
        for x in pep:
            if x in amino_acid_mass_dict:
                peptide_mass += amino_acid_mass_dict[x]
        mass_peptide_map[pep] = peptide_mass
    
    return mass_peptide_map



def calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007):
    """
    Calculates the mass to charge ratio.  

    Parameters
    ----------
    peptide_mass_map : dict 
        Dictonary of peptide sequences with their corresponding molecular mass 
    charge : int
        Charge state used for all peptides?
    proton_mass : float
        The mass of a proton. default = 1.007, should not be changed

    Returns
    -------
    mz_map : dict
        Retruns a new dictonary paairing peptides with their mz-values 

    """
    mz_map = {}

    for pep, peptide_mass in peptide_mass_map.items():
        mz_value = (peptide_mass + charge * proton_mass) / charge
        mz_map[pep] = mz_value # man kann potentiell noch die masse der peptide mit rein machen 

    return mz_map



import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Generates bar charts.

    Parameters
    ----------
    mz_values : dict
        Dictionary mapping peptides to their m/z values
    random_count_range : tuple of int 
        Setting a range of possible intensities (default 0 till 30000)
    seed : int 
        Makes sure that the intensity values are always the same
    Returns
    -------
    none
    """
    np.random.seed(seed)
    mz_list = list(mz_values.values())
    intensities = np.random.randint(random_count_range[0], random_count_range[1], size=len(mz_list))
    plt.bar(mz_list, intensities, width = 1)
    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.show()



def fragment_peptide(peptide):
    """
    Spells the petide sequence in both directions and makes it muitple fragments

    Parameters
    ----------
    peptide : str
        Peptide sequence

    Returns
    -------
    peptide_fragements : list of str
        Returns the b- and y-ion fragments of a peptide sequence
    """
    b_ions = []
    y_ions = []

    for i in range(1, len(peptide)): 
        b_ions.append(peptide[:i])
        y_ions.append(peptide[-i:])

    combinded_list = b_ions + y_ions

    return combinded_list



