import re

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}


def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digest a protein sequence into peptides using a cleavage pattern
    and keep only peptides in a  certain length window.

    Parameters
    ----------
    protein_seq : str
        Sequence of the protein.
    cleave_pattern : str
        Sequence where the Enzyme.
    min_pep_len : int, optional
        Minimum peptide length to keep (default 5).
    max_pep_len : int, optional
        Maximum peptide length to keep (default 30).

    Returns
    -------
    list of str
        Peptides within the length window.
    """
    if not protein_seq:
        return []
    if min_pep_len > max_pep_len:
        raise ValueError("min_pep_len must not be larger than max_pep_len")

    peptides = re.split(cleave_pattern, protein_seq)
    filtered_peptides = [pep for pep in peptides if min_pep_len <= len(pep) <= max_pep_len]
    return filtered_peptides
    # bei sehr großen datasets etc würde des hier viel zu lang dauern -> lambda funktion auf die map 


def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digest every protein in a dictionary into peptides.

    Parameters
    ----------
    protein_map : dict
        Dictionary mapping protein IDs (str) to amino-acid sequences (str).
    cleave_pattern : str
        Regular expression describing where the enzyme cuts.
    min_pep_len : int, optional
        Minimum peptide length to keep (default 5).
    max_pep_len : int, optional
        Maximum peptide length to keep (default 30).

    Returns
    -------
    dict
        Dictionary mapping each protein ID (str) to a list of peptides (list of str).
    """
    digested_map = {}
    for protein_id, sequence in protein_map.items():
        digested_map[protein_id] = digest_protein_sequence(             # bisschen verwirrend ig in 2 teile aufteilen, geht aber auch so 
            sequence, cleave_pattern, min_pep_len, max_pep_len
        )
    return digested_map



def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute the sequence coverage of a protein given a list of peptides.

    A residue counts as covered if it is part of at least one peptide
    that matches the protein sequence.

    Parameters
    ----------
    protein_seq : str
        Amino-acid sequence of the protein.
    peptides : list of str
        Peptides obtained from the digestion of this protein.

    Returns
    -------
    float
        Percentage of residues covered by at least one peptide (0-100).
    """
    if not protein_seq:
        return 0.0

    covered = [False] * len(protein_seq)

    for pep in peptides:
        start = protein_seq.find(pep)
        while start != -1:
            for i in range(start, start + len(pep)):
                covered[i] = True
            start = protein_seq.find(pep, start + 1)

    return sum(covered) / len(protein_seq) * 100