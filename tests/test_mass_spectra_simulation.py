import proteosim_max as ps
amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def test_calculate_mol_mass():
    aa_mass_dict = amino_acid_mass_dalton
    test_peptide = "AA"
    expected_mass = 142.16
    calculated_mass = ps.calculate_mol_mass(test_peptide, aa_mass_dict) 
    assert calculated_mass[test_peptide] == expected_mass

test_calculate_mol_mass()

def test_calculate_mol_mass_collection():
    aa_mass_dict = amino_acid_mass_dalton
    peptides = ['AA', 'PEP']
    expected = {"AA" : 142.16 ,"PEP" :323.36}
    actual = ps.calculate_mol_mass_collection(peptides ,aa_mass_dict)

    assert actual == expected

test_calculate_mol_mass_collection()


def test_calculate_mz_collection():
    peptide_mass_map = {'AA': 142.16}
    actual = ps.calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {"AA" : 72.087}

    assert actual == expected

test_calculate_mz_collection()



def test_fragment_peptide():
    peptide = "PEPT"
    expected = ["P", "PE", "PEP", "T", "PT", "EPT"]
    actual = ps.fragment_peptide(peptide)

    assert set(actual) == set(expected)  #set sortiert fragemnte -> so das es funktioniert egal wie die funktion oben deklariert ist 
    assert "PEPT" not in expected  # which string should not be in the expected?

test_fragment_peptide()