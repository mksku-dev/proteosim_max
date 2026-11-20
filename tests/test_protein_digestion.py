import re
import proteosim_max as ps

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}

def test_digest_protein_sequence():
    dummy_seq_1 = "AAAAAKGGGGGGKCCCCCCC"   # cut sites after both K
    dummy_seq_2 = "AKGGGGGGKCC"            # contains very short peptides
    dummy_seq_3 = "AAAAAKPGGGGRCCCCC"      # K followed by P, then R

    lysc = enzyme_cleavage_patterns['LysC']
    lysn = enzyme_cleavage_patterns['LysN']
    trypsin = enzyme_cleavage_patterns['Trypsin']

    # 1) LysC cuts after every K
    assert ps.digest_protein_sequence(dummy_seq_1, lysc) == ["AAAAAK", "GGGGGGK", "CCCCCCC"]

    # 2) LysN cuts before every K
    assert ps.digest_protein_sequence(dummy_seq_1, lysn) == ["AAAAA", "KGGGGGG", "KCCCCCCC"]

    # 3) Trypsin doesnt cut at K when P follows, but must cut after R
    assert ps.digest_protein_sequence(dummy_seq_3, trypsin) == ["AAAAAKPGGGGR", "CCCCC"]

    # 4) Min-length filter removes the short peptides "AK" and "CC"
    assert ps.digest_protein_sequence(dummy_seq_2, lysc) == ["GGGGGGK"]

    # 5) Max-length filter removes peptides longer than 6
    assert ps.digest_protein_sequence(dummy_seq_1, lysc, max_pep_len=6) == ["AAAAAK"]

    # 6) No cut site: the whole sequence comes back as one peptide
    assert ps.digest_protein_sequence("AAAAAAAAA", lysc) == ["AAAAAAAAA"]

    # 7) Edge case: empty sequence gives an empty list
    assert ps.digest_protein_sequence("", lysc) == []

    # 8) Edge case: min > max must raise a ValueError
    try:
        ps.digest_protein_sequence(dummy_seq_1, lysc, min_pep_len=10, max_pep_len=5)
        raise AssertionError("ValueError was expected but not raised")
    except ValueError:
        pass

test_digest_protein_sequence()


def test_digest_protein_collection():
    dummy_proteins = {
        "PROT_A": "AAAAAKGGGGGGKCCCCCCC",
        "PROT_B": "AAAAAKPGGGGRCCCCC",
        "PROT_C": "AKCC",
    }

    trypsin = enzyme_cleavage_patterns['Trypsin']
    result = ps.digest_protein_collection(dummy_proteins, trypsin)

    #assert set(result.keys()) == {"PROT_A", "PROT_B", "PROT_C"} # makes sure that only the given IDs are in the result 
    assert result["PROT_A"] == ["AAAAAK", "GGGGGGK", "CCCCCCC"]
    assert result["PROT_B"] == ["AAAAAKPGGGGR", "CCCCC"]
    assert result["PROT_C"] == []

test_digest_protein_collection()


def test_compute_sequence_coverage():
    dummy_prot_seq = "AAAAAKGGGGGGKCCCCCCC"   # Length 20
    dummy_peps = ["AAAAAK", "CCCCCCC"]         # 6 + 7 = 13 postitions which are covered 

    coverage = ps.compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 13 / 20 * 100          # 65.0

test_compute_sequence_coverage()