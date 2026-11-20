def test_fasta_reader():
    tmp_fasta_path = './data/dummy-proteins.fasta'
    import proteosim_max as ps
    protein_map = ps.read_fasta(tmp_fasta_path)
    
    
    # Replace the strings with your fasta content
    # which you expect to be now available as a dictionary
    assert protein_map["PROTEIN_ID_1"].startswith("MSISDNSADASBTASN")
    assert protein_map["PROTEIN_ID_2"].startswith("MEDEEKAVEILGNTEAAHPPS")

# test_fasta_reader()