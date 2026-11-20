def read_fasta(filepath : str) -> dict:
    """
    Add a short description here.
    Reads a FASTA file and returns its sequence as a dicctonary 

    Parameters
    ----------
    filepath : str
        Path to the FASTA file

    Returns
    -------
    dict 
        Dictonary mapping each protein IDs as keys and sequences as values
    """
    

    protein_map = {}                                                    # empty dictonary  
    current_id = None                                                   # gives the starting id (since no protein, no id)
    current_sequence = []                                               #empty list 

    with open(filepath, 'r') as fasta_handle:                            #opens the file (fasta_path) as fasta_handle
        for line in fasta_handle:                                        # makes for loop which goes through every line in fasta_handle
            stripped = line.strip()                                      # deletes spaces at the start of the line and also deltes \n
            if not stripped:                                             # if  the line is not stripped (means empty) than it skips this line and goes to the next 
                continue
            if stripped.startswith('>'):                                 # when stripped and its starts with ">" -> checks if its a header 
                if current_id is not None:                               # looks if the id is None 
                    protein_map[current_id] = ''.join(current_sequence)  # gives the dictonary the sequence with the current id 
                    current_sequence = []                                # emptys the list again 
                current_id = stripped.split('|')[1]                      # gives the new id (1 ups in each round)
            else:
                current_sequence.append(stripped)                        # when line is not a header then the sequence is added to the list

    if current_id is not None:                                           # is needed for the last header since the last header would be ignored in the other if function 
        protein_map[current_id] = ''.join(current_sequence)              # test ob letzte zeile drin ist -> vllt auch für tests die funktionieren müssen damit keine 

    return protein_map


