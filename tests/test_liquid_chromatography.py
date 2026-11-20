import proteosim_max as ps

def test_predict_lc_retention_times():

    peptides = ['PEPTIDE']
    expected = {'PEPTIDE': 7.8}

    actual = ps.predict_lc_retention_times(peptides)
    assert actual == expected  # mit py test there can be annäherungen kommen 

test_predict_lc_retention_times()
#circular !!! muss ich ändern os bringt der test nichts 

# hier eher testen ob die struktur der antwort richtig is -> also ob es aussieht wie das expected dict 



def test_select_retention_time_window():
    peptide_rt_map = {'A': 5, 'B': 10, 'C': 15}
    selected = ps.select_retention_time_window(peptide_rt_map, lower_ret_time=8, upper_ret_time=12)

    assert selected == ["B"]

test_select_retention_time_window()
# better make two ddifferent test fordifferent situation and also for edge/corner cases 