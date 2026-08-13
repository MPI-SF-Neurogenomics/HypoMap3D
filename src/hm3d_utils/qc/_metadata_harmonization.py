# import re

# def fix_obs_names(s):
#     # Case 1 and 2 are mutually exclusive!
#     intermediate = None

#     # Case 1: '../data/v2-04-SRR26100757_SC2-04' → SRR26100757_SC2-04. each cell has a different SRR in smartseq, so this is the cell barcode
#     match1 = re.match(r".*/v\d+-\d+-([^-_]+)_([^/]+)$", s)
#     if match1:
#         intermediate = f"{match1.group(1)}_{match1.group(2)}"

#     # Case 2: 'AAACCCACAACCAGAG-1_SC1-01-SAMN18754353' → 'AAACCCACAACCAGAG_SC1-01-SAMN18754353'
#     match2 = re.match(r"^([ATGC]+)-1_(.+)$", s)
#     if match2:
#         intermediate = f"{match2.group(1)}_{match2.group(2)}"

#     if intermediate:
#         s = intermediate

#     # Harmonize dataset name format without dropping prefix
#     match_sc = re.search(r'(SC\d+)[-_](\d+)[-_](\S+)', s)
#     if match_sc:
#         sc_fixed = f"{match_sc.group(1)}-{match_sc.group(2)}_{match_sc.group(3)}"
#         s = s[:match_sc.start()] + sc_fixed

#     return s



# def harmonize_sample_names(value):
#     # Match SC followed by digits, separator (either - or _), digits, then anything else
#     match = re.match(r'SC(\d+)[-_](\d+)[-_](.+)', value)
#     if match:
#         return f'SC{match.group(1)}-{match.group(2)}_{match.group(3)}'
#     return value  # Return original if it doesn't match
