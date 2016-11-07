from Bio.Data.IUPACData import protein_letters_1to3
import operator

# The possible values for the 'aa' field (in Protein and elsewhere)
AA_CHOICES = sorted(protein_letters_1to3.items(), key=operator.itemgetter(1))

# The possible values for the 'ss' field (in Protein and elsewhere)
'''SS_CHOICES = (
    ('G', 'G'),
    ('g', 'g'),
    ('H', 'H'),
    ('h', 'h'),
    ('I', 'I'),
    ('i', 'i'),
    ('T', 'T'),
    ('t', 't'),
    ('E', 'E'),
    ('e', 'e'),
    ('B', 'B'),
    ('S', 'S'),
)'''
SS_CHOICES = (
    ('H', '&alpha; helix'),
    #('h', '&alpha; helix (terminal)'),
    ('G', '3<sub>10</sub> helix'),
    #('g', '3<sub>10</sub> helix (terminal)'),
    ('E', '&beta; sheet'),
    #('e', '&beta; sheet (terminal)'),
    ('T', 'Turn'),
    #('t', 'H-bonded turn (terminal)'),
    ('S', 'Bend'),
    ('B', '&beta;-bridge'),
    ('I', '&pi; helix'),
    #('i', '&pi; helix'),
)

PLOT_PROPERTY_CHOICES = [
    ("phi","phi"),
    ("psi","psi"),    
    ("L1","L1"),
    ("L2","L2"),
    ("L3","L3"),
    ("L4","L4"),
    ("L5","L5"),
    ("a1","a1"),
    ("a2","a2"),
    ("a3","a3"),
    ("a4","a4"),
    ("a5","a5"),
    ("a6","a6"),
    ("a7","a7"),
    ("ome","ome"),
    ("chi1","chi1"),
    ("chi2","chi3"),
    ("chi3","chi3"),
    ("chi4","chi4"),
    ("chi5","chi5"),
    ("zeta","zeta")
    #("h_bond_energy","H Bond"),
    ]
