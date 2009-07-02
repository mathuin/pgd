from django.db import models
from constants import AA_CHOICES, SS_CHOICES, AA_CHOICES_DICT

# Protein model
# (was 'protein_info')
# Contains the general information about a protein
class Protein(models.Model):
    code        = models.CharField(max_length=4, primary_key=True, unique=True)
    threshold   = models.IntegerField() # new type; this should probably be a boolean type
    resolution  = models.FloatField(db_index=True)
    rfactor     = models.FloatField()

    def __unicode__(self):
        return self.code

# Chain model
# Contains information about chains within a protein
# 
# NOTE: This model is in a temporary state that allows easier conversion from the current database
# the primary key will be a 5 character string derived from the chainID and protein code.
# this will all a very simple query to generate both the table of chains and set the FK
# in the residue table.  Once splicer is updated and importing directly to the database the PK
# will be changed to an integer.
class Chain (models.Model):
        id      = models.CharField(max_length=5, primary_key=True)
        protein = models.ForeignKey(Protein, related_name='chains')
        code    = models.CharField(max_length=1)

# Residue model
# (was 'protein')
# Contains information regarding each amino acid and its geometry
# (Note: fields need to be commented)
class Residue(models.Model):
    protein         = models.ForeignKey(Protein, related_name='residues')
    chain           = models.ForeignKey(Chain, related_name='residues')
    aa              = models.CharField(max_length=1, choices=AA_CHOICES) # new type
    chainID         = models.CharField(max_length=1)
    # Is oldID necessary? The correct type?
    oldID           = models.CharField(max_length=5, null=True)
    chainIndex      = models.PositiveIntegerField()
    a1              = models.FloatField(null=True)
    a2              = models.FloatField()
    a3              = models.FloatField()
    a4              = models.FloatField()
    a5              = models.FloatField()
    a6              = models.FloatField(null=True)
    a7              = models.FloatField(null=True)
    L1              = models.FloatField(null=True)
    L2              = models.FloatField()
    L3              = models.FloatField()
    L4              = models.FloatField()
    L5              = models.FloatField()
    ss              = models.CharField(max_length=1, choices=SS_CHOICES) # new type (was blob, but all entries 1 char)
    phi             = models.FloatField()
    psi             = models.FloatField()
    ome             = models.FloatField()
    chi             = models.FloatField()
    bm              = models.FloatField()
    bs              = models.FloatField()
    bg              = models.FloatField(null=True)
    h_bond_energy   = models.FloatField()
    zeta            = models.FloatField()
    terminal_flag   = models.BooleanField() 
    xpr             = models.BooleanField() # this field may not be necessary; it has never been implemented

    #def __str__(self):
    #    return '%d' % self.chainIndex

    def __getattribute__(self,name):
        if name == 'aa_full':
            return AA_CHOICES_DICT[self.aa]
    
        # normal attribute
        else:
            return object.__getattribute__(self, name)
