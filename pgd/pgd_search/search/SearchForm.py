import math
from django import forms

from pgd_core.models import Protein
from pgd_search.models import searchSettings
from pgd_search.views import RESIDUE_INDEXES
from constants import AA_CHOICES, SS_CHOICES

"""
Search form used by search handler.
"""

class SearchFormBase(forms.Form):
    threshold       = forms.ChoiceField(choices=[(25,25),(90,90)], required=False)
    resolutionMin   = forms.FloatField(required=False, min_value=0, initial=0, widget=forms.TextInput(attrs={'size':3}))
    resolutionMax   = forms.FloatField(required=False, min_value=0, initial=1.5, widget=forms.TextInput(attrs={'size':3}))
    proteins        = forms.ModelMultipleChoiceField(queryset=Protein.objects.all().order_by('code'), required=False)
    residues        = forms.ChoiceField(choices=[(i,i) for i in range(1, searchSettings.segmentSize+1)], initial=5)

# Build a dict for the fields of variable number
form_dict = {'__module__' : 'pgd_search.views'}

for i in RESIDUE_INDEXES:
    form_dict["aa_%i" % i]      = forms.MultipleChoiceField(choices=AA_CHOICES, required=False, widget=forms.SelectMultiple(attrs={'class':'field'}))
    form_dict["aa_i_%i" % i]    = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))

    FORM_SS_CHOICES = [('','')]
    for choice in SS_CHOICES:
        FORM_SS_CHOICES.append(choice)
    form_dict["ss_%i" % i]      = forms.ChoiceField(choices=FORM_SS_CHOICES, required=False, widget=forms.Select(attrs={'class':'field'}))
    form_dict["ss_i_%i" % i]    = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))

    # the loops here are just to save on space/typing
    for j in range(1,8):
        form_dict["a%i_%i" % (j,i)]     = forms.ChoiceField(choices=AA_CHOICES, required=False, widget=forms.TextInput(attrs={'class':'field', 'size':8}))
        form_dict["a%i_i_%i" % (j,i)]   = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))
    for j in range(1,6):
        form_dict["L%i_%i" % (j, i)]    = forms.FloatField(required=False, widget=forms.TextInput(attrs={'class':'field', 'size':8}))
        form_dict["L%i_i_%i" % (j, i)]  = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))
    for j in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta"):
        form_dict["%s_%i" % (j, i)]     = forms.FloatField(required=False, widget=forms.TextInput(attrs={'class':'field', 'size':8}))
        form_dict["%s_i_%i" % (j, i)]   = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))

# Create the Search Form with the fields from the dict
SearchForm = type('SearchForm', (SearchFormBase,), form_dict)