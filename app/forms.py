from django import forms

OPERATION_CHOICES = (
    ("S2H", "Get haplotype from sequence"),
    ("H2S", "Get sequence from haplotype"),
)

NOMENCLATURE_CHOICES = (
    ("POP", "Population genetics"),
    ("FOR", "Forensic genetics"),
)


class HaploSearchForm(forms.Form):
    inputfile = forms.FileField()
    operation = forms.ChoiceField(choices=OPERATION_CHOICES)
    nomenclature = forms.ChoiceField(choices=NOMENCLATURE_CHOICES)
