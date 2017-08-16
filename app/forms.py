from django import forms

OPERATION_CHOICES = (
    ("S2H", "Sequences -> Haplotypes"),
    ("H2S", "Haplotypes -> Sequences"),
)

NOMENCLATURE_CHOICES = (
    ("POP", "Population genetics"),
    ("FOR", "Forensic genetics"),
)


class HaploSearchForm(forms.Form):
    inputfile = forms.FileField()
    operation = forms.ChoiceField(choices=OPERATION_CHOICES)
    nomenclature = forms.ChoiceField(choices=NOMENCLATURE_CHOICES)
    hvri = forms.BooleanField(
        label="I am only analysing HVRI",
        initial=False,
        required=False,
        help_text="HVRI haplotypes with three-digits mutation"
    )
