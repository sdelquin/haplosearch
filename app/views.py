from django.shortcuts import render
from .forms import HaploSearchForm


def index(request):
    return render(request, "app/index.html")


def start(request):
    form = HaploSearchForm()
    return render(request, "app/start.html", {"form": form})
