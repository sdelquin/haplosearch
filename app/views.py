from django.shortcuts import render
from django.http import HttpResponse
from tempfile import mktemp
import os
import time
from .forms import HaploSearchForm
from .utils import handle_uploaded_file
from .haploutils import manage_haplosearch
from .exceptions import HaploException


def index(request):
    return render(request, "index.html")


def start(request):
    error, msg, outputfile_path, elapsed_time = False, None, None, None
    if request.method == "POST":
        form = HaploSearchForm(request.POST, request.FILES)
        if form.is_valid():
            inputfile = form.cleaned_data["inputfile"]
            operation = form.cleaned_data["operation"]
            nomenclature = form.cleaned_data["nomenclature"]
            inputfile_path, outputfile_path = mktemp(), mktemp()
            handle_uploaded_file(inputfile, inputfile_path)
            try:
                start = time.time()
                manage_haplosearch(
                    inputfile_path, outputfile_path, nomenclature, operation
                )
                end = time.time()
                elapsed_time = (end - start)
            except (HaploException, Exception) as e:
                error = True
                msg = e
                os.remove(outputfile_path)
            finally:
                os.remove(inputfile_path)
    else:
        form = HaploSearchForm()
    return render(
        request,
        "start.html",
        {
            "form": form,
            "error": error,
            "msg": msg,
            "outputfile_path": outputfile_path,
            "elapsed_time": elapsed_time
        }
    )


def download(request):
    filepath = request.POST["filepath"]
    filecontent = open(filepath)
    os.remove(filepath)
    response = HttpResponse(filecontent, content_type="text/plain")
    response["Content-Disposition"] = "attachment; filename=output.txt"
    return response


def help(request):
    return render(request, "help.html")
