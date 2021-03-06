from django.conf.urls import url
from . import views

urlpatterns = [
    url(r"^$", views.index, name="index"),
    url(r"^start/$", views.start, name="start"),
    url(r"^download/$", views.download, name="download"),
    url(r"^help/$", views.help, name="help"),
]
