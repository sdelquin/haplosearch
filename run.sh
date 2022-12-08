#!/bin/bash

source ~/.pyenv/versions/haplosearch/bin/activate
cd $(dirname $0)
exec gunicorn -c gunicorn.conf.py main.wsgi:application
