#!/bin/bash

source ~/.virtualenvs/haplosearch/bin/activate
cd $(dirname $0)
exec gunicorn -c gunicorn.conf.py main.wsgi:application
