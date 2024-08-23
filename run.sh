#!/bin/bash

source .venv/bin/activate
cd $(dirname $0)
exec gunicorn -c gunicorn.conf.py main.wsgi:application
