#!/bin/bash
# Master script.

cd "$(dirname "$0")"
source ~/.virtualenvs/haplosearch/bin/activate
exec uwsgi --ini uwsgi.cfg
