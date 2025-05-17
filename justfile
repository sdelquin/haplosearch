# Launch Django development server
runserver:
    uv run ./manage.py runserver

# Make migrations for single app or whole project
makemigrations app="":
    uv run ./manage.py makemigrations {{ app }}

# Commit migrations for single app or whole project
migrate app="":
    uv run ./manage.py migrate {{ app }}

alias mm := mmigrate

# Make migrations & Commit migrations (all in one)
mmigrate app="":
    uv run ./manage.py makemigrations {{ app }} && uv run ./manage.py migrate {{ app }}

# Show migrations for single app or whole project
showmigrations app="":
    uv run ./manage.py showmigrations {{ app }}

# Create superuser
su:
    uv run ./manage.py createsuperuser

# Sync uv
[macos]
sync:
    uv sync --no-group prod

# Sync uv
[linux]
sync:
    uv sync --no-dev --group prod

# Deploy project to production
deploy:
    #!/usr/bin/env bash
    git pull
    just sync
    uv run ./manage.py migrate
    npm install
    uv run ./manage.py collectstatic --no-input
    supervisorctl restart haplosearch

# Open a Django shell
@sh:
    uv run ./manage.py shell

# Test project
test:
    uv run pytest
