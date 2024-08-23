from fabric.api import cd, env, local, prefix, run

env.hosts = ['haplosearch.com']


def deploy():
    local('git push')
    with prefix('source .venv/bin/activate'):
        with cd('~/code/haplosearch'):
            run('git pull')
            run('pip install -r requirements.txt')
            with cd('app/static'):
                run('npm install')
            run('python manage.py collectstatic --noinput')
            run('supervisorctl restart haplosearch')
