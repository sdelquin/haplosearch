from fabric.api import env, local, cd, run

env.hosts = ['cloud']


def deploy():
    local('git push')
    with cd('~/haplosearch'):
        run('git pull')
        run('pipenv install')
        with cd('app/static'):
            run('npm install')
        run('pipenv run python manage.py collectstatic --noinput')
        run('supervisorctl restart haplosearch')
