from fabric.api import env, local, cd, run

env.hosts = ['haplosearch.com']


def deploy():
    local('git push')
    with cd('~/haplosearch'):
        run('git pull')
        run('pipenv install')
        run('bower install')
        run('pipenv run python manage.py collectstatic --noinput')
        run('supctl restart haplosearch')
