from fabric.api import env, local, cd, run, prefix

env.hosts = ['cloud']


def deploy():
    local('git push')
    with prefix('source ~/.virtualenvs/haplosearch/bin/activate'):
        with cd('~/code/haplosearch'):
            run('git pull')
            run('pip install -r requirements.txt')
            with cd('app/static'):
                run('npm install')
            run('python manage.py collectstatic --noinput')
            run('supervisorctl restart haplosearch')
