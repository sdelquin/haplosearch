from fabric.api import env, local, prefix, cd, run

env.hosts = ["haplosearch.com"]


def deploy():
    local("git push")
    with prefix("source ~/.virtualenvs/haplosearch/bin/activate"):
        with cd("~/haplosearch"):
            run("git pull")
            run("pip install -r requirements.txt")
            run("bower install")
            run("python manage.py collectstatic --noinput")
            run("supctl restart haplosearch")
