# ISR Exercise

> Solution to exercise in FYS-3002.

## Install

Clone with

```sh
git clone https://github.com/engeir/fys3002-isr-task.git isr-task
```

Built using `pytohn3.9`.


<details><summary><i><b>Install python3.9</b></i></summary><br><ul>

This assumes that [pyenv](https://github.com/pyenv/pyenv-installer) is installed.

```sh
pyenv install 3.9.2
```

</ul></details>

It is good practice to use virtual environments. Creating one can be done in different
ways, but with `pyenv` installed this is as easy as

```sh
pyenv virtualenv 3.9.2 my-new-venv
pyenv activate my-new-venv
```

In the above code snippet, the virtual environment uses python3.9.2, which has to be
installed before you run the commands. The name of the environment is `my-new-venv`, which
is activated on the second line.

<details><summary><i><b>Install poetry</b></i></summary><br><ul>

```sh
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
```

</ul></details>

Install dependencies using poetry from folder `isr-task`:

```sh
poetry install
```

If a virtual environment is not activated when the `poetry install` command is run, a new
one is created by poetry inside the project folder.

## Run

Run as a package or as single scripts:

```sh
poetry run isr-task
poetry run python src/isr_task/soln.py
poetry run python src/isr_task/timing.py
.
.
.
```

## About

Uses [`black`](https://black.readthedocs.io/en/stable/),
[`isort`](https://pycqa.github.io/isort/) and [`flake8`](https://flake8.pycqa.org/en/latest/) for formatting at commit
with [`pre-commit`](https://pre-commit.com/).
