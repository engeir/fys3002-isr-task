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

<details><summary><i><b>Install poetry</b></i></summary><br><ul>

```sh
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
```

</ul></details>

Install dependencies using poetry from folder `isr-task`:

```sh
poetry install
```

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
