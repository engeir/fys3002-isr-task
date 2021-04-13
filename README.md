# ISR Exercise

> Solution to exercise in FYS-3002.

## Install

Clone with

```sh
git clone https://github.com/engeir/fys3002-isr-task.git isr-task
```

Dependencies are handled using poetry:

<details><summary><i><b>Install poetry</b></i></summary><br><ul>

```sh
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
```

</ul></details>

Install dependencies using poetry from folder `isr-task`:

```sh
poetry install
```

Run as a package or as single scripts:

```sh
poetry run isr-task
poetry run python soln.py
poetry run python timing.py
.
.
.
```
