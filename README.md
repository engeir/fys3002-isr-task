# ISR Exercise

> Solution to exercise in FYS-3002.

## Install

Clone with

```sh
git clone https://github.com/engeir/fys3002-isr-task.git
```

Dependencies are handled using poetry:

<details><summary><i><b>Install</b></i></summary><br><ul>

```sh
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
```

</ul></details>

Install and create a virtualenv using poetry from folder `isr-task`:

```sh
poetry install
poetry shell
python -m isr-task
```
