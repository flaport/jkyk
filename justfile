dist:
	uv run python -m build --wheel

dev:
	uv sync --all-extras
	uv pip install -e .
	uv run pre-commit install

uv:
	curl -LsSf https://astral.sh/uv/install.sh | sh

inits:
	cd python/jkyk && uv run mkinit --relative --recursive --write && uv run ruff format __init__.py

ipykernel:
	uv run python -m ipykernel install --user --name jkyk --display-name jkyk

test:
	uv run pytest -s -n logical

docs:
	uv run mkdocs build

serve:
	uv run mkdocs serve -a localhost:8080

nbrun:
	find nbs -maxdepth 1 -mindepth 1 -name "*.ipynb" -not -path "*/.ipynb_checkpoints/*" -not -path "./.venv/*" | xargs parallel -j `nproc --all` uv run papermill {} {} -k jkyk :::

nbdocs:
	find nbs -maxdepth 1 -mindepth 1 -name "*.ipynb" -not -path "*/.ipynb_checkpoints/*" -not -path "./.venv/*" | xargs parallel -j `nproc --all` uv run jupyter nbconvert --to markdown --embed-images {} --output-dir docs/nbs ':::'

tree:
	@tree -a -I .git --gitignore

clean:
	rm -rf site
	rm -rf .venv
	rm -rf docs/nbs/*
	find python -name "*.c" | xargs rm -rf
	find python -name "*.pyc" | xargs rm -rf
	find python -name "*.so" | xargs rm -rf
	find python -name "*.pyd" | xargs rm -rf
	find . -name "*.egg_info" | xargs rm -rf
	find . -name ".ipynb_checkpoints" | xargs rm -rf
	find . -name ".mypy_cache" | xargs rm -rf
	find . -name ".pytest_cache" | xargs rm -rf
	find . -name ".ruff_cache" | xargs rm -rf
	find . -name "__pycache__" | xargs rm -rf
	find . -name "build" | xargs rm -rf
	find . -name "builds" | xargs rm -rf
	find . -name "dist" -not -path "*node_modules*" | xargs rm -rf
