PYTHON ?= python
ENSEMBLE ?= outputs/v3/demo/ensemble
FIGURES ?= outputs/v3/demo/figures

.PHONY: test figures clock demo
test:
	PYTHONPATH=src $(PYTHON) -m pytest tests/test_v2_*.py tests/test_v3_*.py -q

figures:
	PYTHONPATH=src $(PYTHON) analysis/figures.py --ensemble "$(ENSEMBLE)" --output "$(FIGURES)"

clock:
	PYTHONPATH=src $(PYTHON) tools/validate_clock_v3.py --output outputs/clock-validation

demo:
	PYTHONPATH=src $(PYTHON) examples/v3/make_demo.py --output outputs/v3/demo
