PYTHONPATH := $(CURDIR)/c_python_ext/bin

.PHONY: test

test:
	PYTHONPATH=$(PYTHONPATH) python tests/tests.py
