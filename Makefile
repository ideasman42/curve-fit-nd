PYTHONPATH := $(CURDIR)/c_python_ext/build/bin

.PHONY: test

test:
	PYTHONPATH=$(PYTHONPATH) python tests/tests.py
