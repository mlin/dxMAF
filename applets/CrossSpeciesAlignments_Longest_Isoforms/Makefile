all: default

default:
	$(MAKE) -C src all

clean:
	$(MAKE) -C src clean

test: all
	test/test.sh

.PHONY: all default clean test
