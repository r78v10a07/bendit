SHELL := /bin/bash
.PHONY : test all

all:	bendit test
	echo OK!!!!!!!!!

bendit:
	(mkdir -p bin && cd src && $(MAKE) all)

test:
	echo "Testing bendIt"
	./bin/bendIt -s test/t.fasta -o test/s
	diff test/text_out test/s > /dev/null || (echo "BendIt Test failed" >&2 && exit 1)
	rm test/s*
	echo "BendIt Test OK"

clean:
	rm -vfr bin test/s*
	(cd src && $(MAKE) clean )
