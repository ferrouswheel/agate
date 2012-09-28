MAKE=make

PRGMS=\
	repeat-finder \
	ppt-mismatch \
	expected
	
MAKE_SRC=$(MAKE) -C src; \
	for i in $(PRGMS); do \
		mv src/$$i$$EXT .; \
	done

default: binary

binary:
	EXT=
	$(MAKE_SRC)
	
srcpkg: clean
	CWD=pwd
	cd .. && tar --wildcards -X repeatfinder/pkg/exclude -cvz -f repeatfinder/repeatfinder.tar.gz repeatfinder/*
	cd $$CWD
	mv repeatfinder.tar.gz pkg/.

dist: x86 x64
	
	
	
x86:
	EXT=
	OLD_CFLAGS="$(CFLAGS)"
	CFLAGS="$(CFLAGS) -m32"
	$$MAKE_SRC
	CFLAGS=$(OLD_CFLAGS)
	
	if [ ! -d dist.x32 ]; then \
		mkdir dist.x32; \
	fi
	for i in $(PRGMS); do \
		mv $$i dist.x32/$$i; \
	done
	
x64:
	EXT=.64
	OLD_CFLAGS=$$CFLAGS
	CFLAGS=$(CFLAGS) -m64
	$$MAKE_SRC
	CFLAGS=$(OLD_CFLAGS)
	
	if [ ! -d dist.x64 ]; then \
		mkdir dist.x64; \
	fi
	for i in $(PRGMS); do \
		mv $$i.64 dist.x64/$$i; \
	done

clean:
	$(MAKE) -C src clean
	for i in $(PRGMS); do \
		rm -f $$i; \
	done
	rm -f pkg/repeatfinder.tar.gz
