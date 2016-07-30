include ${capdMake}/make/common_makefile.mkf

SUBDIRS=include src tests

if WITH_CAPD_EXAMPLES
  SUBDIRS+= examples
endif

.PHONY: lib

lib:
	$(MAKE) $(AM_MAKEFLAGS) -C src

EXTRA_DIST = docs

dist-hook:
	rm -rf `find $(distdir) -type d -name .svn`
	rm -rf `find $(distdir) -type d -name .git`
