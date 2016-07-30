include ${capdMake}/make/common_makefile.mkf

COMMON_LDFLAGS=-version-info ${CAPD_VERSION_INFO}
COMMON_LIBS=

if LINK_WITH_LIBSTDCPP
   COMMON_LDFLAGS+= -no-undefined
   COMMON_LIBS+= -lstdc++
endif


pkgconfigdir = @PKG_CONFIG_DIR@

bin_SCRIPTS=
pkgconfig_DATA=
lib_LTLIBRARIES=


if HAVE_CAPD

noinst_LTLIBRARIES=libcapdweaksymbols.la
libcapdweaksymbols_la_SOURCES=main_entry.cpp
libcapdweaksymbols_la_CXXFLAGS= ${AM_CXXFLAGS} ${capd_CFLAGS} -I${top_srcdir}/../../capdAlg/include -I${top_srcdir}/../../capdAux/include
libcapdweaksymbols_la_LDFLAGS=

bin_SCRIPTS+= capd-config
pkgconfig_DATA+= capd.pc
lib_LTLIBRARIES+=libcapd.la

libcapd_la_SOURCES=
libcapd_la_LIBADD=${capd_LIBS} ${COMMON_LIBS}
libcapd_la_LDFLAGS=${COMMON_LDFLAGS} -weak_library libcapdweaksymbols.la
libcapd_la_DEPENDENCIES=libcapdweaksymbols.la ${filter %.a %.la,${libcapd_la_LIBADD}}


# Below is hard-coded rule for cross building with mxe. There we do not have shared libraries and libtool do not want to build dll's.
capd_dll: libcapd.la  libcapdweaksymbols.la
	$(CXX) -shared -nostdlib /opt/mxe/mxe/usr/lib/gcc/$(host_alias)/4.8.1/../../../../$(host_alias)/lib/dllcrt2.o /opt/mxe/mxe/usr/lib/gcc/$(host_alias)/4.8.1/crtbegin.o -L. -Wl,--allow-multiple-definition  -Wl,--whole-archive .libs/libcapd.a  -Wl,-no--whole-archive -o .libs/libcapd.dll -L/opt/mxe/mxe/usr/lib/gcc/$(host_alias)/4.8.1 -L/opt/mxe/mxe/usr/lib/gcc/$(host_alias)/4.8.1/../../../../$(host_alias)/lib -lstdc++ -lmingw32 -lgcc -lmoldname -lmingwex -lmsvcrt -ladvapi32 -lshell32 -luser32 -lkernel32 -lmingw32 -lgcc -lmoldname -lmingwex -lmsvcrt /opt/mxe/mxe/usr/lib/gcc/$(host_alias)/4.8.1/crtend.o  -Wl,--enable-auto-image-base -Xlinker --out-implib -Xlinker libcapd.dll.a -lboost_regex-mt -lboost_system-mt -lboost_filesystem-mt -llog4cxx -lapr-1 -laprutil-1 -lapr-1 -lexpat -liconv -lstdc++  -lwsock32 -lws2_32 -lrpcrt4 .libs/libcapdweaksymbols.a


endif

if HAVE_CAPD_GUI

bin_SCRIPTS+= capd-gui-config
pkgconfig_DATA+= capd-gui.pc
lib_LTLIBRARIES+=libcapd-gui.la

libcapd_gui_la_SOURCES=
libcapd_gui_la_LIBADD=${capd_gui_LIBS} ${COMMON_LIBS}
libcapd_gui_la_LDFLAGS= ${COMMON_LDFLAGS} -weak_library libcapdweaksymbols.la
libcapd_gui_la_DEPENDENCIES=libcapdweaksymbols.la ${filter %.a %.la,${libcapd_gui_la_LIBADD}}

endif

if HAVE_MPCAPD

bin_SCRIPTS+= mpcapd-config
pkgconfig_DATA+= mpcapd.pc
lib_LTLIBRARIES+=libmpcapd.la

libmpcapd_la_SOURCES=
libmpcapd_la_LIBADD=${mpcapd_LIBS} ${COMMON_LIBS}
libmpcapd_la_LDFLAGS=${COMMON_LDFLAGS} -weak_library libcapdweaksymbols.la
libmpcapd_la_DEPENDENCIES=libcapdweaksymbols.la ${filter %.a %.la,${libmpcapd_la_LIBADD}}

endif

if HAVE_MPCAPD_GUI

bin_SCRIPTS+= mpcapd-gui-config
pkgconfig_DATA+= mpcapd-gui.pc
lib_LTLIBRARIES+=libmpcapd-gui.la

libmpcapd_gui_la_SOURCES=
libmpcapd_gui_la_LIBADD=${mpcapd_gui_LIBS} ${COMMON_LIBS}
libmpcapd_gui_la_LDFLAGS=${COMMON_LDFLAGS}  -weak_library libcapdweaksymbols.la
libmpcapd_gui_la_DEPENDENCIES=libcapdweaksymbols.la ${filter %.a %.la,${libmpcapd_gui_la_LIBADD}}

endif

all: chmodScripts

chmodScripts:
	chmod +x ${top_capd_build}/bin/*

.PHONY: chmodScripts

.PHONY: lib

lib: chmodScripts
	    
