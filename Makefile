PREFIX    := /usr/local
BINDIR    := $(PREFIX)/bin

.PHONY: all build debug clean doc

all: build

build:
	cmake --preset=release
	$(MAKE) -C build-release
	cp build-release/ClonalFrame2 .

debug:
	cmake --preset=debug
	$(MAKE) -C build-debug

clean:
	rm -rf build-release build-debug docs cmake-build-debug generated

doc:
	doxygen

install:
	install ClonalFrame2 $(BINDIR)/ClonalFrame2

uninstall:
	rm -f $(BINDIR)/ClonalFrame2

