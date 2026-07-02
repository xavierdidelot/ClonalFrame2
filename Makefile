.PHONY: all build debug clean

all: build

build:
	cmake --preset=release
	$(MAKE) -C build-release
	cp build-release/ClonalFrame2 .

debug:
	cmake --preset=debug
	$(MAKE) -C build-debug

clean:
	rm -rf build-release build-debug
