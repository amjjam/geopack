# Build will compile objects
# install will install the .so file inside ./lib
# uninstall will remove .so file from ./lib
# clean will remove .so file and .o file from ./src

build:
	$(MAKE) -C src build
install:
	$(MAKE) -C src install
clean:
	$(MAKE) -C src clean
uninstall:
	$(MAKE) -C src uninstall
