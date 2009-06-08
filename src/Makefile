build:
	$(MAKE) -C enzohop
	$(MAKE) -C inits
	$(MAKE) -C ring
	$(MAKE) -C enzo
clean:
	$(MAKE) -C enzohop clean
	$(MAKE) -C inits clean
	$(MAKE) -C ring  clean
	$(MAKE) -C enzo  clean
install: build
	$(MAKE) -C enzohop install
	$(MAKE) -C inits install
	$(MAKE) -C ring  install
	$(MAKE) -C enzo  install


