INSTALL_DIR = /prod/kweh/bin/
BIN_DIR = /prod/bin

.PHONY: install clean fullclean

dist/kweh: kweh.py output.py read_mesh.py data_structures.py
	pyinstaller --onefile kweh.py

install: dist/kweh
	cp dist/kweh $(INSTALL_DIR)
	chmod +x $(INSTALL_DIR)/kweh
	ln -sf $(INSTALL_DIR)/kweh $(BIN_DIR)/kweh


clean:
	-rm kweh.spec
	-rm -r build
	-rm -r __pycache__

fullclean:
	-rm kweh.spec
	-rm -r build
	-rm -r __pycache__
	-rm -r dist
