EXES = jb-attribute jb-columns jb-counter jb-delete jb-derive jb-global jb-header jb-merge jb-print jb-promote jb-region jb-rename

all: $(EXES)

jb-attribute: jb-attribute.o jb.o
	$(CXX) -o jb-attribute jb-attribute.o jb.o
jb-columns: jb-columns.o jb.o
	$(CXX) -o jb-columns jb-columns.o jb.o
jb-counter: jb-counter.o jb.o
	$(CXX) -o jb-counter jb-counter.o jb.o
jb-delete: jb-delete.o jb.o
	$(CXX) -o jb-delete jb-delete.o jb.o
jb-derive: jb-derive.o jb.o sort.o
	$(CXX) -o jb-derive jb-derive.o jb.o sort.o
jb-global: jb-global.o jb.o
	$(CXX) -o jb-global jb-global.o jb.o
jb-header: jb-header.o jb.o
	$(CXX) -o jb-header jb-header.o jb.o
jb-merge: jb-merge.o jb.o
	$(CXX) -o jb-merge jb-merge.o jb.o
jb-print: jb-print.o jb.o
	$(CXX) -o jb-print jb-print.o jb.o
jb-promote: jb-promote.o jb.o
	$(CXX) -o jb-promote jb-promote.o jb.o
jb-region: jb-region.o jb.o
	$(CXX) -o jb-region jb-region.o jb.o
jb-rename: jb-rename.o jb.o
	$(CXX) -o jb-rename jb-rename.o jb.o

install: $(EXES)
	-@cp $(EXES) ../../bin
clean:
	-@rm -f $(EXES)
	-@rm -f *.o
