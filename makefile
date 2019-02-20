LIBS=-lgslcblas -lgsl -lblas -L/usr/lib -lm  -L/usr/lib
CFLAGS=-g -pg
sct_smart_boxes.o: sct_smart_boxes.c
	gcc $< -c $(CFLAGS) $(LIBS) -o $@

sct_smart_boxes.so: sct_smart_boxes.c
	gcc $< -fPIC $(LIBS) -shared -O3 -o $@

sct_run: sct_run.c sct_smart_boxes.o makefile
	gcc sct_run.c sct_smart_boxes.o $(CFLAGS) $(LIBS) -o $@

test: test.c sct_smart_boxes.o
	gcc test.c sct_smart_boxes.o $(CFLAGS) $(LIBS) -o $@
