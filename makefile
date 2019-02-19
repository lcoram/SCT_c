sct_smart_boxes.so: sct_smart_boxes.c
	gcc $< -fPIC -lgslcblas -lgsl -lblas -L/usr/lib -lm -shared -o $@

sct_test: sct_test.c sct_smart_boxes.c makefile
	gcc sct_test.c -g -pg -lgslcblas -lgsl -lblas -L/usr/lib -lm -o $@
