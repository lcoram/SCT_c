#SCT_rewrite.o: SCT_rewrite.c
#	gcc $< -fPIC -lgslcblas -lgsl -lblas -L/usr/lib -lm -o $@

SCT_wrapper.so: SCT_wrapper.c
	gcc $< -fPIC -lgslcblas -lgsl -lblas -L/usr/lib -lm -shared -o $@
