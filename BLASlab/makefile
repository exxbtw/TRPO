CC = gcc
CFLAGS = -I./OpenBLAS -L./OpenBLAS
LIBS = -lopenblas

all: run_tests

build:
	$(CC) tests.c $(CFLAGS) $(LIBS) -o tests

run_tests: build
	@echo "Running tests with 1 thread..."
	export OPENBLAS_NUM_THREADS=1 && ./tests
	@echo "Running tests with 4 threads..."
	export OPENBLAS_NUM_THREADS=4 && ./tests

clean:
	rm -f tests