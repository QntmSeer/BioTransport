.PHONY: all test-python test-rust test clean build-rust sync-remote run-remote

all: test

build-rust:
	cd 02_continuum_transport/biotransport-rs && cargo build --release

test-rust:
	cd 02_continuum_transport/biotransport-rs && cargo test

test-python:
	pytest 01_microscale_md/tests 02_continuum_transport/python/tests -v

test: test-python

clean:
	rm -rf build/ dist/ *.egg-info .pytest_cache
	cd 02_continuum_transport/biotransport-rs && cargo clean

sync-remote:
	rsync -avz --exclude 'target' --exclude '__pycache__' --exclude '*.xtc' . agni@192.168.1.112:~/multiscale-bioparticle-transport/

run-remote:
	ssh agni@192.168.1.112 "cd ~/multiscale-bioparticle-transport && make test-rust"
