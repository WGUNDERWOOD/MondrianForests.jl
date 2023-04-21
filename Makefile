all: format

.PHONY: all format

format:
	julia -e 'using JuliaFormatter; format(".")'
