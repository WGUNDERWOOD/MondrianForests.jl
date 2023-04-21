all: format

.PHONY: format

format:
	julia -e 'using JuliaFormatter; format(".")'
