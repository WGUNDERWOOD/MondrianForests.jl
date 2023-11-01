all:
	@rm -f src/*.cov
	@julia --project --threads 6 make.jl
