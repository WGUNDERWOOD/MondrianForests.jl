all:
	@rm src/*.cov
	@julia --project --threads 6 make.jl
