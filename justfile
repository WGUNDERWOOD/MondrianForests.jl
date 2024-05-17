all:
  @rm -f src/*.cov
  @julia --project --color=yes --threads 6 make.jl
