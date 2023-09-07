PREAMBLE="\usepackage[sfdefault,light]{FiraSans}\DeclareMathSizes{10}{11}{9}{9}"

all: svg png

.PHONY: all svg png

svg:
	pgf2svg -p $(PREAMBLE) weather_data.pgf

png:
	inkscape -w 2000 weather_data.svg -o weather_data.png
