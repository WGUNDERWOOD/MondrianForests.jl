{
  description = "MondrianForests.jl";
  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-24.05;
    flake-utils.url = github:numtide/flake-utils;
  };
  outputs = {
    self,
    nixpkgs,
    flake-utils,
  }:
    with flake-utils.lib;
      eachSystem allSystems (system: let
        pkgs = nixpkgs.legacyPackages.${system};
        julia = pkgs.julia_19.withPackages [
          "Aqua"
          "Colors"
          "Coverage"
          "CSV"
          "DataFrames"
          "Dates"
          "Distributions"
          "Documenter"
          "JuliaFormatter"
          "Plots"
          "PyPlot"
          "Random"
          "Suppressor"
          "Test"
        ];
      in {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            julia
          ];
        };
      });
}
