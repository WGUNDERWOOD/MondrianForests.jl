{
  description = "MondrianForests.jl";
  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-23.11;
    nixpkgs-unstable.url = github:NixOS/nixpkgs/nixos-unstable;
    flake-utils.url = github:numtide/flake-utils;
  };
  outputs = {
    self,
    nixpkgs,
    nixpkgs-unstable,
    flake-utils,
  }:
    with flake-utils.lib;
      eachSystem allSystems (system: let
        pkgs = nixpkgs.legacyPackages.${system};
        pkgs-unstable = nixpkgs-unstable.legacyPackages.${system};
        julia = pkgs-unstable.julia.withPackages [
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
