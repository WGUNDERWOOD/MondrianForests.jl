let
pkgs = import <nixpkgs> { };
in pkgs.mkShell {
    buildInputs = with pkgs; [
        julia
        cacert
        python3
        python3Packages.matplotlib
        texlive.combined.scheme-full
        git
    ];
    shellHook = ''
        # run this to link Julia PyCall package to nixpkgs python3
        #julia --color=yes -e 'using Pkg; ENV["PYTHON"]="${pkgs.python3}/bin/python3"; Pkg.build("PyCall")'
    '';
}
