{
  description = "coalescent simulations";
  nixConfig.bash-prompt = "[coalescent]$ ";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    utils = { url = "github:numtide/flake-utils"; };
    slim-runner = { url = "github:EgorLappo/slim-runner"; };

  };

  outputs = { self, nixpkgs, utils, slim-runner }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        py = pkgs.python3;

        py-env = py.withPackages (ps:
          with ps; [
            tqdm
            numpy
            pandas
            polars
            pyarrow
            scipy
            matplotlib
            seaborn
            jupyter
            multiprocess
          ]);

        slim-runner-exe = slim-runner.packages.${system}.default;

      in rec {
        devShells.default = with pkgs;
          mkShell {
            name = "slimsim";
            buildInputs = [ py-env messer-slim slim-runner-exe steel ];

            shellHook = ''
              export PYTHONPATH="${py-env}/bin/python"
            '';
          };
      });
}
