{
  inputs = {
    naersk.url = "github:nix-community/naersk/master";
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils, naersk }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        naersk-lib = pkgs.callPackage naersk { };
      in {
        defaultPackage = naersk-lib.buildPackage {
          nativeBuildInputs = with pkgs; [
            futhark
            jemalloc
            rust-jemalloc-sys
          ];
          src = ./.;
        };
        devShell = with pkgs;
          mkShell {
            buildInputs = [
              cargo
              rustc
              rustfmt
              rustPackages.clippy
              rust-analyzer
              futhark
              jemalloc
              rust-jemalloc-sys
            ];
            RUST_SRC_PATH = rustPlatform.rustLibSrc;
          };
      });
}
