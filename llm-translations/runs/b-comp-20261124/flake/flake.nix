{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    nixpkgs-stable.url = "github:NixOS/nixpkgs/nixos-26.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, nixpkgs-stable, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

        stable = import nixpkgs-stable {
          inherit system;
        };

        pythonWithPackages = pkgs.python3.withPackages (ps: with ps; [
          pandas
          psutil
        ]);
      in
      {
        devShells.default = pkgs.mkShell {
          hardeningDisable = [ "all" ];
          packages = with pkgs; [
            stdenv.cc
            gfortran
            go
            jdk
            nodejs
            julia
            lua
            luajit
            pythonWithPackages
            pypy3
            rustc
            cargo
            mono
            dotnet-sdk
          ] ++ [
            stable.swift
            stable.swiftpm
            stable.swift-corelibs-libdispatch
          ];

          DOTNET_ROOT = "${pkgs.dotnet-sdk}/share/dotnet";
          LD_LIBRARY_PATH = "${stable.swift-corelibs-libdispatch}/lib";
        };
      });
}

