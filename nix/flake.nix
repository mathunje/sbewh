{
    inputs = { 
        nixpkgs.url = "github:nixos/nixpkgs";
    };

    outputs = { self, nixpkgs }:
        let 
             pkgs = import nixpkgs {
                 system="x86_64-linux";
                 config.allowUnfree = true;
             };

             myPython = pkgs.python312.withPackages (p: with p; [
                numpy
                scipy
                matplotlib
                psutil
                ffmpeg-python
             ]);

        in {

            devShell.x86_64-linux = 
                pkgs.mkShell (pkgs.mkShell.override { stdenv = pkgs.cudaPackages.backendStdenv; }) {
                    buildInputs = [ 
                           pkgs.bash-completion
                           pkgs.mpi
                           pkgs.ffmpeg-full
                           pkgs.cudaPackages.cuda_cudart
                           pkgs.cudatoolkit
                           pkgs.cmake
                           pkgs.mkl
                           pkgs.zlib
                           pkgs.fftw
                           pkgs.gdb
                           myPython
                           pkgs.gtest
                           pkgs.boost
                    ];
                };
    };
}
