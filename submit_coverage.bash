#!/bin/bash
REPO_TOKEN='3dda8c89-c14f-47b3-aee8-2425e215a91f' julia -e 'cd(Pkg.dir("MyPkg")); using Coverage; Codecov.submit_token(process_folder())'

