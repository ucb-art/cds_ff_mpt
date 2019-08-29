# cds_ff_mpt
BAG primitives for cds_ff_mpt technology.

## Instructions for setting up BAG workspace using primitives repository
1. Make an empty directory.
2. Go into directory, initialize git:
    - git init
3. Add BAG primitives repo (cds_ff_mpt) as git submodule:
    - git submodule add \<URL>
4. Add BAG_framework as git submodule:
    - git submodule add git@github.com:ucb-art/BAG_framework.git
5. Run install script
    - ./\<primitives folder>/install.sh
6. Source .cshrc for tcsh, or .bashrc for bash
    - source .cshrc
    - source .bashrc
7. Run setup_submodules.py
    - ./setup_submodules.py
8. Commit changes

