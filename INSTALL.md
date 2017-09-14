To use the GAIA repository, you need a few one-time steps to INSTALL this
in your shell environment.

1) get the repo, if you had not done that:

    git clone https://github.com/leemundy/GAIA/


2) set some useful things in your ~/.gitconfig file. We have a template
   in dot_gitconfig, so e.g.


    cp -i dot_gitconfig ~/.gitconfig

   you will need to edit a few things, at leeast your git name and email

3) we control your shell environment via the gaia_start.csh file. Adding a

   source ~/GAIA/gaia_start.csh

   should work for most people on the astrononmy desktops. For those using
   their own laptop environment, just make sure you have python3 (e.g. 
   miniconda3 or anaconda3 with jupyter, and ds9 in your environment)
