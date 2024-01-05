# HoloLattices
A generic code design for  to solve the problems typically arising in applications of holography to condensed matter systems. It consists of
a set of initialization blocks, which include importing the form of equations from Mathematica (or elsewhere), setting up the grid and choosing the initial solution
and the parameters of the problem; and the internal loop, where the solution is found after several iterations. 



[![DOI](https://zenodo.org/badge/504278109.svg)](https://zenodo.org/badge/latestdoi/504278109)


This code has successfully been used for publication purposes, such as for the paper [Isolated Zeros Destroy Fermi Surface in Holographic Models with a Lattice](https://link.springer.com/article/10.1007/JHEP01\(2020\)151)

# Setup Jupyter for the interactive examples
- Fistly, create a new virtual environment outside of the distribution tree with `virtualenv holo-env`

        .
        ├── holo-env
        └── HoloLattices
- Jump to root HoloLattices directory and activate it 

        source ../holo-env/bin/activate

- Then setup the package with `pip install .` in the root direcotry.

- The virtual environment is ready, (you will be able to deactivate it with `deactivate`)

- Now setup the custom IPython kernel for Jupyter `ipython kernel install --user --name=venv`

- This kernel will now be visible among choices in the top-right corner of Jupyter Notebook or VSCode (after restart)

See this thread for the reference (https://stackoverflow.com/questions/58119823/jupyter-notebooks-in-visual-studio-code-does-not-use-the-active-virtual-environm)

(for conda see https://janakiev.com/blog/jupyter-virtual-envs/#add-virtual-environment-to-jupyter-notebook):

# Tutorial
There are a set of tutorial codes available



# For developers
The repository includes Jupyter notebook. In order to track it correctly, make sure to include the local .gitconfig into your Git configuration with 

    git config --local include.path ../.gitconfig

See https://stackoverflow.com/questions/18329621/how-to-store-a-git-config-as-part-of-the-repository and https://git-scm.com/docs/git-config#_includes

Make sure to update jupyter in order to avoid errors in filter implementation.



