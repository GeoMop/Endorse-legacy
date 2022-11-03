The repository scale model. It describes water flow and transport of a tracer along a single storage horizontal well
with random large fractures, impact of the neighboring wells is considered. Boundary conditions are given by a piezometric head field that is reault of 
a larger model.

The result ot the repository model is the maximal concentration on the boundary as a function of the time.

Files:

- config.yaml: vstupni parametry geometrie
- mesh.py: tvorba geometrie a meshe; nazvy entit jsou zmenene oproti vzorovym yaml souborum!
- process.py: zakomentovane spusteni flowa, vysledkem je mesh
- 03_th.yaml: testovaci uloha, nazvy entit odpovidaji mesh.py
- setup_python_environment.sh: Shell script to setup a virtual environment.


Create virtual environment 'env' with preinstalled packages:

    python3 setup_venv.py
    

Run the calculation:

    env/bin/python process.py


Further notes:

number of tractable fractures: about 50
approx. time of the calculation:

