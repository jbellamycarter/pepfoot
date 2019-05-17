PepFoot
=======

PepFoot is intended for analysis and discovery in peptide footprinting, however 
it may be expanded to encompass more in future releases. PepFoot is released under the [LGPL-3.0 license](https://choosealicense.com/licenses/lgpl-3.0/)

PepFoot is currently maintained by [Jedd Bellamy-Carter, University of Nottingham](jeddidiah.bellamy-carter@nottingham.ac.uk). Any queries or improvements to the
software should be directed there or by submitting an issue on [GitHub](https://github.com/jbellamycarter/pepfoot/issues).

Installing
----------

Download the appropriate file from [Releases](https://github.com/jbellamycarter/pepfoot/releases) and follow the instructions below.

### Windows

To install PepFoot simply run the `pepFoot_1_1_WinOS.exe` installer and follow the wizard. The full PepFoot GUI should then run without problem.

### MacOSX

To install PepFoot simply mount the `pepFoot_1_1_MacOSX.dmg` file and drag `pepFoot.app` into your `Applications`. The full PepFoot GUI can then be accessed from this app.


### Linux/Python Users

It is recommended to run PepFoot through your local Python3 distribution for security. To install PepFoot simply extract this repository and run `python setup.py install`. This will add the command `pepfoot` to your local Python distribution as well as handle the package dependencies described below. Launching `pepfoot` from a terminal will launch the full PepFoot GUI.

Linux users can add `pepFoot.desktop` to your local `applications` directory and place a copy of `pepFoot.png` in your local `icons` directory. This `.desktop` file can now be used to launch the full PepFoot GUI.

####Requirements

* [Python 3](https://www.python.org)
* [PyQt >= 5.11](https://www.riverbankcomputing.com/software/pyqt/)
* [matplotlib](https://matplotlib.org/)
* [numpy](https://www.numpy.org/)
* [scipy](https://scipy.org/)
* [pyteomics](https://pyteomics.readthedocs.io/)
* [h5py](https://www.h5py.org/)

[NGL viewer](https://github.com/arose/ngl) is provided by the minified file `ngl.js` that is included in this directory.


Project Schema
--------------

| Key            	|                                                                                  Description 	|
|------------------	|----------------------------------------------------------------------------------	|
| name             	| Project name                                                                     	|
| creation date    	| Creation datetime of project                                                     	|
| data files       	| List of data files [filepath1, filepath2, ...]                                   	|
| sequence         	| Sequence of protein                                                              	|
| length range     	| Peptide length range parameters                                                  	|
| charge range     	| Charge state range parameters                                                    	|
| enzyme           	| Enzyme used for cleavage                                                         	|
| missed cleave    	| Number of missed cleavages                                                       	|
| fixed mods       	| List of fixed modifications applied to protein                                   	|
| differential mod 	| Differential modification applied to peptides                                    	|
| peptides         	| 1xN array of peptide ids e.g. 1-15                                               	|
| m/z array        	| 2xN array of m/z ranges for analysis [[unmod m/z],[mod m/z]]                     	|
| charge array     	| 1xN array of peptide charges                                                     	|
| rt array         	| 2xN array of rt ranges for analysis [[unmod rt],[mod rt]]                        	|
| areas            	| 2xNxM array of area values from analysis [[unmod area],[mod area], ...]          	|
| fractional mod   	| NxM array of fractional modification values from analysis                        	|
| treatment        	| Nested list with indices for data files grouped by treatment [[#, ...],[#, ...]] 	|
| pdb file         	| PDB file associated with project                                                 	|
