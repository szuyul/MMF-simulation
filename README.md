# MMF-simulation
Compute the vectorial model of guided modes in an optical multimode fiber (MMF) and simulate fiber transmission in different representations.

Functions:
1. MMF_simTM_PIM computes the propagation invariant modes (PIMs) of a straight MMF with user defined fiber specification. It also provides a transmission matrix (TM) of the MMF in PIM representation with user defined fiber shape.
2. MMF_simTM_LP computes the linearly polarized (LP) modes of a straight MMF with user defined fiber specification. The LP modes are the approximated PIMs under the weakly guiding modes condition. It also provides a TM of the MMF in LP representation with user defined fiber shape.
3. MMF_simTM_camera computes the full TM of the MMF in experimental recording representation (e.g., camera pixel grid) with user defined fiber specification and shape. The full TM accounts for both horizontal and vertical polarization states.

Scripts (examples):
1. MMF_show_PIM displays the computed PIMs of a user defined MMF.
2. MMF_show_LP displays the computed LP modes as well as Laguerre-Gaussian modes of a user defined MMF.
3. MMF_focusing simulates the experiments of MMF TM measurements and focusing through the MMF with obtained TM.


created by Szu-Yu Lee (szuyul@mit.edu) 2017-2020
The Wellman Center for Photomedicine
