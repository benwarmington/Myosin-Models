**Revised DCR Model**

Revised DCR model is the single filament DCR model. Option for StS is equivalent to `c' within the body of work. Other options are the same. You will have to add the differential equations, GLMakie, Sundials, ProgressLogging and FFTW packages for it to run. 
The outputs of the program are UVALS, a 20000 by N+1 matrix giving the states of the rotors and backbone; vel, a matrix of the same size giving the angular velocities of rotors and the linear velocity of the filament; bound, a matrix of the same size with a boolean 
dictating wether a motor is in contact or not; phase_diff, a matrix showing the phase difference between neighboring rotors.  

**Sarc_Model**

This is the sarcomere model. Once again, variables at top are the same as the ones found in the body of work with the exception of StS == c. Outputs of this model are the same as the above, however the arrangement of the output matrix is slightly unweildy as
row 1 is the main filament/ z-disk displacement, then row 2 is the displacement at the basal end of the first filament, then there is the states of N rotors, then the next filament etc.

**Simulate Myosin System**

Clearly a lot of options here, I'd advise against touching h, fid and dev_time. The other options are as they appear in chapter 5 with the exception of F_back, the isotonic force and F_spring_const, which is kappa. The function returns x, the locations of the end of the lever arms; x_state_tracker, which shows which state (1-4) each motor is
in; F_elast, the elastic force each motor imparts on the filament; t, N, X0, F_spring_const; these are just the final time, number of motors, initial conditions and kappa respectively.
