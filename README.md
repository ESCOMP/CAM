# ZM Kinetic Energy Parcel Calculation
## Summary
This branch includes the estimate of the parcel kinetic energy in the Zhang-McFarlane deep convection for the calculation of the entraining buoyancy profile which determines the intergrated CAPE. The convection top is no longer determined by the top most negative buyancy regions (depending on the value zmconv_num_cin (=1 in CAM6, =5 in CAM5). It is determined by the KE=0 level. KE is dtermined by the transfer of energy from the buoyancy PE, with a sepcific efficiency and initial parcel KE.
