#ifndef _WEIGHT_H
#define _WEIGHT_H

double compute_point_weight(const double3& point_position, double wrad, uint atom, uint point);

void assign_cube_weights(LittleCube& cube);

#define WEIGHT_GPU 1

/* pesos en CPU */
#define BECKE_CUTOFF 0    // cota minima para el peso
#define BECKE 0           // formula de Becke o de Stratman
#define WEIGHT_CUTOFFS 0  // acota vecinos por nucleo

#endif