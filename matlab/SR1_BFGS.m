function [pk] = SR1_BFGS(gradxk, BK)
% Implementación del método de SR1 y BFGS para la resolución del pk a través 
% de la resolución de un sistema de ecuaciones
pk = linsolve(BK, -gradxk); 