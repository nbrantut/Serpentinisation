function [value,isterminal,direction] = grainsplit(c,cmax)

value = c-cmax;
isterminal = 1;
direction = 0;