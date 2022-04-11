# Solving Hidden Subset Sum Problem and Hidden Linear Combination Problem

This repository the known attacks for solving HSSP and HLCP.

## HSSP
To work with HSSP instance open sage and load the hssp.sage file

> load("hssp.sage")

For HSSP_n use:
>H=hssp(n) 
>H.gen_instance()

For HSSP_n^kappa use:
>H=hssp(n,kappa)
>H.gen_instance()
This generate the modulus H.x0, the weights vector H.a, the matrix H.x, the  sample vector H.b, kappa=-1 is HSSP_n by construction

