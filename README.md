# Solving Hidden Subset Sum Problem and Hidden Linear Combination Problem

This repository the known attacks for solving HSSP and HLCP.

## HSSP
To work with HSSP instance open sage and load the hssp.sage file

> load("hssp.sage")

For HSSP_n use:
>H=hssp(n) 
>
>H.gen_instance()

For HSSP_n^kappa use:
>H=hssp(n,kappa)
>
>H.gen_instance()


This generate the modulus H.x0, the weights vector H.a, the matrix H.x, the  sample vector H.b; kappa=-1 is HSSP_n by construction.

To run the attacks use hssp_attack.    

>hssp_attack(H,alg='default')

H is the instance to be attacked and alg is the algorithm to use :

<ul>
<li>if alg='default' or alg='multi' runs the multivariate attack </li>
       
<li>if alg='ns_original' runs the original Nguyen-Stern attack </li>
       
<li>if alg='ns' runs the Nguyen-Stern attack with the improved orthogonal lattice attack </li>
       
<li>if alg='statistical' runs the heuristic statistical attack  </li>
 </ul>
