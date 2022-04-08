load("hlcp.sage")

n=60
B=5
H=hlcp(n,B,-1)
H.gen_instance()
Ls=hlcp_attack(H)
assert Ls[-2]==n
H.gen_instance(m=n*2)
Lnso=hlcp_attack(H,'ns')
assert Lnso[-2]==n
H.gen_instance(m=n*4)
Lns=hlcp_attack(H,'ns')
assert Lns[-2]==n

print "\n\n--------------------------->HLCP_n Test: success!\n"
kappa=B*n/(B+1)
H=hssp(n,B,kappa)
H.gen_instance()
Lm=hssp_attack(H)
assert Lm[-2]==n
assert Lnso[-2]==n
Lns=hssp_attack(H,'ns')
assert Lns[-2]==n
print "\n\n--------------------------->HLCP_n^kappa Test: success!\n"
