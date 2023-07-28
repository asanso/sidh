load('richelot_aux.sage')

def isogeny(base, degree, E, R,Q):    
    Ea = E
    Ra = R
    Qa =Q
    l = []
    while degree > 0:
        phi = Ea.isogeny(base^(degree-1)*Ra)
        Ea = phi.codomain()
        Ra = phi(Ra)
        Qa = phi(Qa)
        degree -= 1
        l.append(phi)
    return Ea,Qa,l


def apply_isogeny_list(isogenies, P, dual=False):
    PP = P
    isogenies = reversed(isogenies) if dual else isogenies
    for phi in isogenies:
        if dual:
            if phi.degree() == 1:
                continue
            phi = phi.dual()

        PP = phi(PP)

    return PP

p = 2^8*3^5 - 1
N = p+1
_.<I> = GF(p)[]
K.<i> = GF(p^2, modulus=I^2+1)
E = EllipticCurve(K, [1, 0])

cofactor = N/2^8
Q = cofactor*E.gens()[0]
phi =E.isogeny((N/3^2)*E.random_point())
phiQ = phi(Q)
E1 = phi.codomain()

cofactor2 = N/2^4
E0,Q0,chain= isogeny(2, 4, E, (2^4)*Q,Q)
R0 = cofactor2*E0.random_point()
while (2^4*R0).is_zero() and R0.weil_pairing(Q0, 2^4)^(2^3) == 1:
    R0 = cofactor2*E0.random_point()
Qd = apply_isogeny_list(chain,Q0,true)
Rd = apply_isogeny_list(chain,R0,true)
alpha = discrete_log(Rd,Qd,Qd.order(),operation='+')
R0 = R0 - alpha*Q0;
assert E0.isogeny(R0).codomain().j_invariant() == E.j_invariant()

E01,Q01,chain1= isogeny(2, 4, E1, (2^4)*phiQ,phiQ)

R01 = cofactor2*E01.random_point()
while (2^4*R01).is_zero() and R01.weil_pairing(Q01, 2^4)^(2^3) == 1:
    R01 = cofactor2*E01.random_point()
Qd1 = apply_isogeny_list(chain1,Q01,true)
Rd1 = apply_isogeny_list(chain1,R01,true)
alpha1 = discrete_log(Rd1,Qd1,Qd1.order(),operation='+')
R01 = R01 - alpha1*Q01;
assert E01.isogeny(R01).codomain().j_invariant() == E1.j_invariant()

f = 2^4 - 3^2
phi_f = E0.isogenies_prime_degree(f)[0].dual()
Ef= phi_f.domain()

w0 = Q0.weil_pairing(R0,2^4)
w1 = Q01.weil_pairing(R01,2^4)

lam = discrete_log(w1,w0,2^4,operation='*')
a = (lam^-1) % 2^4
b = (a*3^2) % 2^4

# Algorithm 1 step 1
c = f^-1 % 2^4

# Algorithm 1 step 2
Pb2 = c*phi_f.dual()(Q0)
Qb2 = c*phi_f.dual()(R0)

# Algorithm 1 step 3
chain2, codomain = Does22ChainSplit(Ef, E01, 3^2*Pb2,  3^2*Qb2, Q01,  b*R01, 4)

index = -1
if codomain[0].j_invariant() == E0.j_invariant():
    index = 0
else:
    index = 1
assert index != -1


# random 3^2 basis in E01
PA,QA =  E01.gens()
PA  = (2^8*3^3)*PA
QA  = (2^8*3^3)*2^8*QA

# Algorithm 1 step 4
# Map the points using the computed isogeny chain
imPAA = (Ef(0), PA)
for f in chain2:
    imPAA = f(imPAA)
# Algorithm 1 step 5
if not (3 * imPAA[index]).is_zero() and (3^2 * imPAA[index]).is_zero():
    assert imPAA[index].curve().j_invariant() == E0.j_invariant()
    assert imPAA[index].curve().isogeny(imPAA[index]).codomain().j_invariant() == E01.j_invariant()
    print("A generator of ker(phiA) is %",imPAA[index])

psi_list = imPAA[index].curve().isomorphisms(E0)
K = apply_isogeny_list(chain,psi_list[0](imPAA[index]),true)
assert E.isogeny(K).codomain().j_invariant() == E1.j_invariant()
